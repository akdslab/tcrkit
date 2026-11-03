"""Membrane topology of the bundled MHC sequences, by DeepTMHMM.

Internal maintenance script, not part of the API::

    python -m tcrkit._build.annotate_mhc_deeptmhmm                   # run DeepTMHMM, update the CSV
    python -m tcrkit._build.annotate_mhc_deeptmhmm --raw-dir out/    # also keep its output files
    python -m tcrkit._build.annotate_mhc_deeptmhmm --parse out/      # reuse a finished run
    python -m tcrkit._build.annotate_mhc_deeptmhmm --missing-only    # only entries without topology

Runs DeepTMHMM (Hallgren et al. 2022) through BioLib on every sequence in
``data/mhc_uniprot.fasta``, pinned to version :data:`VERSION` - the authors' recommended
release; later BioLib versions are a wrapper around a separate service. Needs ``pybiolib``
(``pip install pybiolib``) and network access; no account.

Each residue gets one of DeepTMHMM's labels - ``S`` signal peptide, ``O`` outside
(extracellular), ``M`` transmembrane helix, ``I`` inside (cytoplasmic), ``B`` beta
strand, ``P`` periplasm - and the runs are written to ``data/mhc_annotations.csv``, by
accession, as 1-based inclusive span lists: ``signal_spans``,
``outside_spans``, ``tmhelix_spans``, ``inside_spans``. Also ``deeptmhmm_prediction``
(``SP+TM``, ``TM``, ``SP``, ``GLOB``), ``mature_start`` (the signal peptide's length;
the mature chain, and IMGT position 1 of the first G domain, starts at
``mature_start + 1``) and ``deeptmhmm_version``.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / 'data'
TABLE = DATA / 'mhc_annotations.csv'
FASTAS = (DATA / 'mhc_uniprot.fasta',)
VERSION = '1.0.24'
APP = f'DTU/DeepTMHMM:{VERSION}'
RESULT = 'predicted_topologies.3line'

#: DeepTMHMM label -> span column, the names tcrkit's MHC table uses.
LABELS = {'S': 'signal_spans', 'O': 'outside_spans', 'M': 'tmhelix_spans',
          'I': 'inside_spans'}


def accession(header: str) -> str:
    """The accession in a FASTA header: ``sp|P61769|...`` -> ``P61769``,
    ``prf||1006281B ...`` -> ``1006281B``."""
    first = header.lstrip('>').split()[0]
    if first.count('|') >= 2 and first.split('|')[0] in ('sp', 'tr'):
        return first.split('|')[1]
    return first.rstrip('|').split('|')[-1]


def runs(topology: str, label: str) -> list[tuple[int, int]]:
    """1-based inclusive ranges where `topology` has `label`."""
    return [(m.start() + 1, m.end()) for m in re.finditer(f'{re.escape(label)}+', topology)]


def parse_3line(path: Path) -> pd.DataFrame:
    """DeepTMHMM's ``predicted_topologies.3line`` -> one row per sequence."""
    lines = [l.rstrip('\n') for l in path.read_text().splitlines() if l.strip()]
    rows = []
    for i in range(0, len(lines), 3):
        header, seq, topo = lines[i:i + 3]
        if len(seq) != len(topo):
            raise ValueError(f'{header}: sequence and topology lengths differ')
        name, _, prediction = header[1:].partition(' | ')
        row = {'accession': accession(name), 'deeptmhmm_prediction': prediction.strip(),
               'topology_length': len(seq)}
        for label, col in LABELS.items():
            spans = runs(topo, label)
            row[col] = str(spans) if spans else None
        sig = runs(topo, 'S')
        row['mature_start'] = sig[0][1] if sig and sig[0][0] == 1 else 0
        rows.append(row)
    return pd.DataFrame(rows)


def run_deeptmhmm(fasta: Path, out_dir: Path) -> Path:
    """Run DeepTMHMM on `fasta` via BioLib; save its output files in `out_dir`."""
    import biolib

    app = biolib.load(APP)
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copy(fasta, Path(tmp) / 'input.fasta')
        cwd = os.getcwd()
        os.chdir(tmp)                     # BioLib uploads input files by relative path
        try:
            job = app.cli(args=['--fasta', 'input.fasta'])
        finally:
            os.chdir(cwd)
    if job.get_exit_code() != 0:
        raise SystemExit(f'DeepTMHMM failed:\n{job.get_stderr().decode()[-2000:]}')
    out_dir.mkdir(parents=True, exist_ok=True)
    job.save_files(str(out_dir))
    found = list(out_dir.rglob(RESULT))
    if not found:
        raise SystemExit(f'no {RESULT} in the DeepTMHMM output ({out_dir})')
    return found[0]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--raw-dir', type=Path,
                   help='keep DeepTMHMM output files here (default: a temporary dir)')
    p.add_argument('--parse', type=Path, metavar='DIR',
                   help=f'skip the run; parse {RESULT} from a previous --raw-dir')
    p.add_argument('--missing-only', action='store_true',
                   help='run DeepTMHMM only on entries that have no topology yet')
    args = p.parse_args(argv)

    table = pd.read_csv(TABLE)
    wanted = None
    if args.missing_only:
        no_topo = table['deeptmhmm_prediction'].isna() if 'deeptmhmm_prediction' in table \
            else pd.Series(True, index=table.index)
        wanted = set(table.loc[no_topo, 'accession'])
        if not wanted:
            print('every entry already has its topology')
            return 0

    if args.parse:
        result = next(args.parse.rglob(RESULT))
        topo = parse_3line(result)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = args.raw_dir or Path(tmp)
            joined = Path(tmp) / 'all.fasta'
            records = [r for f in FASTAS if f.exists()
                       for r in re.split(r'(?m)^(?=>)', f.read_text()) if r.strip()]
            if wanted is not None:
                records = [r for r in records if accession(r.split('\n', 1)[0]) in wanted]
            joined.write_text(''.join(records))
            result = run_deeptmhmm(joined, out_dir)
            topo = parse_3line(result)

    fasta_ids = wanted if wanted is not None else {
        accession(l) for f in FASTAS if f.exists()
        for l in f.read_text().splitlines() if l.startswith('>')}
    missing = fasta_ids - set(topo['accession'])
    if missing:
        raise SystemExit(f'DeepTMHMM returned no topology for {sorted(missing)}')
    cols = [c for c in topo.columns if c != 'accession']
    for col in cols:
        if col not in table:
            table[col] = None
        table[col] = table[col].astype(object)
    rows = table['accession'].isin(topo['accession'])
    new = topo.set_index('accession')
    for col in cols:
        table.loc[rows, col] = table.loc[rows, 'accession'].map(new[col])
    table.loc[rows, 'deeptmhmm_version'] = VERSION
    table.to_csv(TABLE, index=False)
    print(f'{len(topo)} sequences annotated: '
          f"{topo['deeptmhmm_prediction'].value_counts().to_dict()}; wrote {TABLE.name}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
