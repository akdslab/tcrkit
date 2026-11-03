"""Download the bundled MHC sequences from UniProt: ``data/mhc_uniprot.fasta``.

Internal maintenance script, not part of the API::

    python -m tcrkit._build.fetch_mhc_sequences            # refresh from UniProt
    python -m tcrkit._build.fetch_mhc_sequences --check    # compare, write nothing

Which entries: the ``source`` / ``accession`` columns of ``data/mhc_annotations.csv``,
one entry per tcrkit MHC allele (several alleles can share one). ``source`` is
``uniprot`` for a UniProtKB entry, or ``uniparc`` for the few alleles only UniProt's
sequence archive holds (UniParc, e.g. a PRF record). All are written to one FASTA,
exactly as UniProt serves them. ``source_release`` records the UniProt release and
``sequence_version`` each UniProtKB entry's ``SV=``, so a refresh shows what changed.

UniProtKB and UniParc data are CC BY 4.0 (https://www.uniprot.org/help/license).

After a refresh, re-run ``python -m tcrkit._build.annotate_mhc_deeptmhmm``: the membrane
topology belongs to the sequences.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / 'data'
TABLE = DATA / 'mhc_annotations.csv'
FASTA = DATA / 'mhc_uniprot.fasta'
SOURCES = ('uniprot', 'uniparc')
UNIPROT = 'https://rest.uniprot.org/uniprotkb/accessions'
UNIPARC = 'https://rest.uniprot.org/uniparc/{}.fasta'
BATCH = 100


def accession(header: str) -> str:
    """The accession in a FASTA header: ``sp|P61769|...`` -> ``P61769``,
    ``UPI00005FF33A status=active`` -> ``UPI00005FF33A``."""
    first = header.lstrip('>').split()[0]
    if first.count('|') >= 2 and first.split('|')[0] in ('sp', 'tr'):
        return first.split('|')[1]
    return first.rstrip('|').split('|')[-1]


def parse_fasta(text: str) -> list[tuple[str, str, str]]:
    """``[(accession, header, record text), ...]`` in file order."""
    out = []
    for rec in re.split(r'(?m)^(?=>)', text):
        if rec.strip():
            header = rec.split('\n', 1)[0]
            out.append((accession(header), header, rec if rec.endswith('\n') else rec + '\n'))
    return out


def _get(url: str, params: dict) -> tuple[str, dict]:
    with urllib.request.urlopen(f'{url}?{urllib.parse.urlencode(params)}', timeout=120) as r:
        return r.read().decode(), dict(r.headers)


def download(source: str, accessions: list[str]) -> tuple[str, str]:
    """The source's FASTA for `accessions`, in that order, and its release."""
    texts, release = {}, None
    for i in range(0, len(accessions), BATCH):
        chunk = accessions[i:i + BATCH]
        if source == 'uniprot':
            text, headers = _get(UNIPROT, {'accessions': ','.join(chunk), 'format': 'fasta'})
            release = headers.get('X-UniProt-Release', release)
        else:                                    # UniParc: one request per entry
            text = ''
            for upi in chunk:
                part, headers = _get(UNIPARC.format(upi), {})
                text += part
                release = headers.get('X-UniProt-Release', release)
        for acc, _, rec in parse_fasta(text):
            texts[acc] = rec
    missing = [a for a in accessions if a not in texts]
    if missing:
        raise SystemExit(f'{source} returned no entry for {missing} (obsolete or merged?). '
                         'Pick replacements in mhc_annotations.csv first.')
    return ''.join(texts[a] for a in accessions), release


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--check', action='store_true',
                   help='report differences from the bundled FASTAs, write nothing')
    args = p.parse_args(argv)

    table = pd.read_csv(TABLE)
    unknown = set(table['source']) - set(SOURCES)
    if unknown:
        raise SystemExit(f'unknown source(s) {sorted(unknown)} in {TABLE.name}')
    text, records, release = '', [], {}
    for source in SOURCES:
        accessions = list(dict.fromkeys(table.loc[table['source'] == source, 'accession']))
        if accessions:
            part, release[source] = download(source, accessions)
            text += part
    records = parse_fasta(text)
    if FASTA.exists():
        old = {acc: rec for acc, _, rec in parse_fasta(FASTA.read_text())}
        changed = [a for a, _, rec in records if old.get(a) != rec]
        print(f'UniProt {release.get("uniprot")}: {len(records)} entries, {len(changed)} differ '
              f'from {FASTA.name}{": " + ", ".join(changed) if changed else ""}')
    else:
        print(f'UniProt {release.get("uniprot")}: {len(records)} entries')
    if args.check:
        return 0

    FASTA.write_text(text)
    table['source_release'] = table['source'].map(release)
    version = {acc: int(m.group(1)) if (m := re.search(r'SV=(\d+)', h)) else None
               for acc, h, _ in records}
    table['sequence_version'] = table['accession'].map(version).astype('Int64')
    table.to_csv(TABLE, index=False)
    print(f'wrote {FASTA.name} and updated {TABLE.name}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
