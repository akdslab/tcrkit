"""Peptide-MHC: an MHC allele to the protein sequence of its chain(s).

One function, :func:`mhc_sequence`. It looks an allele up in a reference of 186 class-I
and class-II alleles - 160 human, 24 mouse, plus ``B2M_human`` and ``B2M_mouse`` for the
class-I light chain::

    mhc_sequence('HLA-A*02:01')                 # the extracellular domain
    mhc_sequence('HLA-A*02:01', region='full')  # ... with its leader and tails
    mhc_sequence(['A*0201', 'DRB1*04:01'], paired=True)   # both chains of each

The allele is standardized with :func:`tcrkit.normalize_mhc` first, so any spelling that
understands works as a key.

The reference ships with the package and also backs :func:`tcrkit.annotate_mhc` and its
gene matching (:func:`tcrkit.mhc.best_mhc`). Its sequences are UniProt entries
(``data/mhc_uniprot.fasta``, exactly as UniProt serves them), one standing in for each
allele, and its annotations - membrane topology and the rest - are in
``data/mhc_annotations.csv``. :func:`load_mhc_table` joins the two into one DataFrame.

Sequences can come from another source: ``source='decodertcr'`` uses DecoderTCR's HLA
reference instead (downloaded once and cached; see :mod:`tcrkit._decodertcr`), and
``table=`` takes your own. Whatever the source, the regions mean the same thing:
DecoderTCR's chains get their topology from the matching bundled allele, carried over by
alignment, so ``region='extracellular'`` cuts the same stretch of the chain either way.

:func:`tcrkit.build_tcr_pmhc` is the whole-complex builder that uses this alongside the
TCR capabilities; this module is only the pMHC half.
"""

from __future__ import annotations

import ast
from functools import lru_cache
from pathlib import Path

import pandas as pd

from ._util import ERROR_COL, collect, squeeze as _squeeze
from .normalize import _parse_allele, normalize_mhc

__all__ = ['mhc_sequence', 'load_mhc_table', 'list_mhc', 'bundled_mhc_table_path',
           'MHC_SOURCES',
           'region_of', 'REGIONS', 'TOPOLOGY', 'TOPOLOGY_LETTERS',
           'topology_labels']


# ==========================================================================================
# MHC allele -> protein sequence (bundled reference)
# ==========================================================================================

_DATA = Path(__file__).resolve().parent / 'data'

#: Where MHC sequences can come from: ``'uniprot'`` (the bundled reference, default)
#: and ``'decodertcr'`` (DecoderTCR's HLA reference, fetched on first use).
MHC_SOURCES = ('uniprot', 'decodertcr')


def bundled_mhc_table_path() -> Path:
    """Path to the bundled MHC annotation table, ``mhc_annotations.csv``.

    Its sequences are in ``mhc_uniprot.fasta`` beside it; :func:`load_mhc_table` joins
    the two, and passing it this path does the same.
    """
    try:
        from importlib.resources import files
        return Path(str(files('tcrkit') / 'data' / 'mhc_annotations.csv'))
    except Exception:
        return _DATA / 'mhc_annotations.csv'


def _read_fasta(path: Path) -> dict:
    """``{accession: sequence}`` - ``sp|P61769|...`` -> ``P61769``, ``UPI...`` as is."""
    seqs, acc = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            first = line[1:].split()[0]
            acc = first.split('|')[1] if first.count('|') >= 2 else first
            seqs[acc] = []
        elif acc is not None:
            seqs[acc].append(line.strip())
    return {a: ''.join(parts) for a, parts in seqs.items()}


def _bundled_table() -> pd.DataFrame:
    """The bundled reference: ``mhc_annotations.csv`` joined to ``mhc_uniprot.fasta``."""
    table_path = bundled_mhc_table_path()
    fasta_path = table_path.with_name('mhc_uniprot.fasta')
    for path in (table_path, fasta_path):
        if not path.exists():
            raise FileNotFoundError(
                f'No MHC reference at {path}. It ships with tcrkit, so the install is '
                'incomplete - reinstall it. Or pass your own table with '
                'mhc_sequence(..., table=...).')
    table = pd.read_csv(table_path)
    seqs = _read_fasta(fasta_path)
    missing = sorted(set(table['accession']) - set(seqs))
    if missing:
        raise ValueError(f'{fasta_path.name} lacks {missing}, listed in {table_path.name}')
    table.insert(table.columns.get_loc('species') + 1, 'sequence',
                 table['accession'].map(seqs))
    table.insert(table.columns.get_loc('sequence') + 1, 'length',
                 table['sequence'].str.len())
    table['prediction'] = table['deeptmhmm_prediction']
    table['topology_from'] = 'deeptmhmm ' + table['deeptmhmm_version'].astype(str)

    def helix_ref(row):
        span = _spans(f"[{row['helix_ref_span']}]") if isinstance(row['helix_ref_span'], str) else []
        return row['sequence'][span[0][0] - 1:span[0][1]] if span else None

    table['helix_ref_sequence'] = table.apply(helix_ref, axis=1)
    table['in_detect'] = table['in_detect'].astype(bool)
    return table


@lru_cache(maxsize=None)
def load_mhc_table(path=None, source: str = 'uniprot') -> pd.DataFrame:
    """The MHC sequence reference as one DataFrame, cached for the session.

    The bundled one (``source='uniprot'``, the default) has 186 alleles - 160 human and
    24 mouse, class I and class II, plus ``B2M_human`` and ``B2M_mouse`` - keyed on the
    same compact allele name :func:`normalize_mhc` produces, so the two compose
    directly. Each allele's sequence is the UniProt entry standing in for it
    (``accession``), identical to the allele - alleles UniProt has no exact entry for
    are not bundled; :func:`list_mhc` shows what each source has.

    Columns: ``allele``; ``gene`` / ``mhc_class`` / ``species`` / ``mhc_prefix``;
    ``sequence`` (the full UniProt sequence, signal peptide and tails included) and its
    ``length``; ``source`` (``'uniprot'`` or ``'uniparc'``) and ``accession``; the
    membrane topology as 1-based span lists - ``signal_spans``, ``outside_spans``,
    ``tmhelix_spans``, ``inside_spans`` - with ``prediction`` / ``topology_from`` saying
    where they came from, and ``mature_start`` (the signal peptide's length);
    ``in_detect`` / ``detect_rank`` / ``detect_class``, the subset MHC genes are matched
    on (:func:`tcrkit.mhc.best_mhc`, behind :func:`tcrkit.annotate_mhc`);
    ``helix_ref_sequence`` / ``helix_ind``, the MHC Motif Atlas helix ranges for
    ``helix='tpm-bench'``; and the provenance columns of ``mhc_annotations.csv``.

    Args:
        path: read another table of the same shape (a CSV with ``allele`` and
            ``sequence`` columns) instead. The bundled table's own path
            (:func:`bundled_mhc_table_path`) gives the bundled reference.
        source: ``'uniprot'`` (default) or ``'decodertcr'`` - DecoderTCR's HLA
            reference, 87 human chains, downloaded on first use and cached
            (:func:`tcrkit._decodertcr.load_decodertcr_table`). Ignored with `path`.
    """
    if path is not None and Path(path).resolve() != bundled_mhc_table_path().resolve():
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f'No MHC sequence reference at {path}.')
        return pd.read_csv(path)
    if source not in MHC_SOURCES:
        raise ValueError(f'source must be one of {MHC_SOURCES}, got {source!r}')
    if source == 'decodertcr':
        from ._decodertcr import load_decodertcr_table
        return load_decodertcr_table()
    return _bundled_table()


def list_mhc(what: str = 'alleles', source: str = 'uniprot') -> pd.DataFrame:
    """What a source has: its MHC alleles, genes, or (DecoderTCR) pairs. For looking at.

    >>> list_mhc('genes').query("gene == 'DRB1'")[['gene', 'n_alleles']].values.tolist()
    [['DRB1', 30]]

    Args:
        what: ``'alleles'`` - one row per allele: ``allele``, ``gene``, ``mhc_class``,
            ``species``; ``'genes'`` - one row per gene, with ``n_alleles`` and the
            ``alleles``; ``'pairs'`` - the heterodimers DecoderTCR stores, one row each:
            ``pair`` (DecoderTCR's key), ``mhc_a`` / ``mhc_b`` (tcrkit's allele names),
            ``mhc_class``. The bundled reference stores single chains, so pairs need
            ``source='decodertcr'``.
        source: ``'uniprot'`` (default) or ``'decodertcr'`` (downloaded on first use).
    """
    table = load_mhc_table(source=source)
    if what == 'alleles':
        return table[['allele', 'gene', 'mhc_class', 'species']].reset_index(drop=True)
    if what == 'genes':
        genes = (table.groupby(['species', 'gene', 'mhc_class'], sort=True)['allele']
                 .agg(n_alleles='size', alleles=lambda a: ', '.join(a)).reset_index())
        return genes[['gene', 'mhc_class', 'species', 'n_alleles', 'alleles']]
    if what == 'pairs':
        if 'pairs' not in table:
            raise ValueError("the bundled reference stores single chains, not pairs; "
                             "pairs come from source='decodertcr'")
        allele_of = dict(zip(table['decodertcr_name'], table['allele']))
        keys = sorted({k for ks in table['pairs'] for k in ks.split(', ')})
        rows = []
        for key in keys:
            a, b = key.split('__')
            rows.append({'pair': key, 'mhc_a': allele_of[a], 'mhc_b': allele_of[b],
                         'mhc_class': 'I' if b == 'B2M' else 'II'})
        return pd.DataFrame(rows)
    raise ValueError(f"what must be 'alleles', 'genes' or 'pairs', got {what!r}")


#: The topology regions of a membrane protein, and the span column each is recorded in,
#: in order along the chain. The names are the standard ones: a **signal peptide** (the
#: cleaved N-terminal leader), the **extracellular** domain (the ectodomain, the part
#: outside the cell), the **transmembrane** helix, and the **cytoplasmic** tail.
TOPOLOGY = {
    'signal_peptide': 'signal_spans',
    'extracellular': 'outside_spans',
    'transmembrane': 'tmhelix_spans',
    'cytoplasmic': 'inside_spans',
}

#: What ``region=`` accepts: the four topology regions, plus ``'mature'`` (the chain
#: after its signal peptide is cleaved, which is what is actually expressed) and
#: ``'full'`` (the precursor as deposited). Common synonyms are accepted too.
REGIONS = (*TOPOLOGY, 'mature', 'full')

#: One letter per region in the ``topology_seq`` string, as DeepTMHMM writes them.
TOPOLOGY_LETTERS = {'signal_peptide': 'S', 'extracellular': 'O',
                    'transmembrane': 'M', 'cytoplasmic': 'I'}

_REGION_ALIASES = {
    'signal': 'signal_peptide', 'leader': 'signal_peptide',
    'ectodomain': 'extracellular', 'outside': 'extracellular',
    'tm': 'transmembrane', 'tmhelix': 'transmembrane',
    'intracellular': 'cytoplasmic', 'inside': 'cytoplasmic', 'tail': 'cytoplasmic',
    'precursor': 'full',
}


def _spans(value):
    """``'[(1, 286)]'`` -> ``[(1, 286)]``; anything unusable -> ``[]``."""
    if isinstance(value, str) and value.strip() not in ('', '[]'):
        try:
            return [tuple(s) for s in ast.literal_eval(value)]
        except (ValueError, SyntaxError):
            return []
    return list(value) if isinstance(value, (list, tuple)) else []


def topology_labels(row) -> list:
    """One region name per residue of a reference row's sequence, ``'-'`` where unknown."""
    seq = row.get('sequence')
    labels = ['-'] * (len(seq) if isinstance(seq, str) else 0)
    for name, col in TOPOLOGY.items():
        for start, end in _spans(row.get(col)):
            for i in range(max(0, start - 1), min(end, len(labels))):
                labels[i] = name
    return labels


def _runs(labels) -> dict:
    """``['extracellular', 'extracellular', 'transmembrane']`` ->
    ``{'extracellular': '1-2', 'transmembrane': '3'}``, 1-based."""
    out: dict[str, list[str]] = {}
    start = 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            if labels[start] != '-':
                out.setdefault(labels[start], []).append(
                    f'{start + 1}' + (f'-{i}' if i > start + 1 else ''))
            start = i
    return {k: ', '.join(v) for k, v in out.items()}


def _named(spans) -> str | None:
    """``{'extracellular': '1-286'}`` -> ``'extracellular'``; the regions present, in
    order along the chain."""
    present = [name for name in TOPOLOGY if spans.get(name)]
    return ', '.join(present) if present else None


def _region_index(row, region: str) -> list:
    """The 0-based positions of `region` in a reference row's sequence."""
    seq = row.get('sequence')
    n = len(seq) if isinstance(seq, str) else 0
    if region == 'full':
        return list(range(n))
    if region == 'mature':
        signal = _spans(row.get('signal_spans'))
        return list(range(max((e for _, e in signal), default=0), n))
    spans = _spans(row.get(TOPOLOGY[region]))
    if not spans:
        # no topology recorded: the whole chain is the best answer for the part that is
        # normally all of it, and no answer at all for the others
        return list(range(n)) if region == 'extracellular' else []
    return [i for start, end in spans for i in range(max(0, start - 1), min(end, n))]


def region_of(row, region: str = 'extracellular'):
    """One reference row's sequence cut down to `region`, and that cut's topology.

    Returns:
        ``(sequence, spans, topology_seq)`` - the residues of `region`; where each
        topology region falls *within them*, 1-based
        (``{'extracellular': '25-307', ...}``); and one letter per residue,
        :data:`TOPOLOGY_LETTERS`.
    """
    region = _REGION_ALIASES.get(str(region).lower(), region)
    if region not in REGIONS:
        raise ValueError(f'region must be one of {sorted(REGIONS)}, not {region!r}')
    seq = row.get('sequence')
    if not isinstance(seq, str):
        return None, {}, None
    index = _region_index(row, region)
    if not index:
        return None, {}, None
    all_labels = topology_labels(row)
    labels = [all_labels[i] for i in index]
    return (''.join(seq[i] for i in index), _runs(labels),
            ''.join(TOPOLOGY_LETTERS.get(x, '-') for x in labels))


def _sequence_map(table=None, region: str = 'extracellular', source: str = 'uniprot') -> dict:
    """``{allele -> sequence}`` for a reference, bundled or your own.

    Values are ``(sequence, topology)``. Your own table needs only ``allele`` and
    ``sequence`` columns - with topology spans too if you want ``region=`` to cut
    anything - and its allele spellings are
    normalized for you, so ``'HLA-A*02:01'`` in the file is found by ``'A*0201'`` in a
    query. Both the spelling as written and the normalized one are keys, so neither way
    round misses.
    """
    if table is None:
        tab, normalize_keys = load_mhc_table(source=source), False
    else:
        tab = pd.read_csv(table) if isinstance(table, (str, Path)) else table
        normalize_keys = True

    missing = {'allele', 'sequence'} - set(tab.columns)
    if missing:
        raise ValueError(f'An MHC table needs {sorted(missing)} column(s); '
                         f'got {list(tab.columns)}.')

    seqs = {}                            # allele -> (sequence, spans, topology_seq)
    for _, r in tab.iterrows():
        seqs[str(r['allele']).strip()] = region_of(r, region)
    if normalize_keys:
        written = list(seqs)
        for raw, key in zip(written, normalize_mhc(written, errors='coerce')['allele']):
            if isinstance(key, str):
                seqs.setdefault(key, seqs[raw])
    return seqs


def mhc_sequence(alleles, paired: bool = False, region: str = 'extracellular',
                 mark_regions: bool = False, col_mapping=None, table=None,
                 squeeze: bool = False, source: str | None = None,
                 **kwargs) -> pd.DataFrame:
    """Look up the protein sequence of MHC/HLA alleles.

    >>> out = mhc_sequence('HLA-A*02:01')
    >>> out[['allele', 'mhc_class']].to_dict('records')
    [{'allele': 'A0201', 'mhc_class': 'I'}]
    >>> out['sequence'][0][:20]
    'GSHSMRYFFTSVSRPGRGEP'

    The allele is normalized with :func:`normalize_mhc` first, so any spelling that
    function understands works as a key - ``'HLA-A*02:01'``, ``'A*0201'`` and
    ``'A0201'`` all find the same row.

    Args:
        alleles: a single allele name or any collection of them - a list, Series, dict,
            list of records, or a DataFrame.
        paired: look up *both* chains of the heterodimer. Class I gets its heavy chain
            plus the invariant B2M light chain; class II gets the alpha/beta pair
            (DRB* pairing with DRA0101). Returns ``mhc_a_*`` / ``mhc_b_*`` columns.
        region: which part of the chain to return (see :data:`REGIONS`):
            ``'extracellular'`` (the default - past the signal peptide, before the
            membrane, i.e. the part outside the cell), ``'full'`` (the chain as
            deposited), ``'mature'`` (everything after the signal peptide),
            ``'signal'``, ``'transmembrane'`` or ``'cytoplasmic'``. An allele with no
            topology annotation falls back to its whole sequence; the ``topology``
            column says which were annotated.
        mark_regions: add ``topology_seq``, a string as long as ``sequence`` with one
            letter per residue saying which region it is in - ``S`` signal peptide,
            ``O`` extracellular, ``M`` transmembrane, ``I`` cytoplasmic, ``-``
            unknown (:data:`TOPOLOGY_LETTERS`, DeepTMHMM's alphabet).
        col_mapping: where the allele lives when the input is a frame or records, e.g.
            ``{'allele': 'mhc'}``.
        table: your own reference instead of the bundled one - a DataFrame or a CSV
            path with ``allele`` and ``sequence`` columns. Allele spellings in it are
            normalized, so they need not already be in the compact form.
        squeeze: return the plain result instead of a frame.
        source: where the sequences come from when `table` is not given:
            ``'uniprot'`` (default, the bundled reference) or ``'decodertcr'``
            (DecoderTCR's HLA reference: downloaded on first use, then cached). Regions
            mean the same thing for both - see :data:`MHC_SOURCES`.
        **kwargs: forwarded to :func:`normalize_mhc` (e.g. ``restrict_allele_fields``).

    Returns:
        A DataFrame indexed like the input. ``allele_input`` echoes what you passed,
        then the normalized ``allele`` and its ``sequence`` and ``length``, plus
        ``gene`` / ``mhc_class`` / ``species``. With ``paired=True`` the allele and
        sequence columns are doubled as ``mhc_a_allele`` / ``mhc_a_sequence`` /
        ``mhc_b_allele`` / ``mhc_b_sequence``.

        ``regions`` names the topology regions the returned sequence actually
        contains, in order - ``'extracellular'`` for the default,
        ``'signal_peptide, extracellular, transmembrane, cytoplasmic'`` for
        ``region='full'``.

        One column per region (:data:`TOPOLOGY`) then says where that region falls
        *inside the sequence returned*, 1-based - ``signal_peptide`` ``'1-24'``,
        ``extracellular`` ``'25-307'`` - and is empty for a region the returned
        sequence does not contain. So ``region='extracellular'`` gives
        ``extracellular '1-286'`` with the other three empty.

        An allele that does not parse, or parses but is not in the reference, leaves
        the sequence empty and says why in ``error`` - the reference is not all of
        IMGT, so misses are expected. Filter with ``out.sequence.notna()``.
    """
    if table is not None and source is not None:
        raise ValueError('give table= or source=, not both')
    seqs = _sequence_map(table, region=region, source=source or 'uniprot')

    inp = collect(('allele',), (alleles,), col_mapping=col_mapping)
    # normalize_mhc does the parsing; errors='coerce' so one bad allele cannot abort
    parsed = normalize_mhc(list(inp['allele']), paired=paired, errors='coerce', **kwargs)

    def look_up(key, raw=None):
        """The sequence for a normalized allele, falling back to the raw input.

        The fallback matters for the species-tagged class-I light chains: the reference
        keys them as ``B2M_human`` / ``B2M_mouse``, which is the form
        :func:`normalize_mhc` *emits* for a pair, but normalizing one of those labels on
        its own goes through mhcgnomes and comes back as a bare ``B2M`` / ``H2-B2M``
        with the species dropped. Trying the raw string too makes
        ``mhc_sequence('B2M_human')`` work, which is how these appear in real pMHC
        tables.
        """
        for candidate in (key, raw):
            if isinstance(candidate, str) and candidate.strip() in seqs:
                return (*seqs[candidate.strip()], None)
        if key is None or (isinstance(key, float) and pd.isna(key)):
            return None, {}, None, None
        return None, {}, None, f'{key!r} not in the reference'

    recs = []
    for (_, p) in parsed.iterrows():
        err = p.get(ERROR_COL)
        if isinstance(err, str):
            recs.append({ERROR_COL: err})
            continue
        if paired:
            a, b = p.get('mhc_a_allele'), p.get('mhc_b_allele')
            a_seq, a_spans, a_marks, a_err = look_up(a, p.get('allele_input'))
            b_seq, b_spans, b_marks, b_err = look_up(b)
            recs.append({
                'mhc_a_allele': a, 'mhc_a_sequence': a_seq,
                'mhc_b_allele': b, 'mhc_b_sequence': b_seq,
                **{k: p.get(k) for k in ('gene', 'mhc_class', 'species')},
                **({'mhc_a_topology_seq': a_marks, 'mhc_b_topology_seq': b_marks}
                   if mark_regions else {}),
                'regions': _named(a_spans),
                ERROR_COL: ' | '.join(e for e in (a_err, b_err) if e) or None,
            })
        else:
            key = p.get('allele')
            seq, spans, marks, miss = look_up(key, p.get('allele_input'))
            recs.append({
                'allele': key, 'sequence': seq,
                'length': len(seq) if isinstance(seq, str) else None,
                **{k: p.get(k) for k in ('gene', 'mhc_class', 'species')},
                'regions': _named(spans),
                **({'topology_seq': marks} if mark_regions else {}),
                **{name: spans.get(name) for name in TOPOLOGY},
                ERROR_COL: miss,
            })

    cols = (['mhc_a_allele', 'mhc_a_sequence', 'mhc_b_allele', 'mhc_b_sequence']
            if paired else ['allele', 'sequence', 'length'])
    if mark_regions:
        cols += (['mhc_a_topology_seq', 'mhc_b_topology_seq'] if paired
                 else ['topology_seq'])
    cols += ['gene', 'mhc_class', 'species', 'regions']
    if not paired:
        cols += list(TOPOLOGY)
    cols += [ERROR_COL]
    out = pd.DataFrame(recs, index=inp.index).reindex(columns=cols)
    out.insert(0, 'allele_input', list(inp['allele']))
    return _squeeze(out, ('allele_input',)) if squeeze else out
