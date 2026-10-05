"""Peptide-MHC: an MHC allele to the protein sequence of its chain(s).

One function, :func:`mhc_sequence`. It looks an allele up in a reference of 251 class-I
and class-II alleles - 226 human, 25 mouse, plus ``B2M_human`` and ``B2M_mouse`` for the
class-I light chain::

    mhc_sequence('HLA-A*02:01')                 # the extracellular domain
    mhc_sequence('HLA-A*02:01', full=True)      # ... plus the TM and cytoplasmic tails
    mhc_sequence(['A*0201', 'DRB1*04:01'], paired=True)   # both chains of each

The allele is standardized with :func:`tcrkit.normalize_mhc` first, so any spelling that
understands works as a key.

The reference is derived from IMGT/HLA and is **not kept in version control**, so it is
absent from a fresh clone or install. :func:`load_mhc_table` raises
:class:`FileNotFoundError` saying so; put ``mhc_sequences.csv`` at
:func:`bundled_mhc_table_path` or pass your own with ``table=``. Nothing else in tcrkit
depends on it.

:func:`tcrkit.build_tcr_pmhc` is the whole-complex builder that uses this alongside the
TCR capabilities; this module is only the pMHC half.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from ._util import ERROR_COL, collect, squeeze as _squeeze
from .normalize import _parse_allele, normalize_mhc

__all__ = ['mhc_sequence', 'load_mhc_table', 'bundled_mhc_table_path']


# ==========================================================================================
# MHC allele -> protein sequence (bundled reference)
# ==========================================================================================

def bundled_mhc_table_path() -> Path:
    """Path to the MHC sequence reference shipped with the package."""
    try:
        from importlib.resources import files
        return Path(str(files('tcrkit') / 'data' / 'mhc_sequences.csv'))
    except Exception:
        return Path(__file__).resolve().parent / 'data' / 'mhc_sequences.csv'


@lru_cache(maxsize=None)
def load_mhc_table(path=None) -> pd.DataFrame:
    """Read the MHC sequence reference and cache it for the session.

    251 alleles - 226 human and 25 mouse, class I and class II, plus ``B2M_human`` and
    ``B2M_mouse`` - keyed on the same compact allele name :func:`normalize_mhc` produces,
    so the two compose directly.

    Columns: ``allele``, ``sequence`` (the extracellular domain actually used),
    ``full_sequence`` (including the transmembrane and cytoplasmic tails), ``len`` /
    ``len_full``, and the descriptive ``gene`` / ``mhc_class`` / ``species`` /
    ``mhc_prefix`` / ``source``.
    """
    path = Path(path) if path else bundled_mhc_table_path()
    if not path.exists():
        raise FileNotFoundError(
            f'No MHC sequence reference at {path}. It is derived from IMGT/HLA and is '
            'not distributed with the source, so a fresh clone or install does not have '
            'it: copy mhc_sequences.csv to that path, or pass your own table with '
            "mhc_sequence(..., table=load_mhc_table('/path/to/your.csv')). Every other "
            'capability works without it.')
    return pd.read_csv(path)


def _sequence_map(table=None, full: bool = False) -> dict:
    """``{allele -> sequence}`` for a reference, bundled or your own.

    Your own table needs only ``allele`` and ``sequence`` columns (plus
    ``full_sequence`` if you want ``full=True``), and its allele spellings are
    normalized for you, so ``'HLA-A*02:01'`` in the file is found by ``'A*0201'`` in a
    query. Both the spelling as written and the normalized one are keys, so neither way
    round misses.
    """
    if table is None:
        tab, normalize_keys = load_mhc_table(), False
    else:
        tab = pd.read_csv(table) if isinstance(table, (str, Path)) else table
        normalize_keys = True

    missing = {'allele', 'sequence'} - set(tab.columns)
    if missing:
        raise ValueError(f'An MHC table needs {sorted(missing)} column(s); '
                         f'got {list(tab.columns)}.')

    col = 'full_sequence' if full and 'full_sequence' in tab.columns else 'sequence'
    seqs = {str(a).strip(): seq for a, seq in zip(tab['allele'], tab[col])}
    if normalize_keys:
        written = list(seqs)
        for raw, key in zip(written, normalize_mhc(written, errors='coerce')['allele']):
            if isinstance(key, str):
                seqs.setdefault(key, seqs[raw])
    return seqs


def mhc_sequence(alleles, paired: bool = False, full: bool = False, col_mapping=None,
                 table=None, squeeze: bool = False, **kwargs) -> pd.DataFrame:
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
        full: return ``full_sequence`` (with the transmembrane and cytoplasmic tails)
            instead of the extracellular domain.
        col_mapping: where the allele lives when the input is a frame or records, e.g.
            ``{'allele': 'mhc'}``.
        table: your own reference instead of the bundled one - a DataFrame or a CSV
            path with ``allele`` and ``sequence`` columns. Allele spellings in it are
            normalized, so they need not already be in the compact form.
        squeeze: return the plain result instead of a frame.
        **kwargs: forwarded to :func:`normalize_mhc` (e.g. ``restrict_allele_fields``).

    Returns:
        A DataFrame indexed like the input. ``allele_input`` echoes what you passed,
        then the normalized ``allele`` and its ``sequence`` and ``length``, plus
        ``gene`` / ``mhc_class`` / ``species``. With ``paired=True`` the allele and
        sequence columns are doubled as ``mhc_a_allele`` / ``mhc_a_sequence`` /
        ``mhc_b_allele`` / ``mhc_b_sequence``.

        An allele that does not parse, or parses but is not in the reference, leaves
        the sequence empty and says why in ``error`` - the reference covers 251
        alleles, not all of IMGT, so misses are expected. Filter with
        ``out.sequence.notna()``.
    """
    seqs = _sequence_map(table, full=full)

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
                return seqs[candidate.strip()], None
        if key is None or (isinstance(key, float) and pd.isna(key)):
            return None, None
        return None, f'{key!r} not in the reference'

    recs = []
    for (_, p) in parsed.iterrows():
        err = p.get(ERROR_COL)
        if isinstance(err, str):
            recs.append({ERROR_COL: err})
            continue
        if paired:
            a, b = p.get('mhc_a_allele'), p.get('mhc_b_allele')
            a_seq, a_err = look_up(a, p.get('allele_input'))
            b_seq, b_err = look_up(b)
            recs.append({
                'mhc_a_allele': a, 'mhc_a_sequence': a_seq,
                'mhc_b_allele': b, 'mhc_b_sequence': b_seq,
                **{k: p.get(k) for k in ('gene', 'mhc_class', 'species')},
                ERROR_COL: ' | '.join(e for e in (a_err, b_err) if e) or None,
            })
        else:
            key = p.get('allele')
            seq, miss = look_up(key, p.get('allele_input'))
            recs.append({
                'allele': key, 'sequence': seq,
                'length': len(seq) if isinstance(seq, str) else None,
                **{k: p.get(k) for k in ('gene', 'mhc_class', 'species')},
                ERROR_COL: miss,
            })

    cols = (['mhc_a_allele', 'mhc_a_sequence', 'mhc_b_allele', 'mhc_b_sequence']
            if paired else ['allele', 'sequence', 'length'])
    cols += ['gene', 'mhc_class', 'species', ERROR_COL]
    out = pd.DataFrame(recs, index=inp.index).reindex(columns=cols)
    out.insert(0, 'allele_input', list(inp['allele']))
    return _squeeze(out, ('allele_input',)) if squeeze else out
