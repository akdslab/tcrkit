"""Restore the IMGT junction residues on a non-canonical CDR3.

A CDR3 is reported either *with* its IMGT junction residues (104 ``C`` ... 118
``F``/``W``/``C``/``L``/``V``) or as 105-117 only. This module uses a reference table of
the first kind - which junction each V/J gene pair actually uses - to put the flanks back
on CDR3s of the second kind::

    add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1')     # -> 'CASSYSGNTEAFF'
    add_cdr3_junctions('AVMDSNYQLI', 'TRAJ33*01')    # -> 'CAVMDSNYQLIW'
    junction_of('TRAJ35')                        # -> ('C/C', 0.9999, 14633, 'J')

Only the CDR3 and the J gene are required. The V gene and species just make the match
more specific, and species is detected from the gene when not given. Gene names are
looked up both as passed and as tidytcells standardizes them ('TRAV23/6' ->
'TRAV23/DV6'), and allele suffixes are stripped, so 'TRBJ1-1*01' works too.

A table is bundled with the package and read once per process, so nothing needs passing
around. Both loci live in the one table. :func:`build_junction_table` learns a new one from your
own reference frame.

Table columns (one row per gene key):
    level      which key the row is for, most to least specific: 'V+J' = species + V
               gene + J gene, 'J' = species + J gene (v_gene empty), 'J/any-species' =
               J gene alone (species 'any'). A lookup walks them in that order and stops
               at the first hit.
    species    'human', 'mouse', ... or 'any' for the species-blind level
    v_gene     V gene, empty when the level ignores it
    j_gene     J gene
    junction   the junction that group uses most often, written 'start/end' - e.g. 'C/F'
    purity     its share of the group's rows; 1.0 means the group only ever uses this one
    n          number of reference rows behind the group
    locus      'TRA' / 'TRB', taken from the J gene
    alphabet   every junction seen for that locus, e.g. 'C/F,C/V'. This is what tells
               :func:`add_cdr3_junctions` which residues count as flanks, so a CDR3 that
               already carries its junction is returned untouched instead of being given
               a second one.

"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from ._util import ERROR_COL, collect, drop_unset, is_empty, rows as _rows
from ._util import species_or_default, squeeze as _squeeze

__all__ = ['add_cdr3_junctions', 'junction_of', 'load_junction_table',
           'build_junction_table', 'bundled_junction_table_path', 'MIN_PCT']

#: Junctions rarer than this share of a locus are garbled CDR3s, not real ones.
MIN_PCT = 0.01

_LEVELS = (('V+J', True), ('J', False), ('J/any-species', False))


def bundled_junction_table_path() -> Path:
    """Path to the junction lookup table shipped with the package."""
    try:
        from importlib.resources import files
        return Path(str(files('tcrkit') / 'data' / 'cdr3_junction_lookup.csv'))
    except Exception:
        return Path(__file__).resolve().parent / 'data' / 'cdr3_junction_lookup.csv'


@lru_cache(maxsize=None)
def _std_gene(symbol, species='any'):
    """(gene, species) for a V/J symbol via tidytcells: 'TRAV23/6' -> ('TRAV23/DV6', 'human').

    (None, None) when tidytcells cannot place the symbol. species='any' detects it.
    """
    if not isinstance(symbol, str) or not symbol.strip():
        return None, None
    # The single-value core, not the public DataFrame-returning normalize_tcr: this is a
    # per-symbol lookup on a hot path and wants the scalar shape.
    from .normalize import _standardize

    _result, extras = _standardize(symbol, species, return_species=True,
                                   return_gene=True)
    # extras['species'] is already tcrkit's name ('human' / 'mouse'), which is how the
    # lookup table is keyed - no translation needed.
    return extras.get('gene'), extras.get('species')


@lru_cache(maxsize=None)
def load_junction_table(path=None):
    """Read the junction lookup table and cache it for the session.

    Defaults to the table bundled with the package; pass `path` to use your own (e.g.
    one produced by :func:`build_junction_table`).
    """
    path = Path(path) if path else bundled_junction_table_path()
    if not path.exists():
        raise FileNotFoundError(
            f"No junction lookup table at {path}. Pass path=, or build one with "
            f"tcrkit.build_junction_table()."
        )
    return pd.read_csv(path)


def _index(tab):
    """{(level, species, v_gene, j_gene): (junction, purity, n)} + per-locus flanks, cached."""
    if 'index' not in tab.attrs:
        def blank(x):                                # 'any V' survives a CSV round-trip
            return '' if pd.isna(x) else x
        tab.attrs['index'] = {(r.level, r.species, blank(r.v_gene), r.j_gene):
                              (r.junction, r.purity, r.n) for r in tab.itertuples()}
        f = {lo: (frozenset(x[0] for x in a.split(',')),
                  frozenset(x[-1] for x in a.split(',')))
             for lo, a in zip(tab.locus, tab.alphabet)}
        # unknown locus: accept any flank the table knows about
        f[''] = tuple(frozenset().union(*x) for x in zip(*f.values()))
        tab.attrs['flanks'] = f
    return tab.attrs['index']


def _spellings(symbol):
    """((gene as written, standardized gene), detected species) for a V/J call or gene."""
    if not isinstance(symbol, str) or not symbol.strip():
        return (), None
    raw = symbol.split('*')[0].strip()                   # TRAJ26*01 -> TRAJ26
    std, sp = _std_gene(raw)
    return tuple(dict.fromkeys([raw] + ([std] if std else []))), sp


def _junction_lookup(j, v=None, species=None, table=None):
    """The junction a V/J gene pair uses: ``(junction, purity, n, level)``.

    The single-value core behind :func:`junction_of`.

    Tries species + V + J, then species + J, then the J gene alone, and stops at the
    first hit. Returns ``(None, nan, 0, None)`` when the genes are unknown to the table.

    Args:
        j: the J gene or allele call. Required - it is what determines the 118 residue.
        v: the V gene or allele call. Optional; only makes the match more specific.
        species: 'human' / 'mouse'. Detected from the gene names when not given.
        table: a lookup table from :func:`load_junction_table`; the bundled one by default.
    """
    tab = load_junction_table() if table is None else table
    idx = _index(tab)
    js, sp = _spellings(j)
    vs = _spellings(v)[0] or ('',)
    # is_empty, not falsiness: a missing value arriving as NaN is *truthy*, and would
    # otherwise be stringified into the species key as 'nan' and match nothing. Free
    # text is tolerated here because it comes from the caller's own data.
    species = species_or_default(species, default=None)
    species = str(species or sp or 'any').lower()
    for level, use_v in _LEVELS:
        key_sp = 'any' if level.endswith('any-species') else species
        for jg in js:
            for vg in (vs if use_v else ('',)):
                hit = idx.get((level, key_sp, vg, jg))
                if hit:
                    return (*hit, level)
    return (None, float('nan'), 0, None)


def _complete_one(cdr3, j, v=None, species=None, table=None):
    """Put the IMGT junction residues back on one CDR3. The single-value core.

    A CDR3 that already carries both flanks is returned untouched - no lookup needed.
    One that starts on the junction residue keeps it and only gains its 118. One that
    does not is taken to have lost *both* flanks, so a trailing F counts as part of the
    body.

    Args:
        cdr3: the CDR3 amino-acid sequence.
        j: the J gene or allele call. Required.
        v: the V gene or allele call. Optional; only makes the match more specific.
        species: 'human' / 'mouse'. Detected from the gene names when not given.
        table: a lookup table from :func:`load_junction_table`; the bundled one by default.

    Returns:
        The completed CDR3, or ``None`` when a junction is needed but the genes are
        unknown to the table.
    """
    if not isinstance(cdr3, str) or len(cdr3) < 2:
        return None
    tab = load_junction_table() if table is None else table
    _index(tab)
    flanks = tab.attrs['flanks']
    starts, ends = flanks.get(str(j)[:3].upper() if isinstance(j, str) else '', flanks[''])
    if cdr3[0] in starts and cdr3[-1] in ends:
        return cdr3
    jn = _junction_lookup(j, v, species, tab)[0]
    if jn is None:
        return None
    return cdr3 + jn[-1] if cdr3[0] in starts else jn[0] + cdr3 + jn[-1]




def junction_of(j, v=None, species=None, col_mapping=None, table=None,
                squeeze: bool = False) -> pd.DataFrame:
    """The junction a V/J gene pair uses. Returns one row per gene pair.

    >>> junction_of('TRBJ1-1')[['j', 'junction', 'purity']]
             j junction  purity
    0  TRBJ1-1      C/F     1.0

    Tries species + V + J, then species + J, then the J gene alone, and stops at the
    first hit.

    Args:
        j: a single J gene / allele call, or any collection of them. Required - it is
            what determines the 118 residue.
        v: the V gene / allele call(s). Optional; only makes the match more specific.
        species: 'human' / 'mouse'. Detected from the gene names when not given.
        col_mapping: where the fields live when the input is a frame or records.
        table: a lookup table from :func:`load_junction_table`; the bundled one by default.
        squeeze: return the plain result instead of a frame - the junction string for
            one gene pair, a Series of them for several.

    Returns:
        A DataFrame indexed like the input with the genes echoed back, then
        ``junction`` ('C/F', ...; None when the genes are unknown), ``purity`` (the
        share of reference rows agreeing), ``n`` (how many rows back it) and ``level``
        (which lookup level hit: 'V+J', 'J' or 'J/any-species').
    """
    tab = load_junction_table() if table is None else table
    inp = collect(('j', 'v', 'species'), (j, v, species), col_mapping=col_mapping,
                  required=('j',))

    recs = [dict(zip(('junction', 'purity', 'n', 'level'),
                     _junction_lookup(r['j'], r.get('v'), r.get('species'), tab)))
            for r in _rows(inp)]
    out = pd.DataFrame(recs, index=inp.index,
                       columns=['junction', 'purity', 'n', 'level'])
    res = drop_unset(pd.concat([inp, out], axis=1), ('v', 'species'))
    return _squeeze(res[['junction']], ()) if squeeze else res


def add_cdr3_junctions(cdr3, j=None, v=None, species=None, col_mapping=None, table=None,
                 keep_unknown=False, squeeze: bool = False) -> pd.DataFrame:
    """Restore the IMGT junction residues on CDR3s reported as 105-117 only.

    >>> add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1')[['cdr3', 'cdr3_fixed', 'junction']]
              cdr3     cdr3_fixed junction
    0  ASSYSGNTEAF  CASSYSGNTEAFF      C/F

    >>> list(add_cdr3_junctions(['ASSYSGNTEAF', 'AVMDSNYQLI'],
    ...                         ['TRBJ1-1', 'TRAJ33'])['cdr3_fixed'])
    ['CASSYSGNTEAFF', 'CAVMDSNYQLIW']

    A CDR3 that already carries both flanks is returned untouched - no lookup needed.
    One that starts on the junction residue keeps it and only gains its 118. One that
    does not is taken to have lost *both* flanks, so a trailing F counts as part of the
    body.

    Args:
        cdr3: a single CDR3 amino-acid sequence, or any collection of them - a list,
            Series, dict, list of records, or a DataFrame. When a container carries
            everything, `j` / `v` / `species` are read from it too.
        j: the J gene / allele call(s). Required - the J gene is what determines the
            118 residue, and it is not always F (TRBJ1-1 ends in F, TRAJ33 in W,
            TRAJ35 in C).
        v: the V gene / allele call(s). Optional; only makes the match more specific.
            Omitting it applies the broader J-level lookup, which is usually what you
            want for repair.
        species: 'human' / 'mouse'. Detected from the gene names when not given.
        col_mapping: where the fields live when the input is a frame or records, e.g.
            ``{'cdr3': 'cdr3_aa', 'j': 'j_call'}``.
        table: a lookup table from :func:`load_junction_table`; the bundled one by default.
        keep_unknown: where the J gene is unknown to the table, fall back to the
            original CDR3 instead of leaving ``cdr3_fixed`` empty.
        squeeze: return the plain result instead of a frame - the completed CDR3 for one
            input, a Series of them for several.

    Returns:
        A DataFrame indexed like the input, so it can be assigned straight back onto
        your own frame::

            df['cdr3_fixed'] = add_cdr3_junctions(
                df, col_mapping={'cdr3': 'cdr3_aa', 'j': 'j_call'})['cdr3_fixed']

        Columns: the ``cdr3`` / ``j`` / ``v`` / ``species`` you gave, then
        ``cdr3_fixed`` (None where the genes are unknown, unless ``keep_unknown``),
        ``changed`` (whether anything was added), the ``junction`` / ``purity`` / ``n`` /
        ``level`` the lookup used, and ``error`` - None when the lookup hit, else why
        it did not.
    """
    tab = load_junction_table() if table is None else table
    inp = collect(('cdr3', 'j', 'v', 'species'), (cdr3, j, v, species),
                  col_mapping=col_mapping, required=('cdr3', 'j'))

    recs = []
    for r in _rows(inp):
        cd, jj, vv, sp = r['cdr3'], r['j'], r.get('v'), r.get('species')
        fixed = _complete_one(cd, jj, vv, sp, tab)
        junction, purity, n, level = _junction_lookup(jj, vv, sp, tab)
        if fixed is None and keep_unknown:
            fixed = cd
        recs.append({'cdr3_fixed': fixed,
                     'changed': fixed is not None and fixed != cd,
                     'junction': junction, 'purity': purity, 'n': n, 'level': level,
                     ERROR_COL: (None if junction is not None else
                                 f'no junction known for J gene {jj!r}')})

    out = pd.DataFrame(recs, index=inp.index,
                       columns=['cdr3_fixed', 'changed', 'junction', 'purity', 'n',
                                'level', ERROR_COL])
    res = drop_unset(pd.concat([inp, out], axis=1), ('v', 'species'))
    return _squeeze(res[['cdr3_fixed']], ()) if squeeze else res

# ---------------------------------------------------------------------------
# Building a table from your own reference data
# ---------------------------------------------------------------------------

def _chain_table(d):
    """The lookup rows for one already-renamed frame: every gene key, at all levels."""
    d = d.dropna(subset=['cdr3', 'j_gene', 'species'])
    d = d[d.cdr3.str.len() > 4]
    d = d.assign(junction=d.cdr3.str[0] + '/' + d.cdr3.str[-1])

    locus = d.j_gene.str[:3]
    alphabet = {}                                   # locus -> the junctions it really uses
    for lo, g in d.groupby(locus):
        share = g.junction.value_counts(normalize=True) * 100
        alphabet[lo] = ','.join(sorted(share[share >= MIN_PCT].index))
    d = d[[jn in alphabet[lo] for jn, lo in zip(d.junction, locus)]]

    def group(x, level):                            # the junction each group uses most
        keys = ['species', 'v_gene', 'j_gene']
        n = x.value_counts(keys + ['junction']).rename('n').reset_index()
        n['purity'] = n.n / n.groupby(keys).n.transform('sum')
        return n.drop_duplicates(keys).assign(level=level)

    tab = pd.concat([group(d.dropna(subset=['v_gene']), 'V+J'),
                     group(d.assign(v_gene=''), 'J'),
                     group(d.assign(v_gene='', species='any'), 'J/any-species')],
                    ignore_index=True)

    def std(g):                                     # also file genes under tidytcells'
        return (_std_gene(g)[0] or g) if g else g   # spelling

    tab = pd.concat([tab, tab.assign(v_gene=[std(g) for g in tab.v_gene],
                                     j_gene=[std(g) for g in tab.j_gene])],
                    ignore_index=True)
    tab = tab.drop_duplicates(['level', 'species', 'v_gene', 'j_gene'])
    tab = tab.reset_index(drop=True)
    tab['locus'] = tab.j_gene.str[:3]
    tab['alphabet'] = tab.locus.map(alphabet)
    return tab[['level', 'species', 'v_gene', 'j_gene', 'junction', 'purity', 'n',
                'locus', 'alphabet']]


def build_junction_table(ref, cdr3_col, j_col, v_col=None, species_col=None, path=None):
    """Learn a junction lookup table from your own reference data.

    The reference must report CDR3s *with* their junction residues - that is what the
    table learns from. Rows whose CDR3 carries a junction seen in under
    :data:`MIN_PCT` of the locus are dropped as garbled.

    >>> ref = pd.DataFrame({'cdr3': ['CASSYSGNTEAFF'] * 5, 'j': ['TRBJ1-1'] * 5,
    ...                     'sp': ['human'] * 5})
    >>> tab = build_junction_table(ref, 'cdr3', 'j', species_col='sp')
    >>> add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1', table=tab)
    'CASSYSGNTEAFF'

    Args:
        ref: a DataFrame, or a list of them (concatenated), of reference rearrangements.
        cdr3_col: column holding CDR3s *with* their junction residues.
        j_col: column holding the J calls.
        v_col: column holding the V calls. Optional, but without it the 'V+J' level of
            the table is empty and only the J-level lookups will hit.
        species_col: column holding the species. Optional; 'any' is assumed when absent,
            which still populates the species-blind level.
        path: also write the table to this CSV.

    Returns:
        The lookup table, ready to pass as ``table=`` to :func:`add_cdr3_junctions`.
    """
    frames = [ref] if isinstance(ref, pd.DataFrame) else list(ref)
    parts = []
    for f in frames:
        for name, col in (('cdr3_col', cdr3_col), ('j_col', j_col),
                          ('v_col', v_col), ('species_col', species_col)):
            if col is not None and col not in f.columns:
                raise ValueError(f"{name}={col!r} is not a column in the reference frame.")
        d = pd.DataFrame({
            'cdr3': f[cdr3_col].astype('string'),
            'j_gene': f[j_col].astype('string').str.split('*').str[0],
            'v_gene': (f[v_col].astype('string').str.split('*').str[0]
                       if v_col else pd.Series(pd.NA, index=f.index, dtype='string')),
            'species': (f[species_col].astype('string').str.lower()
                        if species_col else pd.Series('any', index=f.index,
                                                      dtype='string')),
        })
        parts.append(d)

    d = pd.concat(parts, ignore_index=True).drop_duplicates()
    tab = _chain_table(d)
    if path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tab.to_csv(path, index=False)
        print(f"junction lookup: {len(d):,} reference rows -> {len(tab):,} keys -> {path}")
    return tab
