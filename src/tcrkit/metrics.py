"""Repertoire metrics: diversity, overlap, generation probability and TCR distances.

    diversity                 richness, Shannon entropy, clonality and friends, via
                              `scikit-bio <https://scikit.bio>`_
    morisita_horn             overlap between two repertoires, or all pairs of them
    generation_probability    how likely V(D)J recombination was to produce a CDR3, via
                              `OLGA <https://github.com/statbiophys/OLGA>`_
    tcrdist                   pairwise TCR distances, via
                              `tcrdist3 <https://tcrdist3.readthedocs.io>`_

Like the rest of tcrkit these are wrappers: the arithmetic belongs to the established
package, and tcrkit's job is the input handling, the naming and the frame that comes
back. The exception is :func:`morisita_horn`, which is computed here because no
maintained Python package provides it - see its docstring.

These packages are heavy (olga pulls in numba, tcrdist3 and scikit-bio pull in scipy and
friends), so they are an optional extra::

    pip install 'tcrkit[metrics]'      # or: uv sync --extra metrics

They are imported only when called, so the rest of tcrkit works without them, and a
missing package produces an instruction rather than a traceback.

Counts vs clones
----------------
:func:`diversity` and :func:`morisita_horn` work on clone abundances. You can hand them
either already-counted abundances or the raw clone labels, and they will tell the
difference: a column of numbers is taken as counts, a column of strings as one row per
cell, counted for you.
"""

from __future__ import annotations

import math
import warnings
from contextlib import contextmanager
from functools import lru_cache

import pandas as pd

from ._util import ERROR_COL, collect, drop_unset, mute_stdout, species_of
from ._util import squeeze as _squeeze

__all__ = ['diversity', 'morisita_horn', 'generation_probability', 'tcrdist']

#: The 20 amino acids OLGA will accept in a CDR3.
_AA = set('ACDEFGHIKLMNPQRSTVWY')

#: OLGA ships a model per species and locus; these are the T-cell ones.
_OLGA_MODELS = {
    ('human', 'alpha'): 'human_T_alpha', ('human', 'beta'): 'human_T_beta',
    ('mouse', 'alpha'): 'mouse_T_alpha', ('mouse', 'beta'): 'mouse_T_beta',
}


@contextmanager
def _quiet_tcrdist():
    """Silence two warnings tcrdist3 emits that say nothing about your data.

    * ``db_file must be 'alphabeta_gammadelta_db.tsv' or ...`` - tcrdist3's own default
      is ``combo_xcr_2024-03-05.tsv``, which its validator was never updated to accept,
      so the check fires on every stock call. Naming an older db to quieten it would
      trade gene coverage for cosmetics, so the warning is dropped instead.
    * ``IProgress not found`` - tqdm looking for ipywidgets when imported under Jupyter.
      Install ipywidgets and it goes away for real; this only stops it appearing in the
      middle of your results.

    Anything else tcrdist3 warns about still comes through.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='db_file must be')
        warnings.filterwarnings('ignore', message='IProgress not found')
        yield


def _missing(package: str, extra: str = 'metrics'):
    return ImportError(
        f"{package} is needed for this metric but is not installed. "
        f"Install the optional extra:  pip install 'tcrkit[{extra}]'   "
        f"(or: uv sync --extra {extra})"
    )


# ==========================================================================================
# Abundances
# ==========================================================================================

def _abundances(data, count=None, by=None, col_mapping=None):
    """Any input shape -> ``{repertoire: pandas Series of clone -> count}``.

    Accepts already-counted abundances (a numeric column) or raw clone labels (a string
    column, counted here), optionally split into repertoires by `by`.
    """
    col_mapping = col_mapping or {}

    if isinstance(data, pd.DataFrame):
        clone_col = col_mapping.get('clone', 'clone')
        count_col = col_mapping.get('count', count or 'count')
        by_col = col_mapping.get(by, by)
        if by_col is not None and by_col not in data.columns:
            raise ValueError(f"by={by_col!r} is not a column in the DataFrame.")
        if clone_col not in data.columns:
            # no clone column: the first non-count, non-grouping column identifies clones
            others = [c for c in data.columns if c not in (count_col, by_col)]
            if not others:
                raise ValueError(
                    "Could not find a clone column. Pass "
                    "col_mapping={'clone': 'your_column'}."
                )
            clone_col = others[0]
        has_counts = count_col in data.columns

        def counted(part):
            if has_counts:
                return part.groupby(clone_col)[count_col].sum()
            return part[clone_col].value_counts()

        if by_col is None:
            return {'all': counted(data)}
        return {name: counted(part) for name, part in data.groupby(by_col)}

    # a bare Series / list / dict
    if isinstance(data, dict):
        return {'all': pd.Series(data)}
    s = data if isinstance(data, pd.Series) else pd.Series(list(data))
    if pd.api.types.is_numeric_dtype(s):
        return {'all': s}                      # already abundances
    return {'all': s.value_counts()}           # clone labels, one row per cell


def diversity(data, count=None, by=None, col_mapping=None,
              squeeze: bool = False) -> pd.DataFrame:
    """Diversity of one or more repertoires. One row per repertoire.

    >>> out = diversity(['A', 'A', 'A', 'B', 'C'])
    >>> out[['n_clones', 'n_cells']].to_dict('records')
    [{'n_clones': 3, 'n_cells': 5}]
    >>> round(float(out['clonality'][0]), 3)
    0.135

    Args:
        data: clone abundances, or the raw clone labels to count. A list/Series of
            numbers is taken as abundances; a list/Series of strings as one entry per
            cell. A DataFrame is grouped by `by` and counted per clone.
        count: name of the count column when `data` is a DataFrame and the column is
            not called ``count``. Absent, rows are counted.
        by: column to split repertoires on - one output row per value.
        col_mapping: ``{'clone': ..., 'count': ...}`` when your columns are named
            differently.
        squeeze: return the metrics without the repertoire column.

    Returns:
        A DataFrame with one row per repertoire:

        ====================  ================================================
        column                scikit-bio metric
        ====================  ================================================
        n_clones              ``sobs`` - distinct clones (richness)
        n_cells               total abundance
        shannon               ``shannon(base=2)`` - entropy in bits
        pielou_evenness       ``pielou_e`` - shannon / log2(n_clones), 0-1
        clonality             1 - pielou_evenness; 0 = even, 1 = monoclonal
        simpson               ``simpson`` - 1 - sum(p^2), the Gini-Simpson index
        dominance             ``dominance`` - sum(p^2), the complement of simpson
        top_clone_fraction    frequency of the most abundant clone
        ====================  ================================================

        ``shannon`` is in **bits** (base 2), which is what makes ``pielou_evenness``
        and so ``clonality`` fall on 0-1 and stay comparable between samples of
        different richness - raw entropy does not, since it rises just from having more
        clones.

        Note scikit-bio's naming: ``simpson`` is *1 - sum(p^2)* (the probability two
        draws differ, often written Gini-Simpson) and ``dominance`` is *sum(p^2)*. Some
        immunology papers use "Simpson" for the latter, so check which one you mean.

    Requires the ``metrics`` extra (``pip install 'tcrkit[metrics]'``).
    """
    try:
        from skbio.diversity import alpha
    except ImportError as e:                                # pragma: no cover
        raise _missing('scikit-bio') from e

    recs = []
    for name, counts in _abundances(data, count=count, by=by,
                                    col_mapping=col_mapping).items():
        c = pd.Series(counts).astype(float)
        c = c[c > 0]
        values = c.values
        n = int(len(values))
        total = float(c.sum())

        # base=2 so entropy is in bits and pielou_e (which divides by log2 of the
        # richness) is the matching evenness on 0-1
        shannon = float(alpha.shannon(values, base=2)) if n else 0.0
        evenness = float(alpha.pielou_e(values)) if n > 1 else 0.0

        recs.append({
            'repertoire': name,
            'n_clones': int(alpha.sobs(values)) if n else 0,
            'n_cells': total,
            'shannon': shannon,
            'pielou_evenness': evenness,
            'clonality': 1.0 - evenness,
            'simpson': float(alpha.simpson(values)) if n else 0.0,
            'dominance': float(alpha.dominance(values)) if n else 0.0,
            'top_clone_fraction': float((c / total).max()) if n and total else 0.0,
        })

    out = pd.DataFrame(recs)
    return _squeeze(out, ('repertoire',)) if squeeze else out


def morisita_horn(data, other=None, count=None, by=None, col_mapping=None,
                  long: bool = False, squeeze: bool = False) -> pd.DataFrame:
    """Morisita-Horn overlap: 0 when two repertoires share no clones, 1 when identical.

    >>> morisita_horn({'A': 10, 'B': 5}, {'A': 10, 'B': 5}, squeeze=True)
    1.0
    >>> morisita_horn({'A': 10, 'B': 5}, {'C': 10, 'D': 5}, squeeze=True)
    0.0

    Args:
        data: one repertoire (with `other` given), or a table of several to compare
            pairwise (with `by` naming the repertoire column).
        other: the second repertoire, when comparing exactly two.
        count, by, col_mapping: as :func:`diversity`.
        long: return one row per pair (``a``, ``b``, ``morisita_horn``) instead of a
            square matrix - easier to filter and join, and far smaller for many
            repertoires. Same switch as :func:`tcrdist`.
        squeeze: return the plain number instead of a frame - a float for one pair, a
            Series for several.

    Returns:
        A square symmetric DataFrame indexed and columned by repertoire, with 1.0 on the
        diagonal; or the long form with ``long=True``. Comparing exactly two gives a
        2x2 matrix, so ``squeeze=True`` is the way to get the single number.

    The index is abundance-weighted and insensitive to sampling depth, which is why it
    is preferred over Jaccard for repertoires sequenced to different depths: it compares
    the *shape* of the two clone-size distributions, not just which clones are present.

    Note:
        This is the one metric here tcrkit computes itself rather than wrapping. No
        maintained Python package provides Morisita-Horn - it is absent from both
        scikit-bio (which has jaccard and braycurtis) and pyrepseq (jaccard, overlap,
        renyi entropy); the usual reference implementation is R's
        ``vegan::vegdist(method = "horn")``, against which this agrees to 4e-16.
    """
    def index(x, y):
        clones = x.index.union(y.index)
        a = x.reindex(clones).fillna(0).astype(float)
        b = y.reindex(clones).fillna(0).astype(float)
        na, nb = a.sum(), b.sum()
        if not na or not nb:
            return float('nan')
        pa, pb = a / na, b / nb
        denom = float((pa ** 2).sum() + (pb ** 2).sum())
        return 0.0 if denom == 0 else float(2 * (pa * pb).sum() / denom)

    if other is not None:
        reps = {'a': next(iter(_abundances(data, count=count,
                                           col_mapping=col_mapping).values())),
                'b': next(iter(_abundances(other, count=count,
                                           col_mapping=col_mapping).values()))}
    else:
        reps = _abundances(data, count=count, by=by, col_mapping=col_mapping)
        if by is None and len(reps) == 1:
            raise ValueError(
                "Give a second repertoire as `other=`, or pass by= to name the column "
                "that splits `data` into repertoires."
            )

    names = list(reps)
    matrix = pd.DataFrame([[index(reps[a], reps[b]) for b in names] for a in names],
                          index=names, columns=names)

    pairs = matrix.stack().rename('morisita_horn').reset_index()
    pairs.columns = ['a', 'b', 'morisita_horn']
    pairs = pairs[pairs.a != pairs.b].reset_index(drop=True)

    if squeeze:
        return _squeeze(pairs[['morisita_horn']].iloc[:1] if len(names) == 2
                        else pairs[['morisita_horn']], ())
    return pairs if long else matrix


# ==========================================================================================
# Generation probability (OLGA)
# ==========================================================================================

@lru_cache(maxsize=None)
def _olga_model(species: str, chain: str):
    """OLGA's generation-probability model for a species/chain, loaded once."""
    try:
        import olga.generation_probability as _pgen
        import olga.load_model as _load
    except ImportError as e:                                # pragma: no cover
        raise _missing('olga') from e
    import os

    import olga

    key = (species_of(species), chain.lower())
    if key not in _OLGA_MODELS:
        raise ValueError(
            f"No OLGA model for species={species!r} chain={chain!r}. "
            f"Available: {sorted(_OLGA_MODELS)}"
        )
    d = os.path.join(os.path.dirname(olga.__file__), 'default_models', _OLGA_MODELS[key])
    params, marginals = (os.path.join(d, 'model_params.txt'),
                         os.path.join(d, 'model_marginals.txt'))
    v_anchors = os.path.join(d, 'V_gene_CDR3_anchors.csv')
    j_anchors = os.path.join(d, 'J_gene_CDR3_anchors.csv')

    # beta/heavy are VDJ (they have a D gene); alpha/light are VJ
    vdj = key[1] in ('beta', 'delta', 'heavy')
    genomic = _load.GenomicDataVDJ() if vdj else _load.GenomicDataVJ()
    genomic.load_igor_genomic_data(params, v_anchors, j_anchors)
    generative = _load.GenerativeModelVDJ() if vdj else _load.GenerativeModelVJ()
    generative.load_and_process_igor_model(marginals)
    return (_pgen.GenerationProbabilityVDJ(generative, genomic) if vdj
            else _pgen.GenerationProbabilityVJ(generative, genomic))


def generation_probability(cdr3, v=None, j=None, chain='beta', species='human',
                           col_mapping=None, use_genes: bool = True,
                           squeeze: bool = False) -> pd.DataFrame:
    """How likely V(D)J recombination was to produce each CDR3, via OLGA.

    >>> out = generation_probability('CASSIRSSYEQYF')
    >>> f"{out['pgen'][0]:.2e}"
    '1.80e-07'

    A low pgen means the CDR3 is hard to generate by chance, so seeing it in two people
    is more surprising - this is the usual way to weight public-ness.

    Args:
        cdr3: a CDR3 amino-acid sequence (**with** its conserved C and F), or any
            collection of them. OLGA wants the junction, which is the form
            :func:`tcrkit.add_cdr3_junctions` produces.
        v, j: optional V and J calls. Each one you give narrows the recombinations OLGA
            sums over, giving a lower and more specific pgen; they are independent, so
            V alone or J alone is fine. Without either, OLGA marginalises over all V/J.
        chain: 'alpha' or 'beta'.
        species: ``'human'`` or ``'mouse'`` (spellings like 'Mus musculus' understood).
        col_mapping: where the fields live when the input is a frame or records.
        use_genes: set False to ignore `v` / `j` even when the input carries them.
        squeeze: return the plain result instead of a frame.

    Returns:
        A DataFrame indexed like the input, echoing ``cdr3`` / ``v`` / ``j`` and adding
        ``pgen`` and ``log10_pgen``. A CDR3 OLGA cannot parse gets pgen 0.0 and the
        reason in ``error``; nothing raises, and nothing is printed - note that a
        pgen of 0.0 with no error means a real CDR3 this model cannot produce, which is
        different from one it could not read.

    Requires the ``metrics`` extra (``pip install 'tcrkit[metrics]'``).
    """
    model = _olga_model(species, chain)
    inp = collect(('cdr3', 'v', 'j'), (cdr3, v, j), col_mapping=col_mapping,
                  required=('cdr3',))

    recs = []
    for row in inp.to_dict('records'):
        # OLGA takes V and J independently - either, both or neither narrow the
        # recombinations it sums over, so pass whichever were given.
        gene_args = [None, None]
        if use_genes:
            for pos, key in enumerate(('v', 'j')):
                val = row.get(key)
                if isinstance(val, str) and val.strip():
                    gene_args[pos] = val.strip()
        seq = row['cdr3']
        if not isinstance(seq, str) or not seq or not set(seq.upper()) <= _AA:
            # OLGA prints a complaint and returns 0.0 for these, which is
            # indistinguishable from a genuinely unproducible CDR3 - so catch them
            # first and say so in the column instead.
            p, err = 0.0, 'not an amino-acid sequence'
        else:
            try:
                # OLGA writes directly to stdout on some inputs; keep it out of the way
                # print_warnings=False is OLGA's own switch; mute_stdout catches the
                # paths that write regardless
                with mute_stdout():
                    p = float(model.compute_aa_CDR3_pgen(seq, *gene_args,
                                                         print_warnings=False))
                err = None
            except Exception as e:
                p, err = 0.0, f'{type(e).__name__}: {e}'
        recs.append({'pgen': p,
                     'log10_pgen': math.log10(p) if p > 0 else float('-inf'),
                     ERROR_COL: err})

    out = pd.concat([inp, pd.DataFrame(recs, index=inp.index)], axis=1)
    out = drop_unset(out, ('v', 'j'))      # no point echoing genes that were not given
    return _squeeze(out, inp.columns) if squeeze else out


# ==========================================================================================
# TCR distances (tcrdist3)
# ==========================================================================================

def tcrdist(df, chains=('beta',), col_mapping=None, species='human',
            long: bool = False, **kwargs):
    """Pairwise TCR distances, via tcrdist3.

    Args:
        df: a DataFrame of TCRs. Per chain it needs a CDR3 and V/J calls, looked for as
            ``{chain}_cdr3`` / ``{chain}_v`` / ``{chain}_j`` - the same names
            :func:`tcrkit.build_tcr_pmhc` produces.
        chains: ``('beta',)``, ``('alpha',)`` or both. With both, the returned distance
            is the sum over chains, which is what tcrdist3 calls the paired distance.
        col_mapping: where those fields live if named differently, e.g.
            ``{'beta_cdr3': 'junction_aa', 'beta_v': 'v_call', 'beta_j': 'j_call'}``.
        species: ``'human'`` or ``'mouse'`` (spellings like 'Mus musculus' understood).
        long: return one row per pair (``i``, ``j``, ``distance``) instead of a square
            matrix - easier to filter and join, much smaller for sparse use.
        **kwargs: forwarded to ``tcrdist.repertoire.TCRrep``.

    Returns:
        A square DataFrame of distances indexed and columned like `df`, or a long-form
        frame with ``long=True``.

    Note:
        tcrdist3 wants V/J calls **with an allele** (``TRBV19*01``, not ``TRBV19``); this
        appends ``*01`` where one is missing, which is what tcrdist3's own examples do.
        Rows whose genes it cannot resolve are dropped by tcrdist3 itself, so check the
        returned shape against your input.

        Two warnings tcrdist3 emits on every stock call are suppressed because neither
        says anything about your data - see :func:`_quiet_tcrdist`. Any other warning it
        raises still reaches you.

    Requires the ``metrics`` extra (``pip install 'tcrkit[metrics]'``).
    """
    with _quiet_tcrdist():          # the tqdm/ipywidgets warning fires on import
        try:
            from tcrdist.repertoire import TCRrep
        except ImportError as e:                            # pragma: no cover
            raise _missing('tcrdist3') from e

    col_mapping = col_mapping or {}
    short = {'alpha': 'a', 'beta': 'b'}
    cells = pd.DataFrame(index=range(len(df)))

    for chain in chains:
        if chain not in short:
            raise ValueError(f"chains must be 'alpha' and/or 'beta', got {chain!r}.")
        s = short[chain]
        for field, target in (('cdr3', f'cdr3_{s}_aa'), ('v', f'v_{s}_gene'),
                              ('j', f'j_{s}_gene')):
            src = col_mapping.get(f'{chain}_{field}', f'{chain}_{field}')
            if src not in df.columns:
                raise ValueError(
                    f"Missing column {src!r} for the {chain} chain. Pass col_mapping, "
                    f"e.g. col_mapping={{'{chain}_{field}': 'your_column'}}."
                )
            values = df[src].astype(str).values
            if field in ('v', 'j'):
                # tcrdist3 matches against allele-level names
                values = [x if '*' in x else f'{x}*01' for x in values]
            cells[target] = values

    cells['count'] = 1
    kwargs.setdefault('deduplicate', False)
    kwargs.setdefault('compute_distances', True)
    with _quiet_tcrdist():
            tr = TCRrep(cell_df=cells, organism=species_of(species),
                    chains=list(chains), **kwargs)

    total = None
    for chain in chains:
        pw = getattr(tr, f'pw_{chain}')
        total = pw if total is None else total + pw

    matrix = pd.DataFrame(total, index=df.index, columns=df.index)
    if not long:
        return matrix
    pairs = matrix.stack().rename('distance').reset_index()
    pairs.columns = ['i', 'j', 'distance']
    return pairs[pairs.i != pairs.j].reset_index(drop=True)
