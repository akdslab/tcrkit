"""Assemble complete TCR-pMHC complexes from raw V/J/CDR3 calls and MHC alleles.

This is the whole-complex builder: both TCR chains *and* the pMHC side, plus the
self-check between them. For the pMHC side on its own - an allele to its sequence - see
:mod:`tcrkit.pmhc`.

One function, :func:`build_tcr_pmhc`. It runs the whole chain of capabilities over a
table of TCR pairs and their MHC restriction, in the order they have to happen::

    1. normalize_tcr        V and J calls          -> IMGT
    2. add_cdr3_junctions   CDR3s                  -> with their 104 C / 118 F
    3. stitch               V + J + CDR3           -> full-length alpha and beta chains
    4. mhc_sequence         MHC allele(s)          -> both chains of the heterodimer
    5. anarci               the stitched sequences -> V/J/CDR annotations

and then checks itself: ANARCI re-reads each stitched sequence with no knowledge of the
input, and its CDR3 is compared against the CDR3 that went to Stitchr. Every step's
complaint is collected into one ``error`` column and one ``ok`` flag, so a row that went
wrong says where.

Columns the caller brought along - an epitope, a study id, an affinity - are carried
through untouched, so the result is a superset of the input.
"""

from __future__ import annotations

import warnings

import pandas as pd

from ._util import ERROR_COL, species_of

__all__ = ['build_tcr_pmhc', 'TCR_FIELDS']

#: The per-chain input fields, looked for as ``{chain}_{field}`` unless `col_mapping`
#: says otherwise.
TCR_FIELDS = ('v', 'j', 'cdr3')

_CHAINS = ('alpha', 'beta')
_ANARCI_ALLOW = {'alpha': 'A', 'beta': 'B'}


def _resolve(df, col_mapping, key, default):
    """The column backing a logical field, or None when the frame does not carry it."""
    col = (col_mapping or {}).get(key, default)
    return col if col in df.columns else None


def _chain_columns(df, col_mapping, chain):
    """``{'v': col, 'j': col, 'cdr3': col}`` for a chain, or None if it is absent."""
    cols = {f: _resolve(df, col_mapping, f'{chain}_{f}', f'{chain}_{f}')
            for f in TCR_FIELDS}
    return None if any(c is None for c in cols.values()) else cols


def build_tcr_pmhc(df, col_mapping=None, species='human', run_anarci=True,
                   keep=None, mhc_table=None, progress_bar=True) -> pd.DataFrame:
    """Build full TCR-pMHC sequences from gene calls and MHC alleles.

    >>> raw = pd.DataFrame([{'alpha_v': 'TRAV1-2', 'alpha_j': 'TRAJ33',
    ...                      'alpha_cdr3': 'AVMDSNYQLI',
    ...                      'beta_v': 'TCRBV20S1', 'beta_j': 'TRBJ2-7',
    ...                      'beta_cdr3': 'ASSLGQAYEQY',
    ...                      'mhc_a': 'HLA-A*02:01', 'epitope': 'GILGFVFTL'}])
    >>> out = build_tcr_pmhc(raw, run_anarci=False, progress_bar=False)
    >>> out[['epitope', 'beta_v', 'beta_cdr3', 'ok']].to_dict('records')
    [{'epitope': 'GILGFVFTL', 'beta_v': 'TRBV20-1', 'beta_cdr3': 'CASSLGQAYEQYF', 'ok': True}]

    Args:
        df: a table of TCRs. Per chain it needs V, J and CDR3, looked for as
            ``alpha_v`` / ``alpha_j`` / ``alpha_cdr3`` and the ``beta_`` equivalents; a
            chain whose columns are absent is skipped, so single-chain tables work. MHC
            restriction comes from ``mhc_a`` (and optionally ``mhc_b`` for an explicit
            class-II pair).
        col_mapping: where those fields live if they are named differently, e.g.
            ``{'beta_v': 'v_call', 'beta_cdr3': 'junction_aa', 'mhc_a': 'mhc.a'}``.
        species: ``'human'`` or ``'mouse'`` (spellings like 'Mus musculus' understood),
            or the name of a column holding it per row.
        run_anarci: annotate the stitched sequences with ANARCI and check that the CDR3
            survives the round trip. This is the slow step - it numbers every stitched
            chain - so turn it off for a quick pass.
        keep: input columns to carry through. ``None`` (the default) keeps them all.
        mhc_table: your own MHC sequence reference - a DataFrame or CSV path with
            ``allele`` and ``sequence`` columns, allele spellings normalized for you.
            Without one the bundled reference is used, and if that is absent the MHC
            sequence columns are left empty (with a warning) rather than raising, so the
            TCR half of the build still happens.
        progress_bar: show progress for the stitching and numbering passes.

    Returns:
        A DataFrame indexed like `df`, with your input columns first, then per chain:

        ======================  ==================================================
        column                  meaning
        ======================  ==================================================
        {chain}_v, _j           the IMGT-standardized calls actually used
        {chain}_cdr3            the CDR3 actually stitched, junction completed
        {chain}_sequence        stitched nucleotides
        {chain}_sequence_aa     its translation
        {chain}_anarci_v, _j    what ANARCI independently read back
        {chain}_cdr3_ok         did ANARCI recover the CDR3 that went in?
        ======================  ==================================================

        plus ``mhc_a_allele`` / ``mhc_a_sequence`` / ``mhc_b_allele`` /
        ``mhc_b_sequence`` / ``mhc_class``, and finally ``error`` (every complaint,
        ' | '-joined, '-' when clean) and ``ok`` (nothing complained and every requested
        chain stitched).

        Nothing raises for bad data - a row that fails says so in ``error`` and the rest
        of the table is still built.
    """
    from .cdr3_junction import add_cdr3_junctions
    from .normalize import normalize_tcr
    from .pmhc import mhc_sequence
    from .stitch import stitch

    out = (df if keep is None else df[list(keep)]).copy()
    problems = [[] for _ in range(len(df))]

    def note(mask, message):
        """Record a complaint against each flagged row. `message` may be per row."""
        messages = message if isinstance(message, (list, tuple, pd.Series)) else None
        for i, bad in enumerate(mask):
            if bad:
                problems[i].append(messages[i] if messages is not None else message)

    species_col = _resolve(df, col_mapping, 'species', species)
    # one species vocabulary everywhere now, so the same value serves every step
    row_species = (list(df[species_col]) if species_col else species)
    tcr_species = row_species

    present = [c for c in _CHAINS if _chain_columns(df, col_mapping, c)]
    if not present:
        raise ValueError(
            "No chain found. Expected alpha_v / alpha_j / alpha_cdr3 or the beta_ "
            "equivalents; pass col_mapping to say where yours live, e.g. "
            "col_mapping={'beta_v': 'v_call', 'beta_j': 'j_call', "
            "'beta_cdr3': 'junction_aa'}."
        )

    for chain in present:
        cols = _chain_columns(df, col_mapping, chain)

        # 1. gene calls -> IMGT. Keep the original where tidytcells cannot place it,
        #    so the row still has something to stitch with and the error explains why.
        for field in ('v', 'j'):
            # log_failures=False: an unrecognised call is already reported per row in
            # `errors`, so tidytcells logging it as well is duplicate noise - and on a
            # frame of any size it buries the result.
            std = normalize_tcr(df, col_mapping={'symbol': cols[field]},
                                species=tcr_species,
                                log_failures=False)['standardized']
            note(std.isna() & df[cols[field]].notna(),
                 f'{chain} {field.upper()} call not recognised')
            out[f'{chain}_{field}'] = std.where(std.notna(), df[cols[field]]).values

        # 2. CDR3 -> with its junction residues, using the standardized J gene
        junc = add_cdr3_junctions(
            pd.DataFrame({'cdr3': df[cols['cdr3']].values,
                          'j': out[f'{chain}_j'].values}),
            keep_unknown=True)
        note(junc['junction'].isna() & df[cols['cdr3']].notna(),
             f'{chain} J gene unknown to the junction table')
        out[f'{chain}_cdr3'] = junc['cdr3_fixed'].values

        # 3. stitch
        st = stitch(pd.DataFrame({'v': out[f'{chain}_v'].values,
                                  'j': out[f'{chain}_j'].values,
                                  'cdr3': out[f'{chain}_cdr3'].values}),
                    chain=chain, species=row_species, progress_bar=progress_bar)
        note(st[ERROR_COL].ne('-'), f'{chain} did not stitch')
        out[f'{chain}_sequence'] = st['sequence'].values
        out[f'{chain}_sequence_aa'] = st['sequence_aa'].values

    # 4. MHC -> both chains of the heterodimer and their sequences
    mhc_a = _resolve(df, col_mapping, 'mhc_a', 'mhc_a')
    mhc_b = _resolve(df, col_mapping, 'mhc_b', 'mhc_b')
    if mhc_a:
        # an explicit alpha__beta pair when both columns are given, else let
        # mhc_sequence infer the partner (B2M for class I, DRA for DRB*)
        alleles = (df[mhc_a].astype(str) + '__' + df[mhc_b].astype(str)
                   if mhc_b else df[mhc_a])
        try:
            mhc = mhc_sequence(list(alleles), paired=True, table=mhc_table)
        except FileNotFoundError as e:
            # No sequence reference to hand. Normalize the alleles anyway and leave the
            # sequences empty, rather than losing the whole build over the MHC step.
            from .normalize import normalize_mhc

            warnings.warn(f'MHC sequences left empty: {e} Pass mhc_table= to supply '
                          'your own.', stacklevel=2)
            mhc = normalize_mhc(list(alleles), paired=True, errors='coerce')
            mhc['mhc_a_sequence'] = None
            mhc['mhc_b_sequence'] = None
        note(list(mhc[ERROR_COL].notna()),
             ['MHC: ' + str(e) for e in mhc[ERROR_COL]])
        for c in ('mhc_a_allele', 'mhc_a_sequence', 'mhc_b_allele', 'mhc_b_sequence',
                  'mhc_class'):
            out[c] = mhc[c].values

    # 5. ANARCI as an independent referee on what we just built
    if run_anarci:
        from .anarci import anarci

        for chain in present:
            aa = out[f'{chain}_sequence_aa']
            ok = aa.ne('-') & aa.notna()
            ann = anarci(list(aa.where(ok, '')), allow={_ANARCI_ALLOW[chain]},
                         progress_bar=progress_bar)
            ann.index = out.index
            out[f'{chain}_anarci_v'] = ann['v_gene'].where(ok)
            out[f'{chain}_anarci_j'] = ann['j_gene'].where(ok)
            # containment, not equality: ANARCI reports IMGT 105-117 while the CDR3 we
            # stitched with carries the conserved C104 / F118
            recovered = [isinstance(a, str) and isinstance(b, str)
                         and a not in ('-', '') and a in b
                         for a, b in zip(ann['CDR3_aa'], out[f'{chain}_cdr3'])]
            out[f'{chain}_cdr3_ok'] = pd.Series(recovered, index=out.index).where(ok)
            note([bool(o) and not r for o, r in zip(ok, recovered)],
                 f'{chain} CDR3 not recovered from the stitched sequence')

    out[ERROR_COL] = [' | '.join(p) if p else '-' for p in problems]
    out['ok'] = [not p for p in problems]
    return out
