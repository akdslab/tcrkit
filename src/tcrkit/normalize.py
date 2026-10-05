"""Standardize TCR gene/allele names and HLA/MHC allele names to canonical nomenclature.

Two capabilities that do the same kind of job - clean up the many spellings sources use
for the same gene or allele - so they live together:

    normalize_tcr   TCR V/J gene and allele symbols -> IMGT, via `tidytcells`
    normalize_mhc   MHC/HLA allele names -> canonical fields, via `mhcgnomes`

This module is about *names*. To go from an allele to its protein sequence see
:func:`tcrkit.mhc_sequence` in :mod:`tcrkit.pmhc`.

Both take a single value or any collection of them and return a DataFrame with the input
echoed back beside the result; ``squeeze=True`` gives the bare value instead::

    normalize_tcr(['TCRBV20S1', 'TRAV23/6'], species='human')
    normalize_mhc(df, col_mapping={'allele': 'mhc'}, errors='coerce')

:func:`genes_match` and :func:`mhc_chain_of` are plain predicates about a single symbol,
so those two return a bool / a string rather than a frame.
"""

from __future__ import annotations

import re
from functools import lru_cache

import pandas as pd

from ._util import ERROR_COL, collect, is_empty, per_row, rows as _rows
from ._util import species_of, species_or_default, squeeze as _squeeze

__all__ = [
    # TCR genes / alleles
    'normalize_tcr', 'to_gene', 'genes_match',
    # MHC / HLA alleles
    'normalize_mhc', 'mhc_chain_of',
]


# ==========================================================================================
# TCR gene / allele names (tidytcells)
# ==========================================================================================

# tidytcells' own spelling. Purely internal: callers say 'human' / 'mouse'.
_TT = {'human': 'homosapiens', 'mouse': 'musmusculus'}
_FROM_TT = {v: k for k, v in _TT.items()}

# Order 'any' searches in.
_ANY_SPECIES = ('human', 'mouse')


def _tidytcells():
    """tidytcells, imported on first use so ``import tcrkit`` stays cheap."""
    import tidytcells as tt
    return tt


def _species_key(species) -> str | None:
    """'human' / 'mouse' -> tidytcells' spelling; None for 'any'. Raises on a typo."""
    name = species_of(species, default='any', allow=('any',))
    return None if name == 'any' else _TT[name]


def _safe_species(value):
    """As :func:`_species_key` but free text falls back to 'any' instead of raising."""
    return species_or_default(value, default='any', allow=('any',))


def _standardize(symbol, species='any', return_species=False, return_gene=False,
                 **kwargs):
    """One symbol -> the standardized symbol, with the original scalar return shape.

    The single-value core, kept scalar because :func:`tcrkit.add_cdr3_junctions` calls it per
    symbol on a hot path. :func:`normalize_tcr` calls it per row.
    """
    tt = _tidytcells()
    sp = _species_key(species)

    def pack(result, resolved_sp):
        extras = {}
        if return_species:
            extras['species'] = resolved_sp
        if return_gene:
            # resolved_sp is tcrkit's name ('human'); tidytcells wants its own
            gene = (tt.tr.standardize(symbol=symbol, species=_TT.get(resolved_sp,
                                                                     resolved_sp),
                                      **{**kwargs, 'precision': 'gene'})
                    if isinstance(symbol, str) and symbol and resolved_sp else None)
            if gene and result and gene not in str(result):
                raise AssertionError(f"gene={gene!r} not a substring of result={result!r}")
            extras['gene'] = gene
        return (result, extras) if extras else result

    if not isinstance(symbol, str) or is_empty(symbol):
        return pack(None, None)

    if sp is None:
        search_kwargs = {**kwargs, 'on_fail': 'reject', 'log_failures': False}
        for candidate in _ANY_SPECIES:
            result = tt.tr.standardize(symbol=symbol, species=_TT[candidate],
                                       **search_kwargs)
            if result is not None:
                return pack(result, candidate)
        # Nothing resolved. Honour an explicit on_fail='keep' by asking once more.
        if kwargs.get('on_fail') == 'keep':
            return pack(tt.tr.standardize(symbol=symbol, species=_TT[_ANY_SPECIES[0]],
                                          **{**kwargs, 'log_failures': False}), None)
        return pack(None, None)

    return pack(tt.tr.standardize(symbol=symbol, species=sp, **kwargs),
                _FROM_TT.get(sp, sp))


def normalize_tcr(symbols, species=None, col_mapping=None, with_gene=False,
                  squeeze: bool = False, **kwargs) -> pd.DataFrame:
    """Standardize TCR gene/allele symbols to IMGT nomenclature.

    >>> normalize_tcr('TCRBV20S1', species='human')[['symbol', 'standardized']]
       symbol standardized
    0  TCRBV20S1     TRBV20-1

    >>> list(normalize_tcr(['TCRBV20S1', 'TRAV23/6'], species='human')['standardized'])
    ['TRBV20-1', 'TRAV23/DV6']

    Args:
        symbols: a single symbol, or any collection of them - a list, Series, dict,
            list of records, or a DataFrame. See :func:`tcrkit._util.collect` for the
            shapes recognised.
        species: 'human', 'mouse' (aliases such as 'Homo sapiens' accepted), or 'any' to
            try human then mouse and use whichever resolves. May also be a sequence with
            one entry per symbol. An unrecognised name raises. Left as ``None`` (the
            default) the input's own ``species`` column is used if it has one, else
            ``'any'`` - values in your own column are treated as free text and fall
            back to ``'any'`` rather than raising.
        col_mapping: where the fields live when the input is a frame or records, e.g.
            ``{'symbol': 'v_call', 'species': 'organism'}``.
        with_gene: also report the symbol truncated to the gene. Costs one extra
            tidytcells call per row.
        squeeze: return the plain result instead of a frame - the standardized string
            for one symbol, a Series of them for several.
        **kwargs: forwarded to ``tidytcells.tr.standardize``, e.g. ``precision='gene'``,
            ``on_fail='keep'``, ``log_failures=False``.

    Returns:
        A DataFrame indexed like the input, with the input echoed back and the result
        alongside it:

        ============  ==========================================================
        column        meaning
        ============  ==========================================================
        symbol        the symbol as given
        standardized  the IMGT symbol, or None where tidytcells could not place it
        species       the species it resolved under ('any' inputs report which won)
        gene          only with ``with_gene=True``
        ============  ==========================================================
    """
    inp = collect(('symbol',), (symbols,), col_mapping=col_mapping,
                  optional=('species',))
    if species is None:
        # read from the caller's own column: free text, so fall back rather than raise
        wanted = [_safe_species(v) for v in per_row(None, inp, 'species', 'any')]
    else:
        # passed explicitly: an unrecognised name is a typo, so say so
        wanted = per_row(species, inp, 'species', 'any', validate=_species_key)

    recs = []
    for r, sp in zip(_rows(inp), wanted):
        result, extras = _standardize(r['symbol'], sp, return_species=True,
                                      return_gene=with_gene, **kwargs)
        rec = {'standardized': result, 'species': extras.get('species'),
               ERROR_COL: (None if result is not None or is_empty(r['symbol'])
                           else 'not recognised by tidytcells')}
        if with_gene:
            rec['gene'] = extras.get('gene')
        recs.append(rec)

    out = pd.DataFrame(recs, index=inp.index,
                       columns=['standardized', 'species']
                               + (['gene'] if with_gene else []) + [ERROR_COL])
    # the caller's species column, if any, is superseded by the resolved one
    res = pd.concat([inp.drop(columns=['species'], errors='ignore'), out], axis=1)
    if squeeze:
        # the species/gene annotations are extras; squeeze to the standardization itself
        keep = ['standardized'] + (['gene'] if with_gene else [])
        return _squeeze(res[keep], ())
    return res


def to_gene(symbols, species=None, **kwargs) -> pd.DataFrame:
    """:func:`normalize_tcr` truncated to the gene; keeps the input on failure.

    >>> to_gene('TRBV20-1*01', species='human')['standardized'][0]
    'TRBV20-1'
    """
    kwargs.setdefault('precision', 'gene')
    kwargs.setdefault('on_fail', 'keep')
    return normalize_tcr(symbols, species=species, **kwargs)


@lru_cache(maxsize=None)
def _roots(symbol: str) -> str:
    """The bare gene root: 'TRAV23/DV6*01' -> 'TRAV23'."""
    s = str(symbol).strip().upper().split('*')[0]   # drop the allele suffix
    return re.split(r'[/\-]', s)[0]


def genes_match(gene, symbol) -> bool:
    """True if `gene` and `symbol` refer to the same gene.

    >>> genes_match('TRAV23/DV6', 'TRAV23/6')
    True
    >>> genes_match('TRAV10', 'TRAV10*01')
    True

    Tolerates allele suffixes on `symbol`, IMGT's TRAV/DV dual naming, and subgroup
    suffixes ('TRAV5-1' vs 'TRAV5*01'). A plain predicate, so this one returns a bool
    rather than a frame.
    """
    symbol = str(symbol)
    if str(gene) in symbol:
        return True
    return _roots(gene) in _roots(symbol)


# ==========================================================================================
# MHC / HLA allele names (mhcgnomes)
# ==========================================================================================

def _mhcgnomes():
    """mhcgnomes, imported on first use so ``import tcrkit`` stays cheap."""
    import mhcgnomes
    return mhcgnomes


# ---------------------------------------------------------------------------
# Input normalization
# ---------------------------------------------------------------------------

# Species-tagged B2M labels (the form this module emits for a class-I light chain)
# mapped back to something mhcgnomes understands, so output can be fed straight back
# in and round-trip.
_B2M_INPUT_LABELS = {
    'B2M_HUMAN': 'HLA-B2M',
    'B2M_MOUSE': 'H2-B2M',
}


def _fix_allele_str(allele_str: str) -> str:
    """Patch known mhcgnomes parsing quirks in the raw allele string."""
    # HLA-DPA*01:03 is missing its gene number; mhcgnomes expects DPA1.
    allele_str = str(allele_str).replace('DPA*', 'DPA1*')
    return _B2M_INPUT_LABELS.get(allele_str.strip().upper(), allele_str)


# ---------------------------------------------------------------------------
# Chain resolution (single vs. paired)
# ---------------------------------------------------------------------------

def _resolve_chains(allele_str: str, paired: bool, chain: str | None):
    """Classify an allele string as a single chain or an alpha/beta pair.

    Returns ``('single', allele_str)`` or ``('pair', alpha_str, beta_str)``.
    """
    mhcgnomes = _mhcgnomes()

    # Explicit "alpha__beta" pairing always wins, except when both halves are the same
    # string: some sources repeat one whole-molecule name in both columns
    # (mhc.a = mhc.b = 'H2-IAb'), which is a single pair, not two copies of one chain.
    # Fall through and let mhcgnomes split it.
    if '__' in allele_str:
        alpha, beta = allele_str.split('__')
        if alpha != beta:
            return 'pair', alpha, beta
        allele_str = alpha

    parsed = mhcgnomes.parse(allele_str)

    # mhcgnomes recognised a native alpha/beta pair.
    if hasattr(parsed, 'alpha') and hasattr(parsed, 'beta'):
        if paired:
            return 'pair', parsed.alpha.to_string(), parsed.beta.to_string()
        if chain == 'alpha':
            return 'single', parsed.alpha.to_string()
        if chain == 'beta':
            return 'single', parsed.beta.to_string()
        raise ValueError(
            f"Allele {allele_str!r} is a pair - specify chain='alpha' or 'beta', "
            f"or use paired=True"
        )

    # A single chain we can pair with a known partner when paired=True.
    if paired:
        if parsed.gene.name.startswith('DRB'):
            return 'pair', 'DRA0101', allele_str
        if (parsed.gene.species.mhc_prefix == 'Gaga'
                and parsed.gene.name.startswith('BLB')):
            return 'pair', 'BLA', allele_str

        mhc_class = getattr(getattr(parsed, 'gene', None), 'mhc_class', '')
        gene_name = getattr(getattr(parsed, 'gene', None), 'name', '?')
        if 'II' in str(mhc_class):
            raise ValueError(
                f"Allele {allele_str!r} is class {mhc_class} (gene {gene_name!r}) but no "
                f"alpha-chain pairing rule exists - pass the pair explicitly as "
                f"'alpha__beta', or use paired=False"
            )

    return 'single', allele_str


# ---------------------------------------------------------------------------
# Field extraction
# ---------------------------------------------------------------------------

_B2M_BY_SPECIES = {'human': 'B2M_human', 'mouse': 'B2M_mouse'}


def _mhc_class_group(cls: str) -> str:
    """Collapse 'Ia'/'Ib' -> 'I' and 'IIa' -> 'II' for class-consistency checks."""
    cls = str(cls)
    if cls.startswith('II'):
        return 'II'
    if cls.startswith('I'):
        return 'I'
    return cls


# Class-II gene names end in the chain letter, optionally followed by a paralog number:
#   DRA / DQA1 / DPA1 / DMA / AA / EA / BLA      -> alpha
#   DRB1 / DQB1 / DPB1 / DOB / AB / AB1 / BLB2   -> beta
_ALPHA_GENE_SUFFIX = re.compile(r'A\d*$')
_BETA_GENE_SUFFIX = re.compile(r'B\d*$')


def _chain_from_fields(gene: str, mhc_class: str) -> str | None:
    """Classify a chain as 'alpha' / 'beta' from its gene name and MHC class.

    Returns None when the chain cannot be classified, so callers can leave an ambiguous
    pair in the order it was given rather than guessing.
    """
    if gene == 'B2M':
        return 'beta'                      # class-I light chain

    cls = _mhc_class_group(mhc_class)
    if cls == 'I':
        return 'alpha'                     # class-I heavy chain

    # Only test the gene suffix for class II: class-I genes are themselves named 'A' and
    # 'B' (HLA-A, HLA-B), which the suffixes would misread.
    if cls == 'II':
        if _ALPHA_GENE_SUFFIX.search(gene):
            return 'alpha'
        if _BETA_GENE_SUFFIX.search(gene):
            return 'beta'
    return None


def mhc_chain_of(allele_str: str) -> str | None:
    """Return 'alpha', 'beta', or None for a single-chain MHC allele string.

    >>> mhc_chain_of('HLA-A*02:01'), mhc_chain_of('HLA-DRB1*04:01')
    ('alpha', 'beta')

    Class I: the heavy chain (A/B/C/E/K/D/...) is 'alpha', B2M is 'beta'.
    Class II: decided by the gene's chain letter (DRA -> alpha, DRB1 -> beta).
    """
    parsed = _parse_single(_fix_allele_str(allele_str), unified=False,
                          enforce_one_allele=False, restrict_allele_fields=2,
                          drop_mutations=False)
    return _chain_from_fields(parsed['gene'], parsed['mhc_class'])


def _parse_single(allele_str: str, unified: bool, enforce_one_allele: bool,
                  restrict_allele_fields: int, drop_mutations: bool) -> dict:
    """Extract structured fields from a single-chain allele string."""
    mhcgnomes = _mhcgnomes()
    allele = mhcgnomes.parse(allele_str).restrict_allele_fields(restrict_allele_fields)

    mutations = None
    if drop_mutations:
        mutations = allele.mutation_string()   # e.g. 'R9H' ('' when none)
        allele = allele.copy(mutations=())

    if enforce_one_allele and hasattr(allele, 'alleles'):
        raise ValueError(f"Allele {allele_str!r} has multiple possible alleles")

    compact = f'{allele.species.mhc_prefix}-{allele.compact_string()}'.replace('HLA-', '')
    gene = allele.gene
    mhc_class = gene.mhc_class.replace('Ia', 'I').replace('IIa', 'II')
    # mhcgnomes says 'Homo sapiens'; tcrkit reports 'human' everywhere
    species = species_or_default(gene.species.name, default=gene.species.name)

    # `unified` expresses every allele in the paired layout (a-chain plus a B2M b-chain
    # for class I); otherwise a flat single-allele layout is used.
    if unified:
        result = {
            'mhc_a_allele': compact,
            'mhc_b_allele': (_B2M_BY_SPECIES.get(species, 'B2M')
                             if _mhc_class_group(mhc_class) == 'I' else None),
            'mhc_prefix': gene.species.mhc_prefix,
            'species': species,
            'gene': gene.name,
            'mhc_class': mhc_class,
        }
    else:
        result = {
            'mhc_prefix': gene.species.mhc_prefix,
            'species': species,
            'gene': gene.name,
            'allele': compact,
            'mhc_class': mhc_class,
        }

    if drop_mutations:
        result['mhc_mutations'] = mutations
    return result


def _parse_pair(alpha_str: str, beta_str: str, enforce_one_allele: bool,
                restrict_allele_fields: int, drop_mutations: bool = False,
                fix_chain_order: bool = True) -> dict:
    """Combine two single-chain parses into a paired-allele dict.

    With fix_chain_order on (the default), the two chains are assigned to the a/b slots
    by gene identity rather than by input order, so a reversed pair such as
    ``'DRB5*0101__DRA*0101'`` comes back in the canonical alpha/beta order. A pair whose
    two halves are the same chain is an error; a pair with an unclassifiable half is
    left in the order it was given.

    Descriptive fields (prefix/species/gene/class) come from the polymorphic chain: beta
    for class II, alpha (the heavy chain) for class I, since a class-I beta is the
    invariant B2M and describes nothing.
    """
    alpha = _parse_allele(alpha_str, paired=False,
                          enforce_one_allele=enforce_one_allele,
                          restrict_allele_fields=restrict_allele_fields,
                          drop_mutations=drop_mutations)
    beta = _parse_allele(beta_str, paired=False,
                         enforce_one_allele=enforce_one_allele,
                         restrict_allele_fields=restrict_allele_fields,
                         drop_mutations=drop_mutations)

    swapped = False
    if fix_chain_order:
        chain_a = _chain_from_fields(alpha['gene'], alpha['mhc_class'])
        chain_b = _chain_from_fields(beta['gene'], beta['mhc_class'])
        if chain_a is not None and chain_a == chain_b:
            raise ValueError(
                f"MHC chain conflict: {alpha_str!r} (gene {alpha['gene']!r}) and "
                f"{beta_str!r} (gene {beta['gene']!r}) are both {chain_a} chains - "
                f"the pair is missing one half"
            )
        if chain_a == 'beta' and chain_b == 'alpha':
            alpha_str, beta_str = beta_str, alpha_str
            alpha, beta = beta, alpha
            swapped = True

    # B2M parses as class 'other', so only compare classes between two real MHC chains;
    # a B2M beta is valid alongside any class-I alpha.
    is_b2m = beta['gene'] == 'B2M'
    cls_a = _mhc_class_group(alpha['mhc_class'])
    cls_b = _mhc_class_group(beta['mhc_class'])
    if is_b2m:
        if cls_a != 'I':
            raise ValueError(
                f"B2M pairs only with a class-I heavy chain, but alpha {alpha_str!r} "
                f"is class {alpha['mhc_class']}"
            )
    elif cls_a != cls_b:
        raise ValueError(
            f"MHC class mismatch: alpha {alpha_str!r} is class {alpha['mhc_class']} "
            f"but beta {beta_str!r} is class {beta['mhc_class']}"
        )

    # The two halves of a heterodimer must come from the same species. An unqualified
    # 'B2M' is exempt: it names no species (mhcgnomes defaults it to human) and inherits
    # the alpha chain's below. A species-tagged B2M is checked like any other chain.
    if not (is_b2m and str(beta_str).strip().upper() == 'B2M'):
        if alpha['species'] != beta['species']:
            raise ValueError(
                f"MHC species mismatch: alpha {alpha_str!r} is {alpha['species']} "
                f"but beta {beta_str!r} is {beta['species']}"
            )

    # Species-tag a bare B2M from the alpha chain: 'B2M' alone parses as human
    # regardless of its partner, which would mislabel mouse class-I pairs.
    beta_allele = beta['allele']
    if is_b2m:
        beta_allele = _B2M_BY_SPECIES.get(alpha['species'], beta_allele)

    descriptive = alpha if is_b2m else beta
    result = {
        'mhc_a_allele': alpha['allele'],
        'mhc_b_allele': beta_allele,
        'mhc_prefix': descriptive['mhc_prefix'],
        'species': descriptive['species'],
        'gene': descriptive['gene'],
        'mhc_class': descriptive['mhc_class'],
    }
    if fix_chain_order:
        result['chain_order_swapped'] = swapped
    if drop_mutations:
        result['mhc_mutations'] = ' '.join(
            m for m in (alpha['mhc_mutations'], beta['mhc_mutations']) if m)
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _parse_allele(allele, paired: bool = False, chain: str | None = None,
                  unified: bool = False, enforce_one_allele: bool = False,
                  restrict_allele_fields: int = 2, drop_mutations: bool = False,
                  fix_chain_order: bool = True) -> dict:
    """One allele -> a dict of structured fields. The single-value core.

    >>> res = normalize_mhc('A*0201')
    >>> res['allele'], res['gene'], res['mhc_class'], res['species']
    ('A0201', 'A', 'I', 'human')

    Args:
        allele: allele name, or an explicit ``"alpha__beta"`` pair.
        paired: return an alpha/beta pair, inferring the partner chain where a rule
            exists (DRB* -> DRA0101, class I -> B2M).
        chain: 'alpha' or 'beta' - pick one chain of a native pair when paired=False.
        unified: express single alleles in the paired layout (a-chain + B2M b-chain)
            without inferring a class-II partner. ``paired=True`` turns this on.
        enforce_one_allele: raise if the string resolves to multiple alleles.
        restrict_allele_fields: how many allele fields to keep (2 -> A*02:01).
        drop_mutations: strip mutations from the allele and record them in an
            ``'mhc_mutations'`` field (e.g. 'R9H'; '' when the allele has none).
        fix_chain_order: assign the halves of a pair to the a/b slots by gene identity
            instead of input order, so a reversed pair such as ``'DRB5*0101__DRA*0101'``
            is returned canonically. Adds a ``'chain_order_swapped'`` field to pair
            results recording whether the halves were exchanged. Set False to keep the
            input order and surface upstream ordering errors as-is (no field is added).

    Returns:
        For a single chain: ``{'mhc_prefix', 'species', 'gene', 'allele', 'mhc_class'}``.
        For a pair (or ``unified=True``): ``{'mhc_a_allele', 'mhc_b_allele',
        'mhc_prefix', 'species', 'gene', 'mhc_class'}``.
    """
    allele_str = _fix_allele_str(allele)

    kind, *chains = _resolve_chains(allele_str, paired, chain)
    if kind == 'pair':
        return _parse_pair(*chains,
                           enforce_one_allele=enforce_one_allele,
                           restrict_allele_fields=restrict_allele_fields,
                           drop_mutations=drop_mutations,
                           fix_chain_order=fix_chain_order)
    # `paired` implies the paired layout. A class-I allele has no second MHC chain to
    # pair with, so _resolve_chains hands it back as a single - but the caller asked for
    # a/b columns, so express it that way (heavy chain + its invariant B2M) rather than
    # silently returning the flat layout for some alleles and the paired one for others.
    return _parse_single(chains[0], unified or paired, enforce_one_allele,
                         restrict_allele_fields, drop_mutations)


def normalize_mhc(alleles, paired: bool = False, chain: str | None = None,
                  errors: str = 'raise', unified: bool = False,
                  enforce_one_allele: bool = False, col_mapping=None,
                  squeeze: bool = False, **kwargs) -> pd.DataFrame:
    """Standardize MHC/HLA allele names. Returns one row per input allele.

    >>> normalize_mhc('A*0201')[['allele_input', 'allele', 'gene', 'mhc_class']]
      allele_input allele gene mhc_class
    0       A*0201  A0201    A         I

    >>> list(normalize_mhc(['A*0201', 'HLA-B57:01'])['allele'])
    ['A0201', 'B5701']

    Args:
        alleles: a single allele name, an explicit ``"alpha__beta"`` pair, or any
            collection of them - a list, Series, dict, list of records, or a DataFrame.
        paired: return the alpha/beta layout (``mhc_a_allele`` / ``mhc_b_allele``),
            inferring the partner chain where a rule exists: DRB* pairs with DRA0101,
            and a class-I heavy chain with its invariant B2M light chain. Implies
            ``unified``, so every allele comes back in the same two-column shape.
        chain: 'alpha' or 'beta' - pick one chain of a native pair when paired=False.
        errors: what to do when an allele cannot be parsed. 'raise' (default) re-raises.
            'coerce' / 'ignore' leaves that row's fields empty and records why in a
            trailing ``error`` column, so one bad allele does not abort a batch::

                out = normalize_mhc(df['mhc'], errors='coerce')
                out[out.error.notna()].error.value_counts()

        unified: express single alleles in the paired layout (a-chain + B2M b-chain)
            without inferring a class-II partner. ``paired=True`` turns this on.
        enforce_one_allele: raise if a string resolves to multiple alleles.
        col_mapping: where the allele lives when the input is a frame or records, e.g.
            ``{'allele': 'mhc'}``.
        squeeze: return the plain result instead of a frame - a dict for one allele, a
            DataFrame without the ``allele_input`` echo for several.
        **kwargs: forwarded per allele - ``restrict_allele_fields`` (how many fields to
            keep, 2 -> A*02:01), ``drop_mutations``, ``fix_chain_order``.

    Returns:
        A DataFrame indexed like the input. ``allele_input`` echoes what you passed;
        the rest are the parsed fields - ``mhc_prefix``, ``species``, ``gene``,
        ``mhc_class``, and either ``allele`` (single chain) or ``mhc_a_allele`` /
        ``mhc_b_allele`` (pairs, or ``unified=True``).

        A ``chain_order_swapped`` column appears whenever any input was parsed as a pair
        with chain-order fixing on; it is False for non-pairs and for rows that failed.
        ``out.chain_order_swapped.sum()`` counts how many arrived with their chains
        reversed.
    """
    if errors not in ('raise', 'coerce', 'ignore'):
        raise ValueError(f"errors must be 'raise', 'coerce' or 'ignore', got {errors!r}.")

    inp = collect(('allele',), (alleles,), col_mapping=col_mapping)
    values = list(inp['allele'])

    def safe_parse(a):
        try:
            return _parse_allele(a, paired=paired, chain=chain, unified=unified,
                                 enforce_one_allele=enforce_one_allele, **kwargs)
        except Exception as e:
            if errors == 'raise':
                raise
            return {ERROR_COL: f'{type(e).__name__}: {e}'}

    df = pd.DataFrame([safe_parse(a) for a in values], index=inp.index)
    if 'mhc_mutations' in df.columns:
        # empty-string, not NaN, for unparsed rows
        df['mhc_mutations'] = df['mhc_mutations'].fillna('').replace('', '-')
    if 'chain_order_swapped' in df.columns:
        # False, not NaN, for the single-chain and unparsed rows, so the column stays a
        # real bool and can be summed.
        df['chain_order_swapped'] = df['chain_order_swapped'].fillna(False).astype(bool)
    # Always present, whatever `errors` is set to, so the column set never depends on
    # the data or the mode - and always last, so it reads as an annotation.
    df[ERROR_COL] = df.pop(ERROR_COL) if ERROR_COL in df.columns else pd.NA
    df.insert(0, 'allele_input', values)
    return _squeeze(df, ('allele_input',)) if squeeze else df
