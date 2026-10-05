r"""Rebuild full-length TCR nucleotide / amino-acid sequences from V, J and CDR3 calls.

One function, :func:`stitch`. It takes a single rearrangement or any collection of them
and always returns a DataFrame::

    stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta')
    stitch([('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF'), ...])
    stitch(df, col_mapping={'v': 'v_call', 'j': 'j_call', 'cdr3': 'junction_aa'})

Stitchr needs the CDR3 to carry its IMGT junction residues (104 ``C`` ... 118
``F``/``W``/``C``/``L``/``V``). A CDR3 reported as 105-117 only will not stitch - repair
it first with :func:`tcrkit.add_cdr3_junctions`, as a separate step::

    df['cdr3'] = tcrkit.add_cdr3_junctions(df, col_mapping={'cdr3': 'cdr3_aa',
                                                      'j': 'j_call'})['cdr3_fixed']
    out = tcrkit.stitch(df, chain='beta')

This module does not repair anything itself: it stitches exactly the V, J and CDR3 it is
given.

``engine='thimble'`` runs the same job through Stitchr's ``thimble`` CLI - one subprocess
and one IMGT load for the whole input rather than a Python call per row. Both engines use
the same constant-region defaults, so their output is directly comparable.

:func:`stitch_pair` is the one genuinely different shape: it stitches *both* chains of a
receptor from one row, in a single thimble call.

Stitchr needs IMGT germline data downloaded once per environment::

    stitchrdl -s human            # and -s mouse, etc.

Supported species: https://jamieheather.github.io/stitchr/inputdata.html#species-covered
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import pandas as pd

from ._util import (
    ERROR_COL,
    collect,
    fmt_warnings,
    is_empty,
    mute_stderr,
    mute_stdout,
    per_row,
    progress,
    species_of,
    rows as _rows,
    squeeze as _squeeze,
)

__all__ = ['stitch', 'stitch_pair', 'DEFAULT_CONSTANT', 'SUPPORTED_CHAINS',
           'ANARCI_COLS']


# ==========================================================================================
# Constants
# ==========================================================================================

# Stitchr names its data directories in upper case. Purely internal: callers say
# 'human' / 'mouse' like everywhere else in tcrkit.
_STITCHR_SPECIES = {'human': 'HUMAN', 'mouse': 'MOUSE'}

_CHAIN_TO_LOCUS = {'alpha': 'TRA', 'beta': 'TRB', 'gamma': 'TRG', 'delta': 'TRD'}
SUPPORTED_CHAINS = tuple(_CHAIN_TO_LOCUS)

#: Constant region used when the caller does not name one - the convention this
#: package was prepared with. Note these are flat, whereas Stitchr's own rule is
#: J-gene aware (``TRBJ2*`` rearrangements get ``TRBC2*01``, not ``TRBC1*01``), so the
#: two differ for beta chains using TRBJ2. Pass ``constant='auto'`` to use Stitchr's
#: rule instead; see :func:`stitch`.
DEFAULT_CONSTANT = {'alpha': 'TRAC*01', 'beta': 'TRBC1*01',
                    'gamma': 'TRGC1*01', 'delta': 'TRDC*01'}

#: ``constant=`` values that name a policy rather than an allele.
CONSTANT_AUTO = 'auto'                       # let Stitchr's autofill_input decide
_CONSTANT_DEFAULT = ('constant', 'default')  # use DEFAULT_CONSTANT, same as None

_ANARCI_CHAIN = {'alpha': 'A', 'beta': 'B', 'gamma': 'G', 'delta': 'D'}

# thimble stitches a *pair* of loci at a time; these are the two receptors it knows.
_THIMBLE_RECEPTORS = {('alpha', 'beta'): ('TRA', 'TRB'),
                      ('gamma', 'delta'): ('TRG', 'TRD')}

# IMGT reference data is a disk read per (chain, species); module level means one load
# per process rather than one per call.
_IMGT_CACHE: dict = {}

#: ANARCI columns kept when re-reading a stitched sequence (``run_anarci=True``).
ANARCI_COLS = ('CDR1_aa', 'CDR2_aa', 'CDR3_aa', 'v_gene', 'j_gene', 'species', 'chain_type',
               'error',
               'CDR1_aa_start', 'CDR1_aa_end', 'CDR2_aa_start', 'CDR2_aa_end',
               'CDR3_aa_start', 'CDR3_aa_end')

_REQUIRED_COLS = ('v', 'j', 'cdr3')

# Values that mean 'not reported' and must reach thimble as blanks.
_EMPTY_TOKENS = ('', '-', 'nan', 'none', 'na', 'n/a', '<na>', 'null')


# ==========================================================================================
# Shared helpers
# ==========================================================================================

def _stitchr():
    """The Stitchr modules, imported on first use so ``import tcrkit`` stays cheap."""
    from Stitchr import stitchr as st
    from Stitchr import stitchrfunctions as fxn
    return st, fxn


def _normalise_species(species) -> str:
    """Any spelling -> the upper-case name Stitchr's data directories use."""
    return _STITCHR_SPECIES[species_of(species)]


def _check_chain(chain: str) -> str:
    if chain not in _CHAIN_TO_LOCUS:
        raise ValueError(f"chain must be one of {SUPPORTED_CHAINS}, got {chain!r}.")
    return chain


def _check_col_mapping(col_mapping, chain: str = '') -> dict:
    """Validate a caller-supplied ``{'v': ..., 'j': ..., 'cdr3': ...}`` mapping.

    There is deliberately no default: tcrkit does not assume any particular column
    naming convention, so the caller always says where V / J / CDR3 live.
    """
    where = f" for chain {chain!r}" if chain else ''
    if not col_mapping:
        raise ValueError(
            f"col_mapping is required{where}: pass e.g. "
            "{'v': 'v_call', 'j': 'j_call', 'cdr3': 'junction_aa'}"
        )
    missing = [k for k in _REQUIRED_COLS if k not in col_mapping]
    if missing:
        raise ValueError(f"col_mapping{where} is missing keys {missing}; "
                         f"needs all of {list(_REQUIRED_COLS)}.")
    return {k: col_mapping[k] for k in _REQUIRED_COLS}


def _get_imgt_cached(chain: str, species: str, cache: dict | None = None):
    """The Stitchr reference data for a chain/species, loaded once per process.

    Returns ``(tcr_dat, functionality, partial_genes, codons, c_motifs, j_residues,
    low_confidence_js)`` - everything ``Stitchr.stitchr.stitch`` needs besides the
    per-rearrangement arguments.
    """
    st, fxn = _stitchr()
    cache = _IMGT_CACHE if cache is None else cache
    key = (chain, species)
    if key not in cache:
        tcr_dat, functionality, partial_genes = fxn.get_ref_data(
            _CHAIN_TO_LOCUS[chain], st.gene_types, species)
        codons = fxn.get_optimal_codons('', species)
        # The conserved CDR3-ending residue per J gene, and the J genes for which it
        # cannot be confidently identified.
        j_residues, low_confidence_js = fxn.get_j_motifs(species)
        # N-terminal C gene translations, used to infer the correct reading frame.
        c_motifs = fxn.get_c_motifs(species)
        cache[key] = (tcr_dat, functionality, partial_genes, codons,
                      c_motifs, j_residues, low_confidence_js)
    return cache[key]


def _resolve_constant(constant, chain):
    """The constant-region allele to stitch with, or '' to let Stitchr choose.

    ``None`` / ``'constant'`` / ``'default'`` -> :data:`DEFAULT_CONSTANT`, this package's
    fixed convention. ``'auto'`` -> '' , which makes :func:`_stitch_one` hand the choice
    to Stitchr's own J-gene-aware rule. A dict selects per chain. Anything else is taken
    as the allele itself.
    """
    if isinstance(constant, dict):
        constant = constant.get(chain)
    if constant is None:
        return DEFAULT_CONSTANT[chain]
    if isinstance(constant, str):
        key = constant.strip().lower()
        if key in _CONSTANT_DEFAULT:
            return DEFAULT_CONSTANT[chain]
        if key == CONSTANT_AUTO:
            return ''
    return constant


def _stitch_one(v, j, cdr3, chain, species, no_leader=False, constant=None,
                mute=True, cache=None):
    """Stitch one rearrangement. Returns (nt, aa, warnings); raises if Stitchr fails.

    `species` must already be normalised. This is the single place the `tcr_bits` dict
    Stitchr expects is built, so every path hands it identical input.
    """
    st, fxn = _stitchr()
    (tcr_dat, functionality, partial_genes, codons,
     c_motifs, j_residues, low_confidence_js) = _get_imgt_cached(chain, species, cache)

    bits = {
        'v': v,
        'j': j,
        'l': '' if no_leader else v,
        'c': _resolve_constant(constant, chain),
        'cdr3': cdr3,
        'no_leader': bool(no_leader),
        'skip_c_checks': False,
        'skip_n_checks': False,
        'species': species,
        'seamless': False,
        'mode': 'BOTH_FA',
        '5_prime_seq': '', '3_prime_seq': '', 'name': 'TCR',
    }

    if not bits['c']:
        # constant='auto': reuse Stitchr's own autofill rather than copying its table,
        # so the two cannot drift. It fills the leader from V as well, which would undo
        # no_leader, so put that back afterwards.
        keep_leader = bits['l']
        try:
            fxn.autofill_input(bits, _CHAIN_TO_LOCUS[chain])
        except IOError as e:
            raise ValueError(
                f"constant='auto' cannot infer a constant region for species "
                f"{species!r}: {e} Pass the allele explicitly, e.g. "
                f"constant={DEFAULT_CONSTANT[chain]!r}."
            ) from e
        bits['l'] = keep_leader

    # record=True collects Stitchr's warnings instead of printing them; on failure they
    # are attached to the exception, since that is where the useful ones turn up.
    with mute_stdout(enabled=mute), mute_stderr(enabled=mute), \
            warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            stitched = st.stitch(bits, tcr_dat, functionality, partial_genes, codons,
                                 3, {}, c_motifs, j_residues, low_confidence_js)
        except Exception as e:
            e.stitchr_warnings = fmt_warnings(caught)
            raise
    if not mute:                             # asked to see them - put them back on stderr,
        for w in {str(x.message): x for x in caught}.values():   # once each
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)

    # translation_offset is the offset into the first codon; pad it out so the frame
    # translates.
    nt = stitched['stitched_nt']
    return (nt, fxn.translate_nt('N' * stitched['translation_offset'] + nt),
            fmt_warnings(caught))


def _anarci_annotate(seqs, chain, cdr3=None, cols=ANARCI_COLS, progress_bar=False,
                     **anarci_kwargs):
    """Re-read stitched sequences with ANARCI. Returns `*_anarci`-suffixed columns.

    Restricted to the matching chain type by default (alpha -> 'A', beta -> 'B').
    `cols=None` keeps everything ANARCI reports. Given `cdr3` (the CDR3 that was
    stitched, aligned with `seqs`), adds `cdr3_match_anarci`: whether ANARCI re-read
    that same CDR3 out of the stitched sequence. It is a containment test, not equality,
    because ANARCI reports IMGT 105-117 while the input CDR3 normally carries the
    conserved C104 / F118.
    """
    # Imported from the submodule rather than `from . import anarci`, because the package
    # __init__ rebinds the name `tcrkit.anarci` to the function.
    from .anarci import anarci as _anarci

    anarci_kwargs.setdefault('allow', (_ANARCI_CHAIN[chain],))
    res = _anarci(seqs, progress_bar=progress_bar, **anarci_kwargs)
    res = res.drop(columns=['sequence'], errors='ignore')   # drop the input echo
    keep = [c for c in (cols or res.columns) if c in res.columns]
    for _pos in ('CDR1_aa_start', 'CDR1_aa_end', 'CDR2_aa_start', 'CDR2_aa_end',
                 'CDR3_aa_start', 'CDR3_aa_end'):
        if _pos in res.columns and _pos not in keep:
            keep.append(_pos)
    out = res[keep].rename(
        columns=lambda c: f"{c[7:] if c.startswith('anarci_') else c}_anarci")

    if cdr3 is not None and 'CDR3_aa_anarci' in out:
        used = pd.Series(cdr3).reindex(out.index).astype(str)
        got = out['CDR3_aa_anarci'].astype(str)
        out['cdr3_match_anarci'] = [g not in ('-', 'nan') and g in u
                                    for g, u in zip(got, used)]
    return out


_RESULT_KEYS = ('sequence', 'sequence_aa', ERROR_COL, 'stitchr_warnings')


# ==========================================================================================
# In-process stitching
# ==========================================================================================

def _stitch_record(v, j, cdr3, chain='beta', species='human', no_leader=False,
                   constant=None, mute=True):
    """Stitch one rearrangement. Returns a dict. The single-value core of :func:`stitch`."""
    res = dict.fromkeys(_RESULT_KEYS, '-')
    if any(is_empty(x) for x in (v, j, cdr3)):
        return {**res, ERROR_COL: 'missing v/j/cdr3 input'}
    try:
        nt, aa, warned = _stitch_one(v, j, str(cdr3).strip(), chain,
                                     _normalise_species(species), no_leader=no_leader,
                                     constant=constant, mute=mute)
    except Exception as e:
        return {**res, ERROR_COL: str(e),
                'stitchr_warnings': getattr(e, 'stitchr_warnings', '-')}
    return {**res, 'sequence': nt, 'sequence_aa': aa, 'stitchr_warnings': warned}


def stitch(v, j=None, cdr3=None, chain=None, species=None, col_mapping=None,
           no_leader=False, constant=None, mute=True, engine='stitchr',
           run_anarci=False, anarci_cols=ANARCI_COLS, anarci_kwargs=None,
           max_workers=1, batch_size=500, progress_bar=True,
           squeeze: bool = False, **engine_kwargs) -> pd.DataFrame:
    """Rebuild full-length TCR sequences from V/J calls + CDR3.

    Returns one row per rearrangement, whether you passed one or a million.

    >>> out = stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta')
    >>> out['sequence_aa'][0][:20]
    'MSNQVLCCVVLCFLGANTVD'

    >>> out = stitch([('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF'),
    ...               ('TRBV6-5', 'TRBJ2-7', 'CASSYSIRGSRGEQYF')], progress_bar=False)
    >>> list(out['error'])
    ['-', '-']

    Args:
        v, j, cdr3: the V call, J call and CDR3 amino-acid sequence. Either three plain
            values (or three parallel collections), or pass everything through `v` as a
            DataFrame, list of records, dict of columns, dict, or list of
            ``(v, j, cdr3)`` triples. A ``chain`` or ``species`` column in such a
            container is used per row, overriding the keyword.

            The CDR3 must already carry its IMGT junction residues; run
            :func:`tcrkit.add_cdr3_junctions` first if it does not.
        chain: 'alpha', 'beta', 'gamma' or 'delta'. May also be a sequence with one
            entry per rearrangement. Left as ``None`` (the default) the input's own
            ``chain`` column is used if it has one, else ``'beta'``.
        species: ``'human'`` or ``'mouse'``. Spellings like ``'Mus musculus'`` or
            ``'HOMO SAPIENS'`` are understood. Same rules as `chain`; defaults to
            ``'human'``.
        col_mapping: where the fields live when the input is a frame or records, e.g.
            ``{'v': 'v_call', 'j': 'j_call', 'cdr3': 'junction_aa'}``. tcrkit does not
            guess your column names - only exact ``v`` / ``j`` / ``cdr3`` matches are
            picked up automatically.
        no_leader: omit the V gene's leader sequence.
        constant: which constant region to stitch on.

            ================  ====================================================
            value             meaning
            ================  ====================================================
            ``None``          :data:`DEFAULT_CONSTANT` - this package's fixed
                              convention (the default)
            ``'constant'``    the same thing, spelled out; useful from the CLI
            ``'auto'``        let Stitchr decide, using its own J-gene-aware rule
            ``'TRBC2*01'``    that allele
            ``{'beta': ...}`` per chain
            ================  ====================================================

            ``None`` and ``'auto'`` differ for beta: :data:`DEFAULT_CONSTANT` is flat
            (always ``TRBC1*01``) while Stitchr picks ``TRBC1*01`` for ``TRBJ1*``
            rearrangements and ``TRBC2*01`` for ``TRBJ2*`` ones - a 5-residue,
            2-length difference in the stitched chain. Alpha is ``TRAC*01`` either
            way. ``'auto'`` only works for human and mouse; other species must name
            the allele.
        mute: keep Stitchr's own stdout/stderr chatter quiet.
        engine: ``'stitchr'`` (default) calls Stitchr in-process, one rearrangement at a
            time. ``'thimble'`` runs Stitchr's ``thimble`` CLI once for the whole input -
            one subprocess and one IMGT load instead of one Python call per row, which
            is faster for large inputs. See :func:`stitch_pair` for stitching both
            chains of a receptor together.
        run_anarci: re-read the stitched sequences with ANARCI (one batched call) and
            join the ``*_anarci`` columns on, including ``cdr3_match_anarci`` - whether
            ANARCI found the CDR3 that went in.
        anarci_cols: which ANARCI columns to keep; ``None`` keeps all.
        anarci_kwargs: extra keyword arguments for the ANARCI pass.
        max_workers: > 1 fans `batch_size` chunks out over processes (``'stitchr'``
            engine only).
        batch_size: rows per chunk when `max_workers` > 1.
        progress_bar: show a tqdm bar (only for more than one row).
        squeeze: return the plain result instead of a frame - a dict for one
            rearrangement, a DataFrame without the input echo for several.
        **engine_kwargs: passed to the chosen engine, e.g. ``thimble_bin=``,
            ``extra_genes=``, ``seamless=`` for ``engine='thimble'``.

    Returns:
        A DataFrame indexed like the input, starting with the ``v`` / ``j`` / ``cdr3``
        you gave, then:

        ================  ======================================================
        column            meaning
        ================  ======================================================
        sequence          stitched nucleotides, '-' if it failed
        sequence_aa       its translation, '-' if it failed
        error             '-' on success, else the error message
        stitchr_warnings  what Stitchr warned about, ' | '-joined, '-' if quiet
        ================  ======================================================

        A row that cannot be stitched reports the reason in ``error`` rather
        than raising, so one bad row never aborts a batch. 'Unable to locate C terminus
        of CDR3' means that CDR3 is missing its junction residues - see
        :func:`tcrkit.add_cdr3_junctions`.
    """
    if engine not in ('stitchr', 'thimble'):
        raise ValueError(f"engine must be 'stitchr' or 'thimble', got {engine!r}.")
    if engine == 'stitchr' and engine_kwargs:
        # Fail loudly rather than silently ignoring a typo - or a `fix=` left over from
        # before junction repair became its own step (see tcrkit.add_cdr3_junctions).
        raise TypeError(
            f"stitch() got unexpected keyword argument(s) {sorted(engine_kwargs)}. "
            f"Extra keywords are only forwarded to engine='thimble'."
            + ("\nTo complete a CDR3's junction residues, call tcrkit.add_cdr3_junctions "
               "first and stitch the result." if 'fix' in engine_kwargs else '')
        )

    inp = collect(('v', 'j', 'cdr3'), (v, j, cdr3), col_mapping=col_mapping,
                  optional=('chain', 'species'))
    chains = per_row(chain, inp, 'chain', 'beta')
    speciess = per_row(species, inp, 'species', 'human',
                       validate=species_of if species is not None else None)
    for c in set(chains):
        _check_chain(c)

    if engine == 'thimble':
        if len(set(chains)) > 1 or len(set(speciess)) > 1:
            raise ValueError(
                "engine='thimble' stitches one chain and one species per call; got "
                f"chains {sorted(set(chains))} and species {sorted(set(speciess))}. "
                "Split the input, or use the default engine='stitchr'."
            )
        out = _stitch_via_thimble(inp, chain=chains[0], species=speciess[0],
                                  constant=constant, **engine_kwargs)
    elif max_workers > 1 and len(inp) > batch_size:
        func = partial(stitch, chain=chain, species=species, no_leader=no_leader,
                       constant=constant, mute=mute, progress_bar=False,
                       run_anarci=False, **engine_kwargs)
        inp = inp.assign(chain=chains, species=speciess)   # carry them into the batches
        batches = [inp.iloc[i:i + batch_size] for i in range(0, len(inp), batch_size)]
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            out = pd.concat(progress(ex.map(func, batches), total=len(batches),
                                     desc='stitch (parallel)',
                                     enabled=progress_bar))
        out = out.drop(columns=[c for c in inp.columns if c in out.columns])
    else:
        label = chains[0] if len(set(chains)) == 1 else 'mixed'
        it = progress(list(zip(_rows(inp), chains, speciess)), total=len(inp),
                      desc=f'stitch {label}',
                      enabled=progress_bar and len(inp) > 1)
        records = [_stitch_record(r['v'], r['j'], r['cdr3'], chain=c, species=sp,
                                  no_leader=no_leader, constant=constant, mute=mute)
                   for r, c, sp in it]
        out = pd.DataFrame(records, index=inp.index).reindex(columns=list(_RESULT_KEYS))

    res = pd.concat([inp, out], axis=1)

    if run_anarci:
        ok = res['sequence_aa'].ne('-')
        ann = _anarci_annotate(res.loc[ok, 'sequence_aa'], chains[0],
                               cdr3=res.loc[ok, 'cdr3'], cols=anarci_cols,
                               progress_bar=progress_bar, **(anarci_kwargs or {}))
        ann = ann.reindex(res.index).fillna('-')
        if 'error_anarci' in ann:
            ann.loc[~ok, 'error_anarci'] = 'not stitched'
        res = res.join(ann)

    return _squeeze(res, inp.columns) if squeeze else res


def _stitch_via_thimble(inp, chain, species, constant=None, **kwargs):
    """Stitch a single-chain input frame through the thimble CLI. Returns result columns.

    thimble stitches a receptor (a chain *pair*) per call; a single-chain input simply
    leaves the partner blank, which thimble accepts. :func:`stitch_pair` is the way to
    get both chains in the one subprocess.
    """
    pair = ('alpha', 'beta') if chain in ('alpha', 'beta') else ('gamma', 'delta')
    frame = inp.rename(columns={'v': f'{chain}_v', 'j': f'{chain}_j',
                                'cdr3': f'{chain}_cdr3'})
    out = stitch_pair(
        frame, chains=pair, species=species,
        col_mapping={f'{chain}_{f}': f'{chain}_{f}' for f in _REQUIRED_COLS},
        constant=None if constant is None else {chain: constant}, **kwargs)
    return pd.DataFrame({
        'sequence': out[f'{chain}_sequence'],
        'sequence_aa': out[f'{chain}_sequence_aa'],
        ERROR_COL: out[ERROR_COL],
        'stitchr_warnings': out['thimble_log'],
    }, index=inp.index)




# ==========================================================================================
# Batch stitching through thimble (Stitchr's CLI)
# ==========================================================================================

def _thimble_bin(thimble_bin=None) -> str:
    """Path to the thimble executable: as given, else on PATH, else next to this python."""
    if thimble_bin:
        return str(thimble_bin)
    found = (shutil.which('thimble')
             or shutil.which('thimble', path=str(Path(sys.executable).parent)))
    if not found:
        raise FileNotFoundError(
            "thimble executable not found. Install stitchr in this environment "
            "(`pip install stitchr` + `stitchrdl -s human`) or pass thimble_bin=."
        )
    return found


def _thimble_species(species: str) -> str:
    """Species for thimble's `-s`: normalised when recognised, passed through otherwise.

    Unlike the in-process path this does not raise on unknown names - thimble supports
    every species with a directory in the Stitchr data dir (`stitchrdl -s <name>`) and
    validates the name itself, so 'Macaca mulatta' goes out as 'MACACAMULATTA'.
    """
    try:
        return _STITCHR_SPECIES[species_of(species)]
    except ValueError:
        # thimble supports any species with a directory in the Stitchr data dir and
        # validates the name itself, so pass an unrecognised one through in that form
        return re.sub(r'[^A-Z0-9]', '', str(species).upper())


def _thimble_col(df: pd.DataFrame, col: str) -> list:
    """One input column as blank-padded strings; a column that isn't there is all blanks."""
    if col not in df.columns:
        return [''] * len(df)
    s = df[col].astype(str).str.strip().fillna('')
    return s.mask(s.str.lower().isin(_EMPTY_TOKENS), '').tolist()


def stitch_pair(df, col_mapping, chains=('alpha', 'beta'), species='human',
                 constant=None, leader=None, extra_genes=False, seamless=False,
                 thimble_bin=None, in_file=None, out_file=None, work_dir=None,
                 keep_files=False, extra_args=(), verbose=False):
    """Stitch a whole frame in one `thimble` subprocess, returned indexed like `df`.

    Both chains of a pair go in together, which is the point of this over
    :func:`stitch`: one subprocess and one IMGT load for the frame rather than a
    per-row Python call. There are no retries - for species / CDR3-junction fixing use
    :func:`tcrkit.add_cdr3_junctions` followed by :func:`stitch`.

    Args:
        df: the frame to read from.
        col_mapping: where V / J / CDR3 live, with the chain as a prefix - the same
            convention :func:`tcrkit.build_tcr_pmhc` and :func:`tcrkit.tcrdist` use::

                {'alpha_v': 'a_v', 'alpha_j': 'a_j', 'alpha_cdr3': 'a_cdr3',
                 'beta_v':  'b_v', 'beta_j':  'b_j', 'beta_cdr3':  'b_cdr3'}

            Required. A chain whose columns are absent is left blank, so a single-chain
            frame works and only the chain present comes back stitched.
        chains: ``('alpha', 'beta')`` or ``('gamma', 'delta')`` - thimble stitches one
            receptor per call and the two cannot be mixed.
        species: passed to thimble's ``-s``.
        constant: per-chain constant regions, defaulting to :data:`DEFAULT_CONSTANT`
            so a ``stitch_pair`` and a :func:`stitch` sequence for the same row match.
            ``{'beta': 'TRBC2*01'}`` sets beta and leaves alpha to thimble's own choice;
            ``'auto'`` or ``{}`` hands both back to thimble's J-gene-aware autofill.
        leader: fills thimble's ``{locus}_leader`` column - a string for every chain, or
            a dict per chain. Left as ``None`` the column is blank, which means the V
            gene's own leader; thimble cannot be told to omit the leader. But Stitchr
            uses a *DNA* leader verbatim, so ``leader='NNN'`` gives a 3 nt placeholder
            instead of a real leader: the rest of the sequence is then exactly what
            ``stitch(no_leader=True)`` produces, recoverable by dropping 3 nt / 1 aa.
            That also rescues V genes with no leader in IMGT, which otherwise fail
            outright. Any multiple of 3 works; a length that is not (e.g. ``'N'``)
            shifts the frame and changes the whole sequence.
        extra_genes, seamless: thimble's ``-xg`` / ``-sl``.
        extra_args: anything else for the thimble command line, e.g.
            ``('-p', '/path/preferred_alleles.tsv')``.
        in_file, out_file, work_dir, keep_files: the TSVs thimble reads and writes go to
            a temp dir that is removed afterwards, unless `in_file` / `out_file` name
            either of them (in which case that one is written where asked and kept),
            `work_dir` puts the unnamed ones in a directory of your choosing, or
            `keep_files=True` keeps the temp dir and prints the two paths.
        verbose: print the command and thimble's output.

    Returns:
        A DataFrame indexed like `df` with, per chain, ``{chain}_sequence`` /
        ``{chain}_sequence_aa`` ('-' when that chain didn't stitch), plus:
            error           '-' when every requested chain stitched, else thimble's
                            message
            thimble_log     thimble's raw Warnings/Errors, kept even on success - it is
                            where the informational notes live ('defaulting to *01 ...'),
                            '-' if silent
    """
    chains = tuple(chains)
    if chains not in _THIMBLE_RECEPTORS:
        raise ValueError(f"chains must be one of {sorted(_THIMBLE_RECEPTORS)}, got {chains}.")
    loci = _THIMBLE_RECEPTORS[chains]

    # one convention across the package: prefix the field with the chain whenever a
    # function handles more than one, exactly as build_tcr_pmhc and tcrdist do
    cols = {}
    for ch in chains:
        got = {f: col_mapping.get(f'{ch}_{f}') for f in _REQUIRED_COLS} \
            if col_mapping else {}
        if all(got.get(f) for f in _REQUIRED_COLS):
            cols[ch] = got
    if not cols:
        raise ValueError(
            f"col_mapping is required and must name v / j / cdr3 for at least one of "
            f"{list(chains)}, prefixed with the chain - e.g. "
            "{'beta_v': 'v_call', 'beta_j': 'j_call', 'beta_cdr3': 'junction_aa'}"
        )

    requested = [ch for ch in chains
                 if ch in cols and all(c in df.columns for c in cols[ch].values())]
    if not requested:
        raise ValueError(
            f"None of {list(chains)} has its V/J/CDR3 columns in the DataFrame: "
            f"looked for {[c for ch in cols for c in cols[ch].values()]}"
        )
    if constant is None:
        const = {ch: DEFAULT_CONSTANT.get(ch, '') for ch in requested}
    elif isinstance(constant, str) and constant.strip().lower() == CONSTANT_AUTO:
        const = {}          # a blank C column makes thimble apply its own autofill
    elif isinstance(constant, str):
        const = {ch: constant for ch in requested}
    else:
        const = dict(constant)
    lead = ({} if leader is None else
            {ch: leader for ch in requested} if isinstance(leader, str) else dict(leader))

    # exactly thimble's expected header for this receptor - it compares it verbatim
    l1, l2 = loci
    headers = ['TCR_name', f'{l1}V', f'{l1}J', f'{l1}_CDR3', f'{l2}V', f'{l2}J', f'{l2}_CDR3',
               f'{l1}C', f'{l2}C', f'{l1}_leader', f'{l2}_leader', 'Linker', 'Link_order',
               f'{l1}_5_prime_seq', f'{l1}_3_prime_seq',
               f'{l2}_5_prime_seq', f'{l2}_3_prime_seq']

    vals = {ch: {k: _thimble_col(df, cols[ch][k] if ch in cols else '')
                 for k in _REQUIRED_COLS} for ch in chains}
    # thimble reads ',' / '%' in any field as 'stitch every combination', which renames
    # the row ('tcr0' -> 'tcr0-0', 'tcr0-1', ...) and returns more rows than went in -
    # there'd be no honest way to line that back up with df, so refuse rather than pick
    # an allele silently.
    ambiguous = {cols[ch][k]: sum(1 for v in vals[ch][k] if ',' in v or '%' in v)
                 for ch in requested for k in _REQUIRED_COLS
                 if any(',' in v or '%' in v for v in vals[ch][k])}
    if ambiguous:
        raise ValueError(
            f"Multi-valued gene calls found in {ambiguous} (counts of affected rows). "
            "thimble expands ',' and '%' into one stitched TCR per combination, so these "
            "rows can't be mapped back to the frame - pick a single V/J/CDR3 per row first."
        )
    a, b = vals[chains[0]], vals[chains[1]]
    # positional names, so a non-unique / oddly typed df index can't collide or break the TSV
    names = [f'tcr{i}' for i in range(len(df))]

    tmp = work_dir or tempfile.mkdtemp(prefix='tcrkit_thimble_')
    in_tsv = Path(in_file) if in_file else Path(tmp) / 'thimble_in.tsv'
    out_tsv = Path(out_file) if out_file else Path(tmp) / 'thimble_out.tsv'
    for p in (in_tsv, out_tsv):
        p.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(in_tsv, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(headers)
            for i, name in enumerate(names):
                w.writerow([name, a['v'][i], a['j'][i], a['cdr3'][i],
                            b['v'][i], b['j'][i], b['cdr3'][i],
                            const.get(chains[0], ''), const.get(chains[1], ''),
                            lead.get(chains[0], ''), lead.get(chains[1], '')] + [''] * 6)

        cmd = [_thimble_bin(thimble_bin), '-in', str(in_tsv), '-o', str(out_tsv),
               '-s', _thimble_species(species), '-r', chains[0][0] + chains[1][0]]
        cmd += (['-xg'] if extra_genes else []) + (['-sl'] if seamless else []) \
            + list(extra_args)
        if verbose:
            print(f"Running thimble command: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if verbose:
            print(f"Thimble output: {proc.stdout}")
        if proc.returncode != 0 or not out_tsv.exists():
            if verbose:
                print(f"Thimble stderr: {proc.stderr}")
            raise RuntimeError(f"thimble failed (exit {proc.returncode}):\n"
                               f"{(proc.stderr or proc.stdout)[-2000:]}")

        with open(out_tsv) as f:
            stitched = {r['TCR_name']: r for r in csv.DictReader(f, delimiter='\t')}
    finally:
        if keep_files:
            print(f"thimble files kept: {in_tsv} -> {out_tsv}")
        elif not work_dir:
            # only ever removes the temp dir, so in_file / out_file elsewhere are left alone
            shutil.rmtree(tmp, ignore_errors=True)

    out = {}
    for ch, locus in zip(chains, loci):
        for suffix, field in (('sequence', f'{locus}_nt'), ('sequence_aa', f'{locus}_aa')):
            out[f'{ch}_{suffix}'] = [(stitched.get(n, {}).get(field) or '').strip() or '-'
                                     for n in names]
    # thimble reports notes and failures in the one field ('[None]' when it has nothing
    # to say), so 'did every requested chain come back' is what decides `error`
    log = [(stitched.get(n, {}).get('Warnings/Errors') or '').strip() for n in names]
    log = ['' if m.lower() in ('[none]', 'none') else m for m in log]
    errors = []
    for i, name in enumerate(names):
        if all(out[f'{ch}_sequence'][i] != '-' for ch in requested):
            errors.append('-')
        elif log[i]:
            errors.append(log[i])
        elif [ch for ch in requested if not all(vals[ch][k][i] for k in _REQUIRED_COLS)]:
            errors.append('missing v/j/cdr3 input')
        else:
            errors.append('not stitched (no message from thimble)'
                          if name in stitched else 'dropped by thimble')
    out[ERROR_COL] = errors
    out['thimble_log'] = [m or '-' for m in log]
    return pd.DataFrame(out, index=df.index)
