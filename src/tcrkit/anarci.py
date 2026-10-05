"""Annotate V / J genes and CDR regions of a protein sequence with ANARCI.

One function, :func:`anarci`. It takes a single sequence or any collection of them and
always returns a DataFrame with the CDR1/2/3 loops, their positions, the assigned
germline V/J genes, species and chain type::

    anarci(seq, allow={'A', 'B'})
    anarci([seq1, seq2], allow={'A', 'B'})
    anarci(df, col_mapping={'sequence': 'sequence_aa'}, allow={'A', 'B'})

Pass ``allow=`` for TCRs
------------------------
ANARCI matches against HMMs for every chain type it knows, and the alpha and delta loci
overlap heavily. Without restricting the search it can report a **delta** assignment for
what is really an alpha chain, giving the wrong J gene::

    anarci(seq)                    # chain_type 'D', j_gene 'TRDJ4*01'   <- wrong
    anarci(seq, allow={'A', 'B'})  # chain_type 'A', j_gene 'TRAJ29*01'  <- right

So pass ``allow={'A', 'B'}`` for alpha/beta TCRs (or ``{'G', 'D'}`` for gamma/delta)
unless you specifically want ANARCI to choose among all chain types. When ``allow`` is
given, tcrkit also makes ANARCI's domain de-duplication chain-type aware, so the wanted
domain is not discarded in favour of an overlapping hit of another type.

ANARCI normally shells out to the ``hmmscan`` binary, which conda/mamba environments
ship but a plain pip/uv environment does not. This module falls back to an in-process
``pyhmmer`` backend (:mod:`tcrkit._anarci_pyhmmer`) when no ``hmmscan`` is on PATH, so
``hmmscan`` is not a hard requirement. The choice is made once per process and can be
forced either way with ``ANARCI_USE_PYHMMER=1`` / ``=0`` in the environment, or with the
``use_pyhmmer=`` argument.

The numbering itself is unchanged: the pyhmmer backend is a faithful port of ANARCI's
own state-vector logic and produces identical output to the ``hmmscan`` path.
"""

from __future__ import annotations

import contextlib
import os
import re
import shutil
import sys
from io import StringIO
from pathlib import Path

import pandas as pd

from ._util import ERROR_COL, collect, is_empty, progress, squeeze as _squeeze

__all__ = ['anarci', 'CDR_IMGT_LOCS', 'ERROR_ROW', 'read_sequences', 'backend']


#: CDR positions in the IMGT numbering scheme.
#: https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
CDR_IMGT_LOCS = {
    'cdr1': (27, 38),
    'cdr2': (56, 65),
    'cdr3': (105, 117),
}

#: What a sequence ANARCI could not number comes back as.
ERROR_ROW = {
    'CDR1_aa_start': -1, 'CDR1_aa_end': -1, 'CDR1_aa': '-',
    'CDR2_aa_start': -1, 'CDR2_aa_end': -1, 'CDR2_aa': '-',
    'CDR3_aa_start': -1, 'CDR3_aa_end': -1, 'CDR3_aa': '-',
    'query_start': -1, 'query_end': -1, 'species': '-',
    'chain_type': '-', 'v_gene': '-', 'v_identity': -1.0,
    'j_gene': '-', 'j_identity': -1.0, ERROR_COL: 'ANARCI error',
}

_DROP_COLS = ['id', 'description', 'evalue', 'bitscore', 'bias', 'scheme']

_backend = None


# ---------------------------------------------------------------------------
# Backend selection
# ---------------------------------------------------------------------------

def _ensure_backend(hmmerpath: str = '', use_pyhmmer: bool | None = None) -> bool:
    """Patch in the pyhmmer backend if ANARCI needs it. Returns True if patched.

    conda/mamba envs ship the ``hmmscan`` binary, so stock ANARCI works there and
    neither the shim nor its pyhmmer dependency is touched. A plain uv/pip env has no
    binary, so the in-process pyhmmer backend is patched in. Decided once per process.
    """
    global _backend
    if _backend is None:
        if use_pyhmmer is None:
            forced = os.environ.get('ANARCI_USE_PYHMMER')
            if forced is not None:
                need = forced.strip().lower() in ('1', 'true', 'yes')
            else:
                # Mirror how anarci.anarci.run_hmmer resolves the binary.
                hmmscan = (os.path.join(hmmerpath, 'hmmscan') if hmmerpath
                           else shutil.which('hmmscan'))
                need = not (hmmscan and os.access(hmmscan, os.X_OK))
        else:
            need = bool(use_pyhmmer)
        if need:
            from ._anarci_pyhmmer import patch_anarci
            patch_anarci()
        _backend = need
    return _backend


def backend() -> str:
    """Which HMMER backend is in use: 'pyhmmer', 'hmmscan', or 'undecided'.

    'undecided' until the first numbering call, since the choice is made lazily.
    """
    return {None: 'undecided', True: 'pyhmmer', False: 'hmmscan'}[_backend]


@contextlib.contextmanager
def _split_chain_types():
    """Make ANARCI's domain dedup also require a matching chain type.

    Without this, restricting to one chain with ``allow=`` can discard the wanted
    domain in favour of an overlapping hit of another chain type.
    """
    import anarci as _anarci_pkg  # noqa: F401  (registers anarci.anarci in sys.modules)
    _a = sys.modules['anarci.anarci']
    orig = _a._domains_are_same
    _a._domains_are_same = (
        lambda x, y: x.hit_id.split('_')[1] == y.hit_id.split('_')[1] and orig(x, y))
    try:
        yield
    finally:
        _a._domains_are_same = orig


# ---------------------------------------------------------------------------
# Turning ANARCI's nested output into frames
# ---------------------------------------------------------------------------

def _numbered_to_df(anarci_output):
    """The IMGT-numbered residues as a frame, one column per IMGT position."""
    sequences, numbered, _alignment_details, _hit_tables = anarci_output

    rows, indices = [], []
    for r in range(len(sequences)):
        indices.append(sequences[r][0])
        if numbered[r] is None:
            rows.append({})
            continue
        try:
            rows.append({''.join(str(j) for j in i[0]).strip(): i[1]
                         for i in numbered[r][0][0]})
        except Exception:
            rows.append({})
    return pd.DataFrame(rows, index=indices).fillna('-')


def _alignment_details_to_df(anarci_output):
    """The per-domain alignment details, with the germline assignments flattened out."""
    sequences, _numbered, alignment_details, _hit_tables = anarci_output

    indices = [sequences[r][0] for r in range(len(sequences))]
    rows = [{} if alignment_details[r] is None else alignment_details[r][0]
            for r in range(len(alignment_details))]
    df = pd.DataFrame(rows, index=indices).fillna('-')

    def germline_info(row):
        try:
            g = row['germlines']
            return pd.Series({
                'v_gene': g['v_gene'][0][1] if 'v_gene' in g else '-',
                'v_identity': g['v_gene'][1] if 'v_gene' in g else '-',
                'j_gene': g['j_gene'][0][1] if 'j_gene' in g else '-',
                'j_identity': g['j_gene'][1] if 'j_gene' in g else '-',
            })
        except Exception:
            return pd.Series({'v_gene': '-', 'v_identity': '-',
                              'j_gene': '-', 'j_identity': '-'})

    if 'germlines' in df.columns:
        df = df.join(df.apply(germline_info, axis=1)).drop(columns=['germlines'],
                                                           errors='ignore')
    return df


def _extract_cdrs(df_numbered, df_alignment_details=None, cdr_locs=None):
    """CDR sequences and their positions in the query, from the numbered frame.

    Positions are 1-based indices into the input sequence, offset by the domain's
    ``query_start`` so they refer to the sequence as it was passed in.
    """
    cdr_locs = CDR_IMGT_LOCS if cdr_locs is None else cdr_locs
    adj_ind = (df_numbered != '-').cumsum(axis=1)

    if (df_alignment_details is not None
            and 'query_start' in df_alignment_details.columns):
        adj_ind = adj_ind.add(df_alignment_details['query_start'], axis=0)

    out = pd.DataFrame(index=df_numbered.index)
    for cdr, (start, end) in cdr_locs.items():
        cdr_columns = df_numbered.loc[:, str(start):str(end)]
        cdr_sequences = cdr_columns.where(cdr_columns != '-', '').agg(''.join, axis=1)
        if cdr == 'junction':
            out[f'{cdr}_aa'] = cdr_sequences
        else:
            out[f'{cdr.upper()}_aa'] = cdr_sequences
            out[f'{cdr.upper()}_aa_start'] = adj_ind.loc[:, str(start)]
            out[f'{cdr.upper()}_aa_end'] = adj_ind.loc[:, str(end)]

    if 'Id' in df_numbered.columns:
        out.insert(0, 'Id', df_numbered['Id'])
    return out


def _prepare(anarci_output, include_imgt_index=False):
    """One sequence's ANARCI output as a single-row frame."""
    df_align = _alignment_details_to_df(anarci_output)
    df_num = _numbered_to_df(anarci_output)
    if len(df_num) == 0:
        raise ValueError('No data found in the numbered DataFrame.')
    if not df_num.index.is_unique:
        raise ValueError('Sequence ids must be unique.')

    df = pd.concat([_extract_cdrs(df_num, df_align), df_align], axis=1)
    df = df.drop(columns=[c for c in _DROP_COLS if c in df.columns])
    df[ERROR_COL] = '-'

    if include_imgt_index:
        df_num = df_num.loc[:, df_num.iloc[0] != '-']
        row_dicts = df_num.to_dict(orient='index')
        df['imgt_dict'] = df.index.map(row_dicts.get)
    return df


def _expand_imgt_index(df):
    """Turn the ``imgt_dict`` column into one column per IMGT position."""
    if 'imgt_dict' not in df.columns:
        return df

    def order(col):
        m = re.match(r'^(\d+)(\D*)$', str(col))
        return (int(m.group(1)), m.group(2)) if m else (10 ** 6, str(col))

    dicts = [d if isinstance(d, dict) else {} for d in df['imgt_dict']]
    cols = sorted({c for d in dicts for c in d}, key=order)
    wide = pd.DataFrame(dicts, index=df.index).reindex(columns=cols)
    return pd.concat([df.drop(columns=['imgt_dict']), wide], axis=1)


def separate_imgt_index(df):
    """Split trailing IMGT position columns into ``_position`` / ``_suffix`` header rows."""
    nums, sufs = {}, {}
    for c in df.columns:
        m = re.match(r'^(\d+)(\D?)$', str(c))
        nums[c] = int(m.group(1)) if m else -1
        sufs[c] = (m.group(2) or ' ') if m else ' '
    position = pd.DataFrame([nums], index=['_position'])
    suffix = pd.DataFrame([sufs], index=['_suffix'])
    return pd.concat([position, suffix, df], axis=0)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _number(items, progress_bar=True, include_imgt_index=False, separate_index=False,
            mute=True, raise_error=False, use_pyhmmer=None, **anarci_kwargs):
    """IMGT-number ``[(id, sequence), ...]``. The multi-row core of :func:`anarci`.

    Kept separate so the public function owns input detection and frame assembly while
    this owns the ANARCI call and its per-sequence error handling.
    """
    import anarci as _anarci_pkg

    @contextlib.contextmanager
    def suppress():
        if mute:
            with open(os.devnull, 'w') as devnull, \
                    contextlib.redirect_stdout(devnull), \
                    contextlib.redirect_stderr(devnull):
                yield
        else:
            yield

    kwargs = {
        'scheme': 'imgt', 'output': False, 'outfile': False, 'csv': False,
        'ncpu': 1, 'hmmerpath': str(Path(sys.executable).parent),
        'assign_germline': True,
    }
    kwargs.update(anarci_kwargs)
    # ANARCI turns on chain-type restriction from the mere presence of `allow`, and
    # then fails on a None. Treat 'no restriction' as 'not passed'.
    if kwargs.get('allow', ()) is None:
        kwargs.pop('allow')

    _ensure_backend(kwargs.get('hmmerpath', ''), use_pyhmmer=use_pyhmmer)

    if not items:
        return pd.DataFrame(columns=list(ERROR_ROW))
    it = progress(items, desc='ANARCI', total=len(items),
                  enabled=progress_bar and len(items) > 5)

    results = []
    for sample in it:
        sid = sample[0]
        try:
            if is_empty(sample[1]):
                raise ValueError('empty sequence')
            with suppress(), (_split_chain_types() if 'allow' in kwargs
                             else contextlib.nullcontext()):
                out = _anarci_pkg.run_anarci([(sid, sample[1])], **kwargs)
                df_sample = _prepare(out, include_imgt_index=include_imgt_index)
        except Exception as e:
            if raise_error:
                raise
            df_sample = pd.DataFrame({**ERROR_ROW,
                                      ERROR_COL: f'{type(e).__name__}: {e}'
                                      if str(e) else ERROR_ROW[ERROR_COL]},
                                     index=[sid])
        results.append(df_sample)

    df = pd.concat(results, ignore_index=False)
    if separate_index:
        df = separate_imgt_index(_expand_imgt_index(df))
    return df.fillna('-')

def anarci(sequences, col_mapping=None, progress_bar=True, include_imgt_index=False,
           separate_index=False, mute=True, raise_error=False, use_pyhmmer=None,
           squeeze: bool = False, **anarci_kwargs) -> pd.DataFrame:
    """IMGT-number protein sequences. Returns one row per input sequence.

    >>> out = anarci('MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGTTAKGE'
    ...              'VPDGYNVSRLKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL')
    >>> out[['chain_type', 'CDR3_aa']].to_dict('records')
    [{'chain_type': 'B', 'CDR3_aa': 'ASSFTDTQY'}]

    For TCRs, pass ``allow`` so ANARCI does not pick an overlapping locus:

    >>> out = anarci('AQKITQTQPGMFVQEKEAVTLDCTYDTSDPSYGLFWYKQPSSGEMIFLIYQGSYDQQNATEG'
    ...              'RYSLNFQKARKSANLVISASQLGDSAMYFCAMSQLNSGNTPLVFGKGTRLSVIA',
    ...              allow={'A', 'B'})
    >>> out[['chain_type', 'j_gene']].to_dict('records')
    [{'chain_type': 'A', 'j_gene': 'TRAJ29*01'}]

    Args:
        sequences: a single amino-acid sequence, or any collection of them - a list,
            Series, dict of ``id -> sequence``, list of ``(id, sequence)`` pairs, list
            of records, or a DataFrame. Whatever ids the input carries become the index.
        col_mapping: where the sequence lives when the input is a frame or records,
            e.g. ``{'sequence': 'sequence_aa'}``.
        progress_bar: show a tqdm bar (only for more than 5 sequences).
        include_imgt_index: add an ``imgt_dict`` column mapping IMGT position -> residue.
        separate_index: spread the numbering into one column per IMGT position
            (implies ``include_imgt_index``) and prepend ``_position`` / ``_suffix``
            header rows carrying each column's position as an int and its insertion
            letter, so the IMGT axis can be sliced numerically. Returns the numbering
            frame alone, without the input echo, since those two header rows do not
            correspond to input sequences.
        mute: keep ANARCI's own stdout/stderr chatter quiet.
        raise_error: re-raise instead of returning an error row for a sequence that
            cannot be numbered.
        use_pyhmmer: force the in-process pyhmmer backend on (True) or off (False).
            ``None`` (default) auto-detects. Only affects the process's first call.
        squeeze: return the plain result instead of a frame - a dict for one sequence,
            a DataFrame without the ``sequence`` echo for several.
        **anarci_kwargs: forwarded to ``anarci.run_anarci``. The important one is
            ``allow``: pass ``allow={'A', 'B'}`` for alpha/beta TCRs (or ``{'G', 'D'}``
            for gamma/delta), otherwise an alpha chain can come back assigned to the
            overlapping delta locus with a TRDJ gene - see the module docstring. Also
            ``scheme='imgt'``, ``assign_germline=True``, ``ncpu=``.

    Returns:
        A DataFrame indexed by sequence id, starting with the ``sequence`` you passed,
        then ``CDR{1,2,3}_aa`` and their ``_start`` / ``_end`` positions (1-based into
        that sequence), ``v_gene`` / ``j_gene`` and their identities, ``species``,
        ``chain_type``, ``query_start`` / ``query_end``, and ``error`` ('-' when
        the sequence numbered). Sequences that failed get :data:`ERROR_ROW`, so the
        frame always has one row per input.
    """
    inp = collect(('sequence',), (sequences,), col_mapping=col_mapping)
    items = list(zip(inp.index, inp['sequence']))

    out = _number(items, progress_bar=progress_bar,
                  include_imgt_index=include_imgt_index or separate_index,
                  separate_index=separate_index, mute=mute, raise_error=raise_error,
                  use_pyhmmer=use_pyhmmer, **anarci_kwargs)
    if separate_index:
        return out
    out = out.drop(columns=['query_name'], errors='ignore')
    res = pd.concat([inp.reindex(out.index), out], axis=1)
    return _squeeze(res, ('sequence',)) if squeeze else res


def read_sequences(path):
    """Read one sequence per line from a plain-text or FASTA file -> [(id, seq), ...].

    FASTA headers are used as ids where present; unlabelled lines get 'id_1', 'id_2', ...
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'{path} not found')
    seqs: list[tuple[str, str]] = []
    label = None
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                label = line[1:].strip() or None
                continue
            seqs.append((label or f'id_{len(seqs) + 1}', line))
            label = None
    return seqs
