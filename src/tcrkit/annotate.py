"""What a protein sequence is - TCR, MHC, beta2m or peptide - and its annotation.

:func:`annotate` runs :func:`tcrkit.anarci` and :func:`tcrkit.annotate_mhc` on every
sequence and reports whichever recognised it, so you do not have to know the chain type
first::

    annotate(seq)['chain_type']                # 'TCR_beta', 'MHC_I', 'B2M', 'peptide', ...
    annotate({'a': tcr_a, 'm': hla_a2})        # a mixed collection, one row each
    annotate(seq, per_residue=True)            # one row per residue, any chain type

For the full output of one kind, call :func:`tcrkit.anarci` or :func:`tcrkit.annotate_mhc`.
"""

from __future__ import annotations

import pandas as pd

from ._util import ERROR_COL, collect, squeeze as _squeeze

__all__ = ['annotate', 'CHAIN_TYPES', 'PEPTIDE_MAX_LENGTH']

#: The values of ``chain_type``.
CHAIN_TYPES = ('TCR_alpha', 'TCR_beta', 'TCR_gamma', 'TCR_delta', 'antibody',
               'MHC_I', 'MHC_II_alpha', 'MHC_II_beta', 'B2M', 'peptide', 'unknown')

#: A sequence no tool recognises is called a peptide up to this length.
PEPTIDE_MAX_LENGTH = 30

_ANARCI_TYPE = {'A': 'TCR_alpha', 'B': 'TCR_beta', 'G': 'TCR_gamma', 'D': 'TCR_delta',
                'H': 'antibody', 'K': 'antibody', 'L': 'antibody'}
_MHC_TYPE = {'heavy': 'MHC_I', 'alpha': 'MHC_II_alpha', 'beta': 'MHC_II_beta', 'B2M': 'B2M'}

_COLUMNS = ['chain_type', 'species', 'v_gene', 'j_gene', 'CDR1_aa', 'CDR2_aa', 'CDR3_aa',
            'gene', 'allele', 'allele_identity', 'mhc_class', 'domains', 'helix', ERROR_COL]


def _ok(value) -> bool:
    return isinstance(value, str) and value not in ('', '-')


def _classify(seq, tcr_type, mhc_gene, why_not_tcr, why_not_mhc):
    """``(chain_type, error)`` for one sequence, from what each tool found."""
    from .mhc import mhc_chain

    tcr = _ANARCI_TYPE.get(tcr_type) if _ok(tcr_type) else None
    mhc = _MHC_TYPE.get(mhc_chain(mhc_gene)) if _ok(mhc_gene) else None
    if tcr and mhc:
        return 'unknown', f'recognised as both {tcr} and {mhc}'
    if tcr or mhc:
        return tcr or mhc, None
    if not seq:
        return 'unknown', 'empty sequence'
    if len(seq) <= PEPTIDE_MAX_LENGTH:
        return 'peptide', None
    reasons = '; '.join(r for r in (why_not_tcr, why_not_mhc) if _ok(r))
    return 'unknown', reasons or 'not a TCR, MHC or B2M chain'


def annotate(sequences, per_residue: bool = False, helix: str = 'imgt',
             col_mapping=None, squeeze: bool = False, progress_bar: bool = True,
             **anarci_kwargs) -> pd.DataFrame:
    """What each sequence is, and its annotation, without saying the chain type first.

    Every sequence goes through :func:`tcrkit.anarci` (TCR and antibody V domains) and
    :func:`tcrkit.annotate_mhc` (MHC chains and B2M); ``chain_type`` says which
    recognised it. Neither, and at most :data:`PEPTIDE_MAX_LENGTH` residues: a peptide.

    >>> out = annotate({'beta': 'MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGT'
    ...                         'TAKGEVPDGYNVSRLKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL',
    ...                 'pep': 'GILGFVFTL'}, progress_bar=False)
    >>> out['chain_type'].tolist()
    ['TCR_beta', 'peptide']

    Args:
        sequences: one amino-acid sequence or any collection (list, dict, DataFrame...).
        per_residue: one row per residue instead, for every chain type: ``position``,
            ``aa``, ``chain_type``, ``domain``, ``imgt``, ``region``, and per kind
            ``gene_segment`` (TCR) or ``gene`` / ``allele`` / ``helix`` / ``topology``
            (MHC). Works with :func:`tcrkit.print_residues`.
        helix: how the MHC groove helices are found, as in :func:`tcrkit.annotate_mhc`.
        col_mapping: where the sequence lives in a frame, e.g. ``{'sequence': 'seq'}``.
        squeeze: return the plain result instead of a frame.
        progress_bar: show progress for the ANARCI pass.
        **anarci_kwargs: passed to :func:`tcrkit.anarci`. ANARCI picks among every chain
            type it knows, so an alpha chain can come back as ``TCR_delta``; pass
            ``allow={'A', 'B'}`` if you are looking for alpha/beta chains.

    Returns:
        Indexed like the input: ``sequence``, ``chain_type`` (one of
        :data:`CHAIN_TYPES`), ``species``, then the TCR columns ``v_gene``, ``j_gene``,
        ``CDR1_aa`` .. ``CDR3_aa``, the MHC columns ``gene``, ``allele``,
        ``allele_identity``, ``mhc_class``, ``helix``, and for both ``domains``
        (``'V, C'`` or ``'G-ALPHA1, G-ALPHA2, C-LIKE'``) and ``error``. Columns that do
        not apply to a chain are empty.
    """
    from .anarci import anarci
    from .mhc import annotate_mhc, mhc_species

    inp = collect(('sequence',), (sequences,), col_mapping=col_mapping)
    seqs = inp['sequence'].map(lambda s: '' if not isinstance(s, str) else s.strip().upper())
    named = pd.Series(seqs.values, index=inp.index)

    if per_residue:
        tcr = anarci(named, per_residue=True, domains=True, progress_bar=progress_bar,
                     **anarci_kwargs)
        mhc = annotate_mhc(named, per_residue=True, domains=True, helix=helix)
        parts = []
        for idx, seq in named.items():
            t, m = tcr[tcr.index == idx], mhc[mhc.index == idx]
            ctype, _ = _classify(seq, t['chain_type'].iloc[0] if len(t) else None,
                                 m['gene'].iloc[0] if len(m) else None, None, None)
            if ctype.startswith('TCR') or ctype == 'antibody':
                part = t
            elif ctype.startswith('MHC') or ctype == 'B2M':
                part = m
            else:
                part = pd.DataFrame({'position': range(1, len(seq) + 1), 'aa': list(seq)},
                                    index=[idx] * len(seq))
            parts.append(part.assign(chain_type=ctype))
        res = pd.concat(parts) if parts else pd.DataFrame(columns=['position', 'aa'])
        lead = ['position', 'aa', 'chain_type', 'domain', 'imgt', 'region']
        res = res[[c for c in lead if c in res] + [c for c in res if c not in lead]]
        return _squeeze(res) if squeeze else res

    tcr = anarci(named, domains=True, progress_bar=progress_bar, **anarci_kwargs)
    mhc = annotate_mhc(named, domains=True, helix=helix)
    rows = []
    for idx, seq in named.items():
        t, m = tcr.loc[idx], mhc.loc[idx]
        ctype, err = _classify(seq, t['chain_type'], m['gene'], t[ERROR_COL], m[ERROR_COL])
        row = dict.fromkeys(_COLUMNS)
        if ctype.startswith('TCR') or ctype == 'antibody':
            row.update({c: t[c] for c in ('v_gene', 'j_gene', 'CDR1_aa', 'CDR2_aa',
                                          'CDR3_aa', 'domains')},
                       species=t['germline_species'])
        elif ctype.startswith('MHC') or ctype == 'B2M':
            row.update({c: m[c] for c in ('gene', 'allele', 'allele_identity',
                                          'mhc_class', 'domains', 'helix')},
                       species=mhc_species(m['allele']))
        row.update(chain_type=ctype, error=err)
        rows.append(row)
    out = pd.DataFrame(rows, index=inp.index, columns=_COLUMNS)
    res = pd.concat([inp[['sequence']], out], axis=1)
    return _squeeze(res, ('sequence',)) if squeeze else res
