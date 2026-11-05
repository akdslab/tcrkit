"""The IMGT position ranges tcrkit labels residues with, and where each comes from."""

from __future__ import annotations

import pandas as pd

from .anarci import IMGT_REGIONS
from .gdomain import GDOMAIN_REGIONS, GDOMAINS, HELIX

__all__ = ['imgt_regions', 'IMGT_SOURCES']

#: Where each set of ranges comes from.
IMGT_SOURCES = {
    'tcr': 'IMGT V-DOMAIN FR/CDR definition: '
           'https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html',
    'mhc': 'IMGT G-DOMAIN: '
           'https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html',
    'helix': "tcrkit helix='imgt': the groove-helix residues TCR docking angles are "
             'measured on; on class I the same as the MHC Motif Atlas '
             '(http://mhcmotifatlas.org)',
}

_SCHEMES = ('tcr', 'mhc')

#: ``(mhc_class, chain)`` of each G domain: class I has two helices on its heavy chain,
#: class II one on each chain.
_DOMAIN_CHAIN = {'G-ALPHA1': ('I', 'heavy'), 'G-ALPHA2': ('I', 'heavy'),
                 'G-ALPHA': ('II', 'alpha'), 'G-BETA': ('II', 'beta')}


def imgt_regions(scheme: str | None = None, all_regions: bool = False) -> pd.DataFrame:
    """The IMGT ranges behind tcrkit's ``region`` and ``helix`` columns, with sources.

    Args:
        scheme: ``'tcr'`` (V domain: FR1-4, CDR1-3, as :func:`tcrkit.anarci` labels
            them), ``'mhc'`` (G domains: strands, turns, helix, as
            :func:`tcrkit.annotate_mhc` labels them), or ``None`` for both.
        all_regions: every region. By default only the important ones: the CDRs and
            the groove helices (the ``helix`` column).

    Returns:
        One row per range: ``scheme``, ``mhc_class`` and ``chain`` (``'-'`` for TCR,
        whose CDRs are at the same positions on every chain), ``domain``, ``region``, ``column`` (the output
        column that carries it: ``region`` or ``helix``), ``imgt_start`` / ``imgt_end``
        (inclusive; insertions such as ``'111A'`` or ``'61A'`` fall in the range of their
        number), ``skipped`` (positions left out) and ``source``.

        The helix positions are fixed per G domain: class I has two (G-ALPHA1 and
        G-ALPHA2, both on the heavy chain), class II one per chain (G-ALPHA, G-BETA).
        For MHC, ``region == 'helix'`` (IMGT 50-92) is IMGT's whole helix region;
        the ``helix`` column rows are the narrower groove-helix residues, which are the
        reference for the helices.
    """
    if scheme is not None and scheme not in _SCHEMES:
        raise ValueError(f'scheme must be one of {_SCHEMES} or None, got {scheme!r}')
    rows = []
    if scheme in (None, 'tcr'):
        rows += [('tcr', '-', '-', 'V-DOMAIN', r, 'region', lo, hi, '', IMGT_SOURCES['tcr'])
                 for r, (lo, hi) in IMGT_REGIONS.items()]
    if scheme in (None, 'mhc'):
        for d in GDOMAINS:
            cls, chain = _DOMAIN_CHAIN[d]
            rows += [('mhc', cls, chain, d, r, 'region', lo, hi, '', IMGT_SOURCES['mhc'])
                     for r, (lo, hi) in GDOMAIN_REGIONS.items()]
            lo, hi, skipped = HELIX[d]
            rows.append(('mhc', cls, chain, d, 'helix', 'helix', lo, hi, ', '.join(skipped),
                         IMGT_SOURCES['helix']))
    out = pd.DataFrame(rows, columns=['scheme', 'mhc_class', 'chain', 'domain', 'region', 'column',
                                      'imgt_start', 'imgt_end', 'skipped', 'source'])
    if not all_regions:
        out = out[out['region'].str.startswith('CDR') | out['column'].eq('helix')]
    return out.reset_index(drop=True)
