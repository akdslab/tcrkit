"""IMGT numbering of constant (C) and C-like domains, such as B2M.

The counterpart of :mod:`tcrkit.gdomain` for IMGT's other domain type outside the V
domain: the C-DOMAIN of IG / TR constant regions and the C-LIKE-DOMAIN of other
proteins (B2M, MHC alpha3 / alpha2 / beta2). Both use one numbering, IMGT positions
1-128 with 46-76 absent, and insertions written ``'45.1'``, ``'84.1'`` ... ``'85.1'``,
``'1.8'`` ... ``'1.1'``::

    number_cdomain(seq, genes=['B2M'])     # one row per residue

As with G domains, the sequence is numbered by local alignment to its closest IMGT
entry: the human and mouse rows of IMGT's C-DOMAIN / C-LIKE-DOMAIN protein display,
bundled as ``data/imgt_cdomains.csv``.

IMGT data: Lefranc et al., Dev. Comp. Immunol. 29, 185-203 (2005), and
https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/C-DOMAIN/CLdomains.html,
under CC BY 4.0.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .gdomain import _MIN_SCORE, _aligner, _number_one

__all__ = ['number_cdomain', 'cdomain_region', 'imgt_cdomain_templates', 'CDOMAIN_REGIONS',
           'C_OFFSET', 'offset_label', 'MHC_CLIKE_GENES']

DATA = Path(__file__).resolve().parent / 'data'

#: IMGT C-domain strands and loops, inclusive position ranges.
CDOMAIN_REGIONS = {'A': (1, 15), 'B': (16, 26), 'BC': (27, 38), 'C': (39, 45),
                   'D': (77, 84), 'E': (85, 96), 'F': (97, 104), 'FG': (105, 117),
                   'G': (118, 128)}

#: Insertions after (``15.1``) or before (``85.1``) these positions are the turns.
_TURNS = {15: 'AB', 45: 'CD', 84: 'DE', 85: 'DE', 96: 'EF'}

#: Positions whose insertions come before them, counting down (``1.8`` ... ``1.1``).
_BEFORE = {1, 85, 112}

_LABEL = re.compile(r'^(\d+)(?:\.(\d+))?$')

#: Added to C / C-LIKE positions when they share a chain's ``imgt`` column with a V or G
#: domain, so the labels never repeat: position 23 becomes ``'1023'``, 45.1 ``'1045.1'``.
C_OFFSET = 1000

#: The IMGT C-LIKE templates of each MHC chain kind's membrane-proximal domain: class I
#: alpha3, class II alpha2 / beta2, and B2M.
MHC_CLIKE_GENES = {
    'heavy': ('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G', 'H2-D1', 'H2-L',
              'H2-M5', 'H2-Q7', 'MH1-K1'),
    'alpha': ('HLA-DRA', 'HLA-DQA1', 'HLA-DMA', 'HLA-DOA', 'H2-EA', 'H2-DMA', 'MH2-AA'),
    'beta': ('HLA-DRB1', 'HLA-DQB1', 'HLA-DPB1', 'HLA-DMB', 'HLA-DOB', 'H2-AB', 'H2-EB'),
    'B2M': ('B2M',),
}


def offset_label(label: str) -> str:
    """`label` plus :data:`C_OFFSET`: ``offset_label('45.1') == '1045.1'``."""
    number, _, decimal = label.partition('.')
    return str(int(number) + C_OFFSET) + (f'.{decimal}' if decimal else '')


@lru_cache(maxsize=1)
def imgt_cdomain_templates() -> pd.DataFrame:
    """IMGT's numbered C / C-like domain sequences: one row per entry.

    Columns ``species``, ``gene``, ``accession``, ``domain`` (the domain number in the
    chain, e.g. ``3`` for an MHC alpha3, or ``CH2``), ``sequence`` and ``imgt``, a list
    with one IMGT label per residue. Rebuilt by
    ``python -m tcrkit._build.fetch_imgt_cdomains``; shipped with the package.
    """
    df = pd.read_csv(DATA / 'imgt_cdomains.csv', keep_default_na=False, dtype=str)
    df['imgt'] = df['imgt'].str.split()
    return df


def _slot_key(label: str) -> tuple[int, int]:
    number, decimal = _LABEL.match(label).groups()
    n, d = int(number), int(decimal or 0)
    return (n, -d) if n in _BEFORE else (n, d)


@lru_cache(maxsize=1)
def _slots() -> tuple[str, ...]:
    """Every position any template uses, in IMGT order."""
    labels = {lab for labs in imgt_cdomain_templates()['imgt'] for lab in labs}
    return tuple(sorted(labels, key=_slot_key))


def cdomain_region(label: str) -> str:
    """The IMGT C-domain strand or loop of a position: ``'A'``, ``'AB'``, ... ``'G'``, or ``''``.

    >>> cdomain_region('23'), cdomain_region('45.2'), cdomain_region('85.1')
    ('B', 'CD', 'DE')
    """
    m = _LABEL.match(str(label))
    if m is None:
        return ''
    n = int(m.group(1))
    if m.group(2) and n in _TURNS:
        return _TURNS[n]
    return next((name for name, (lo, hi) in CDOMAIN_REGIONS.items() if lo <= n <= hi), '')


def number_cdomain(sequence: str, genes=None) -> pd.DataFrame:
    """IMGT numbering of the C or C-like domain in one sequence.

    Args:
        sequence: one-letter amino-acid sequence.
        genes: only number on templates of these IMGT genes (e.g. ``['B2M']``); default
            all. The domain is kept if its best template scores at least 0.4 of a
            perfect match.

    Returns:
        One row per residue: ``aa``, ``imgt`` (``'-'`` outside the domain), ``region``
        (:func:`cdomain_region`) and ``template`` (the IMGT entry it was numbered on).
        ``attrs['score']`` holds the best template's score.
    """
    out = _number_cdomain(sequence, None if genes is None else tuple(genes))
    return out.copy(deep=True).__finalize__(out)


@lru_cache(maxsize=1024)
def _number_cdomain(sequence: str, genes) -> pd.DataFrame:
    templates = imgt_cdomain_templates()
    if genes is not None:
        templates = templates[templates['gene'].isin(genes)]
    if templates.empty:
        raise ValueError(f'no IMGT C-domain template for genes {list(genes)}')
    aligner = _aligner()
    scored = [(aligner.score(sequence, t.sequence) / aligner.score(t.sequence, t.sequence), t)
              for t in templates.itertuples()]
    best, template = max(scored, key=lambda s: s[0])
    out = pd.DataFrame({'aa': list(sequence), 'imgt': '-', 'region': '', 'template': ''})
    out.attrs['score'] = round(best, 3)
    if best < _MIN_SCORE:
        return out
    numbering = _number_one(sequence, template, slots=_slots())
    pos = np.fromiter(numbering, int)
    out.loc[pos, 'imgt'] = [lab for lab, _ in numbering.values()]
    out.loc[pos, 'region'] = [cdomain_region(lab) for lab, _ in numbering.values()]
    out.loc[pos, 'template'] = f'{template.species} {template.gene} {template.accession}'
    return out
