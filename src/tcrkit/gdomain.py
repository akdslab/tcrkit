"""IMGT numbering of MHC groove domains (G-DOMAINs), and the groove helices.

The numbering behind :func:`tcrkit.annotate_mhc`. :func:`number_gdomains` numbers each
MHC groove domain of one sequence on IMGT's G-DOMAIN scheme: strands A-D (positions
1-49) and the helix (50-92B), with insertions such as ``'54A'``. Class I heavy chains
have two G domains, G-ALPHA1 and G-ALPHA2; class II has G-ALPHA on the alpha chain and
G-BETA on the beta chain. An IMGT position means the same place in every allele, so
this is the MHC counterpart of ANARCI's numbering of the TCR V domains::

    number_gdomains(seq)                       # one row per residue
    number_gdomains(seq, ['G-BETA'])           # only look for some domains

The ``helix`` column marks the groove-helix residues that TCR docking angles are
measured on, by one of :data:`HELIX_METHODS`:

- ``'imgt'`` (default): fixed IMGT positions (:data:`HELIX`), so the same residues
  in every allele: HLA-A2 50-86 and 140-176, HLA-DR alpha 46-78, beta 54-64 and
  67-91. They are part of IMGT's helix region (50-92B), not all of it.
- ``'tpm-bench'``: the helix ranges of the closest MHC Motif Atlas gene
  (http://mhcmotifatlas.org), carried over by sequence alignment, as the tpm-bench
  benchmark does. The same residues as ``'imgt'`` on class I; on class II the ends
  differ by a few residues, gene by gene.

Each domain is numbered by local alignment to its closest IMGT entry (human, mouse
and rat entries of IMGT's G-DOMAIN protein display, bundled as
``data/imgt_gdomains.csv``). Positions IMGT prints only as digits are written
``'-10'`` .. ``'-1'`` (before position 1) and ``'49.1'`` .. ``'49.7'`` (between
strand D and the helix). Residues outside the G domains (alpha3 / beta2, leaders,
tails) are ``'-'``.

IMGT data: Lefranc et al., Dev. Comp. Immunol. 29, 917-938 (2005), and
https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html,
under CC BY 4.0.
"""

from __future__ import annotations

import ast
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ['number_gdomains', 'gdomain_region', 'is_helix', 'tpm_bench_helix',
           'imgt_gdomain_templates', 'helix_genes', 'GDOMAINS', 'GDOMAIN_REGIONS', 'HELIX',
           'HELIX_METHODS']

DATA = Path(__file__).resolve().parent / 'data'

#: IMGT's MHC groove domains.
GDOMAINS = ('G-ALPHA1', 'G-ALPHA2', 'G-ALPHA', 'G-BETA')

#: IMGT G-domain regions, inclusive position ranges (insertions such as 61A are in
#: the range of their number). Positions 49.1-49.7 and those before 1 are in none.
GDOMAIN_REGIONS = {'A': (1, 14), 'AB': (15, 17), 'B': (18, 28), 'BC': (29, 30),
                   'C': (31, 38), 'CD': (39, 41), 'D': (42, 49), 'helix': (50, 92)}

#: The groove-helix residues of the ``'imgt'`` method, in IMGT positions:
#: ``{domain: (first, last, skipped)}``. G-BETA skips 61B and 62 (HLA-DR beta 65-66).
HELIX = {'G-ALPHA1': (50, 86, ()), 'G-ALPHA2': (52, 86, ()),
         'G-ALPHA': (50, 86, ()), 'G-BETA': (52, 86, ('61B', '62'))}

#: Ways to find the groove helices; see the module docstring.
HELIX_METHODS = ('imgt', 'tpm-bench')

#: A domain counts as found when the local alignment to its best template scores at
#: least this fraction of the template's score against itself. G domains score
#: 0.73-1.00, other chains and the other G domains at most 0.14.
_MIN_SCORE = 0.4

#: Templates scoring within this of the best are candidates (see :func:`_best_template`).
_TEMPLATE_MARGIN = 0.25

_LABEL = re.compile(r'^(-?\d+)(?:\.(\d+)|([A-Z]))?$')


# ==========================================================================================
# IMGT data
# ==========================================================================================

@lru_cache(maxsize=1)
def imgt_gdomain_templates() -> pd.DataFrame:
    """IMGT's numbered G-domain sequences: one row per (domain, entry).

    Columns ``domain``, ``species``, ``gene``, ``allele``, ``sequence`` and ``imgt``,
    a list with one IMGT label per residue. Rebuilt from IMGT by
    ``python -m tcrkit._build.fetch_imgt_gdomains``; shipped with the package.
    """
    df = pd.read_csv(DATA / 'imgt_gdomains.csv', keep_default_na=False)
    df['imgt'] = df['imgt'].str.split()
    return df


def _slot_key(label: str) -> tuple[int, int]:
    number, decimal, letter = _LABEL.match(label).groups()
    return int(number), int(decimal) if decimal else (ord(letter) - 64 if letter else 0)


@lru_cache(maxsize=None)
def _slots(domain: str) -> tuple[str, ...]:
    """Every position any template of `domain` uses, in IMGT order."""
    labels = {lab for labs in imgt_gdomain_templates().query('domain == @domain')['imgt']
              for lab in labs}
    return tuple(sorted(labels, key=_slot_key))


@lru_cache(maxsize=1)
def _aligner():
    from Bio import Align
    aligner = Align.PairwiseAligner(mode='local', open_gap_score=-10, extend_gap_score=-1)
    aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
    return aligner


def gdomain_region(label: str) -> str:
    """The IMGT G-domain region of a position: ``'A'``, ``'AB'``, ... ``'helix'``, or ``''``.

    >>> gdomain_region('61A')
    'helix'
    """
    m = _LABEL.match(str(label))
    if m is None or m.group(2):
        return ''
    n = int(m.group(1))
    return next((name for name, (lo, hi) in GDOMAIN_REGIONS.items() if lo <= n <= hi), '')


def is_helix(domain: str, label: str) -> bool:
    """Whether IMGT position `label` of `domain` is a groove-helix residue (:data:`HELIX`).

    >>> is_helix('G-BETA', '61A'), is_helix('G-BETA', '62')
    (True, False)
    """
    if domain not in HELIX:
        return False
    lo, hi, skipped = HELIX[domain]
    m = _LABEL.match(str(label))
    return bool(m and not m.group(2) and lo <= int(m.group(1)) <= hi and label not in skipped)


# ==========================================================================================
# Numbering one sequence
# ==========================================================================================

def _number_one(sequence: str, template) -> dict[int, tuple[str, str]]:
    """``{0-based position: (IMGT label, template residue)}`` for one G domain.

    Aligned residues take the template's label. Residues the template has no partner
    for take, in order, the free IMGT positions between their aligned neighbours
    (template residue ``''``); any left over stay unnumbered. Past the first and last
    aligned residues the numbering continues one to one up to the template's ends, so
    a divergent stretch at a domain's edge (where the local alignment stops) is
    still numbered.
    """
    aln = _aligner().align(sequence, template.sequence)[0]
    out: dict[int, tuple[str, str]] = {}
    slots = _slots(template.domain)
    blocks = list(zip(*aln.aligned))
    for k, ((qs, qe), (ts, te)) in enumerate(blocks):
        for i in range(qe - qs):
            out[qs + i] = (template.imgt[ts + i], template.sequence[ts + i])
        if k + 1 < len(blocks):
            (nqs, _), (nts, _) = blocks[k + 1]
            prev, nxt = template.imgt[te - 1], template.imgt[nts]
            free = slots[slots.index(prev) + 1:slots.index(nxt)]
            for i, pos in enumerate(range(qe, nqs)):
                if i < len(free):
                    out[pos] = (free[i], '')
    (q0, _), (t0, _) = blocks[0]
    for i in range(1, min(q0, t0) + 1):
        out[q0 - i] = (template.imgt[t0 - i], template.sequence[t0 - i])
    (_, q1), (_, t1) = blocks[-1]
    for i in range(min(len(sequence) - q1, len(template.sequence) - t1)):
        out[q1 + i] = (template.imgt[t1 + i], template.sequence[t1 + i])
    return out


def _best_template(sequence: str, domain: str):
    """``(template, score, numbering)`` of the IMGT entry `sequence` is numbered on.

    The score is the best template's. Among the templates scoring close to it, the one
    leaving the fewest residues inside the domain unnumbered is used (an allele's
    insertion may have a free IMGT position in one template and not in another), then
    the best-scoring one.
    """
    aligner = _aligner()
    scored = [(aligner.score(sequence, t.sequence) / aligner.score(t.sequence, t.sequence), t)
              for t in imgt_gdomain_templates().query('domain == @domain').itertuples()]
    best = max(score for score, _ in scored)
    choices = []
    for score, t in scored:
        if score < max(_MIN_SCORE, best - _TEMPLATE_MARGIN):
            continue
        numbering = _number_one(sequence, t)
        span = max(numbering) - min(numbering) + 1
        choices.append((span - len(numbering), -score, t, numbering))
    if not choices:
        return None, best, {}
    _, _, template, numbering = min(choices, key=lambda c: c[:2])
    return template, best, numbering


def _split(sequence: str, first: dict, second: dict, shared: list[int]) -> int:
    """Where to cut `shared` positions between the domain before them and the one after.

    The cut is where the residues' alignment to the two templates scores best:
    ``shared[:k]`` stay with `first`, ``shared[k:]`` with `second`.
    """
    from Bio.Align import substitution_matrices
    blosum = substitution_matrices.load('BLOSUM62')

    def score(numbering, pos):
        partner = numbering[pos][1]
        return blosum[sequence[pos], partner] if partner else 0.0

    return max(range(len(shared) + 1),
               key=lambda k: sum(score(first, p) for p in shared[:k])
               + sum(score(second, p) for p in shared[k:]))


def number_gdomains(sequence: str, domains=None) -> pd.DataFrame:
    """IMGT numbering of the G domains in one MHC chain's sequence.

    Args:
        sequence: one-letter amino-acid sequence.
        domains: G domains to look for (default: all of :data:`GDOMAINS`). Each is
            kept if its best template scores at least 0.4 of a perfect match.
            Residues two domains both claim (where class I alpha1 meets alpha2, often
            when residues are missing there) go to the one they align better to; a
            domain that loses more than half its residues that way is dropped.

    Returns:
        One row per residue: ``aa``, ``domain``, ``imgt``, ``region`` (see
        :func:`gdomain_region`), ``helix`` (:func:`is_helix`) and ``template`` (the IMGT entry it
        was numbered on); ``'-'`` / ``''`` / False outside the G domains.
        ``attrs['scores']`` holds each domain's score.
    """
    domains = GDOMAINS if domains is None else tuple(domains)
    unknown = set(domains) - set(GDOMAINS)
    if unknown:
        raise ValueError(f'unknown G domains {sorted(unknown)}; use some of {list(GDOMAINS)}')
    out = _number_gdomains(sequence, domains)
    return out.copy(deep=True).__finalize__(out)


@lru_cache(maxsize=4096)
def _number_gdomains(sequence: str, domains: tuple[str, ...]) -> pd.DataFrame:
    hits = sorted(((score, d, t, numbering) for d in domains
                   for t, score, numbering in [_best_template(sequence, d)]),
                  key=lambda h: -h[0])
    placed: list[tuple[str, object, dict]] = []
    for score, domain, template, numbering in hits:
        if template is None:
            continue
        mine, theirs = set(numbering), []           # positions kept / given up by others
        for _, _, other in placed:
            shared = sorted(set(numbering) & set(other))
            if not shared:
                continue
            if min(numbering) < min(other):
                k = _split(sequence, numbering, other, shared)
                lost, kept = shared[k:], shared[:k]
            else:
                k = _split(sequence, other, numbering, shared)
                lost, kept = shared[:k], shared[k:]
            mine -= set(lost)
            theirs.append((other, kept))
        if len(mine) * 2 < len(numbering):
            continue
        for other, kept in theirs:
            for p in kept:
                other.pop(p)
        placed.append((domain, template, {p: numbering[p] for p in sorted(mine)}))

    out = pd.DataFrame({'aa': list(sequence), 'domain': '-', 'imgt': '-', 'region': '',
                        'helix': False, 'template': ''})
    for domain, template, numbering in placed:
        pos = np.fromiter(numbering, int)
        out.loc[pos, 'domain'] = domain
        out.loc[pos, 'imgt'] = [lab for lab, _ in numbering.values()]
        out.loc[pos, 'template'] = f'{template.species} {template.gene} {template.allele}'
    out['region'] = [gdomain_region(lab) if lab != '-' else '' for lab in out['imgt']]
    out['helix'] = [is_helix(d, lab) for d, lab in zip(out['domain'], out['imgt'])]
    out.attrs['scores'] = {h[1]: round(h[0], 3) for h in hits}
    return out


# ==========================================================================================
# The tpm-bench helices
# ==========================================================================================

@lru_cache(maxsize=1)
def helix_genes() -> pd.DataFrame:
    """One sequence per MHC gene with its groove-helix ranges (``helix_ind``, 1-based,
    inclusive), from the MHC Motif Atlas (http://mhcmotifatlas.org)."""
    from .pmhc import load_mhc_table

    df = load_mhc_table()
    df = df[df['helix_ind'].notna()]
    df = pd.DataFrame({'allele': df['allele'].values,
                       'sequence': df['helix_ref_sequence'].values,
                       'gene': df['gene'].values,
                       'helix_ind': df['helix_ind'].values})
    df['helix_ind'] = df['helix_ind'].map(ast.literal_eval)
    return df.set_index('gene', drop=False)


def _simple_aligner(mode: str, open_gap: float, extend_gap: float):
    from Bio import Align
    return Align.PairwiseAligner(mode=mode, match_score=1, mismatch_score=-1,
                                 open_gap_score=open_gap, extend_gap_score=extend_gap)


def _best_helix_gene(sequence: str) -> str:
    """The Motif Atlas gene whose sequence best matches, by local alignment score per residue."""
    aligner = _simple_aligner('local', -2, -1)
    best_gene, best = None, 0.0
    for gene, seq in zip(helix_genes()['gene'], helix_genes()['sequence']):
        score = aligner.score(sequence, seq) / max(len(sequence), len(seq))
        if score > best:
            best_gene, best = gene, score
    if best_gene is None:
        raise ValueError('no MHC Motif Atlas gene matches the sequence')
    return best_gene


def _map_ranges(ref_seq: str, target_seq: str, ranges) -> list[list[int | None]]:
    """Carry 1-based ranges on `ref_seq` over to `target_seq` by global alignment.

    A reference position facing a gap maps to the target residue before the gap.
    """
    aln = _simple_aligner('global', -3, -0.25).align(ref_seq, target_seq)[0]
    mapping, r, q = {}, 0, 0
    for a, b in zip(aln[0], aln[1]):
        r += a != '-'
        q += b != '-'
        if a != '-':
            mapping[r] = q
    return [[mapping.get(s), mapping.get(e)] for s, e in ranges]


def tpm_bench_helix(sequence: str) -> np.ndarray:
    """The groove-helix residues of one MHC chain by the ``'tpm-bench'`` method.

    The chain is matched to its closest Motif Atlas gene (:func:`helix_genes`), and the
    gene's helix ranges are carried over by global alignment.

    Returns:
        One True / False per residue.

    Raises:
        ValueError: for sequences under 30 residues, too short to match a gene.
    """
    if len(sequence) < 30:
        raise ValueError(f'sequence too short ({len(sequence)} residues) to match an MHC gene')
    gene = helix_genes().loc[_best_helix_gene(sequence)]
    out = np.zeros(len(sequence), bool)
    for start, end in _map_ranges(gene['sequence'], sequence, gene['helix_ind']):
        # a helix that begins before `sequence` does maps to start 0: keep the part the
        # sequence has. (tpm-bench slices from start - 1 = -1, which drops the helix.)
        out[max(start, 1) - 1:end] = True
    return out
