"""V / J / C marking of a TCR chain, behind ``anarci(per_residue=True, domains=True)``.

ANARCI numbers the V domain only. This adds, per residue:

- ``domain``: ``'V'`` (what ANARCI numbered), ``'C'`` (the constant domain) or ``'-'``.
- ``gene_segment``: which germline gene a residue comes from - ``'L'`` (the V gene's
  leader, where it matches Stitchr's L-PART1 + L-PART2), ``'V'``, ``'J'``, ``'N/D'``
  (CDR3 residues matching neither: non-templated additions, or a D gene), ``'C'``, or
  ``'-'`` (anything else before the V domain, the stretch between the domains, membrane
  and tail).
- The constant domain numbered on IMGT's C-DOMAIN scheme (:mod:`tcrkit.cdomain`), with
  :data:`C_OFFSET` added, so its labels never repeat the V domain's: C position 23 is
  ``'1023'``, ``45.1`` is ``'1045.1'``. ``region`` gets the C strand / loop name
  (``'A'``, ``'AB'``, ... ``'G'``).

Where V ends and J starts inside the CDR3 is read from the germline genes ANARCI called.
Everything up to IMGT 104 is V and everything from 118 on is J. Inside the CDR3:

- V: a residue is V while it and all before it equal ANARCI's V germline at the same
  IMGT position (ANARCI's germlines are laid out on IMGT positions).
- J: counting back from IMGT 118 (the F/W of the J's FGxG motif), a residue is J while
  it and all after it equal the germline J. The full J comes from Stitchr's IMGT data
  (its J-REGION sequences and ``J-region-motifs.tsv``, which says where 118 is), when
  that is installed (``stitchrdl -s human``). Otherwise ANARCI's J germlines are used,
  which only reach back to IMGT 115, so J then starts at 115 at the earliest.

That is the longest germline match, so a non-templated residue that happens to equal
the germline counts as germline: the boundary can sit a residue or two off. Insertions
(``111A``) are never V germline.
"""

from __future__ import annotations

import re
from functools import lru_cache

from .cdomain import C_OFFSET, cdomain_region, number_cdomain, offset_label

__all__ = ['C_OFFSET', 'C_GENES', 'mark_domains']

#: The IMGT constant-domain templates each ANARCI chain type is numbered on. IG chains
#: have none bundled, so their constant domain stays unnumbered.
C_GENES = {'A': ('TRAC',), 'B': ('TRBC2',), 'G': ('TRGC1',), 'D': ('TRDC',)}

#: A remainder shorter than this is not searched for a constant domain.
_MIN_C_LENGTH = 30

_INT = re.compile(r'^\d+$')


@lru_cache(maxsize=1)
def _germlines():
    from anarci.germlines import all_germlines
    return all_germlines


def _germline(kind: str, chain: str, species: str, gene: str) -> str | None:
    """The germline as 128 characters, one per IMGT position ('-' = none), or None."""
    try:
        return _germlines()[kind][chain][species][gene]
    except (KeyError, TypeError):
        return None


@lru_cache(maxsize=None)
def _stitchr_entries(species: str) -> tuple[dict, str | None]:
    """``({(region, allele): nucleotides from its reading frame}, data folder)`` for the
    TR genes in Stitchr's IMGT data; ``({}, None)`` when Stitchr or the data is missing."""
    import os
    try:
        from Stitchr.stitchrfunctions import data_dir
    except Exception:
        return {}, None
    folder = os.path.join(data_dir, str(species).upper())
    if not os.path.isdir(folder):
        return {}, None
    out = {}
    for name in sorted(os.listdir(folder)):
        if not re.fullmatch(r'TR[ABGD]\.fasta', name):
            continue
        with open(os.path.join(folder, name)) as f:
            header, seq = None, []
            for line in [*f, '>']:
                if line.startswith('>'):
                    if header is not None:
                        bits = header.split('|')
                        if len(bits) > 7 and bits[7].strip().isdigit():
                            out[(bits[4], bits[1])] = ''.join(seq).upper()[int(bits[7]) - 1:]
                    header, seq = line[1:].strip(), []
                else:
                    seq.append(line.strip())
    return out, folder


def _translate(nt: str) -> str:
    from Bio.Seq import Seq
    return str(Seq(nt[:len(nt) // 3 * 3]).translate())


@lru_cache(maxsize=None)
def _stitchr_j(species: str) -> dict:
    """``{J allele: (amino acids, index of the IMGT-118 residue)}`` from Stitchr's data.

    Empty when Stitchr or its data for `species` is not installed.
    """
    import os
    entries, folder = _stitchr_entries(species)
    motifs_path = os.path.join(folder, 'J-region-motifs.tsv') if folder else ''
    if not os.path.isfile(motifs_path):
        return {}
    motifs = {}
    with open(motifs_path) as f:
        next(f)
        for line in f:
            bits = line.rstrip('\n').split('\t')
            if len(bits) >= 5 and bits[4].lstrip('-').isdigit():
                motifs[bits[0]] = int(bits[4])
    out = {}
    for (region, gene), nt in entries.items():
        if region == 'J-REGION' and gene in motifs:
            aa = _translate(nt)
            k = len(aa) + motifs[gene]
            if 0 <= k < len(aa):
                out[gene] = (aa, k)
    return out


@lru_cache(maxsize=None)
def _stitchr_leaders(species: str) -> dict:
    """``{V allele: leader amino acids}`` (L-PART1 + L-PART2, with the codon it shares
    with the V-REGION) from Stitchr's data; empty when it is not installed."""
    entries, _ = _stitchr_entries(species)
    out = {}
    for (region, gene), nt in entries.items():
        if region == 'L-PART1+L-PART2':
            v = entries.get(('V-REGION', gene), '')
            n = -(-len(nt) // 3)                         # codons, the split one included
            aa = _translate(nt + v[:3])[:n]
            if aa:
                out[gene] = aa
    return out


def _allele(table: dict, gene):
    """`gene`'s entry in `table`, else another allele of the same gene, else None."""
    if gene in table:
        return table[gene]
    base = str(gene).split('*')[0]
    same = [v for k, v in table.items() if k.split('*')[0] == base]
    return same[0] if same else None


def _leader_match(seq: str, first: int, leader: str) -> tuple[int, int]:
    """``(matched residues, gap)``: the longest run of `seq` before `first` (the first
    V-domain residue) equal to the end of `leader`, counting back. The leader may end
    up to 3 residues before `first` (V residues ANARCI leaves unnumbered: the `gap`), or
    one residue into it (the codon the leader shares with the V-REGION)."""
    best = (0, 0)
    for gap in range(0, 4):
        for lead in (leader, leader[:-1]):
            end = first - gap
            n = 0
            while n < min(end, len(lead)) and seq[end - 1 - n] == lead[-1 - n]:
                n += 1
            if n > best[0]:
                best = (n, gap)
    return best


def _mark_leader(seq: str, first: int, leaders, segment: list) -> None:
    """Mark ``'L'`` on the residues before `first` that match a V-gene leader in
    `leaders` (the best of them), and ``'V'`` on any V residues between it and `first`."""
    if first == 0:
        return
    n, gap = max((_leader_match(seq, first, x) for x in leaders if x), default=(0, 0))
    if n < 3:                       # too short to tell from chance
        return
    end = first - gap
    for i in range(end - n, end):
        segment[i] = 'L'
    for i in range(end, first):
        segment[i] = 'V'


def _j_germline(species, gene):
    """``(amino acids, index of IMGT 118)`` of a J allele from Stitchr, or None."""
    return _allele(_stitchr_j(species), gene)


def _matches(germ: str | None, label: str, aa: str) -> bool:
    if germ is None or not _INT.match(label):
        return False
    p = int(label)
    return 1 <= p <= len(germ) and germ[p - 1] == aa


def mark_domains(seq: str, imgt: list[str], region: list[str], row) -> tuple[list, ...]:
    """``(domain, gene_segment, imgt, region)`` lists for one sequence.

    `imgt` / `region` are the V-domain labels ANARCI gave each residue (``'-'`` / ``''``
    elsewhere); `row` is that sequence's ANARCI summary (``chain_type``, ``species``,
    ``v_gene``, ``j_gene``), or None if it did not number.
    """
    n = len(seq)
    imgt, region = list(imgt), list(region)
    domain, segment = ['-'] * n, ['-'] * n
    inside = [i for i, lab in enumerate(imgt) if lab != '-']
    if not inside or row is None:
        return domain, segment, imgt, region
    chain, species = row.get('chain_type'), row.get('species')

    # the V domain, and which germline gene each of its residues comes from
    v_germ = _germline('V', chain, species, row.get('v_gene'))
    j_germ = _germline('J', chain, species, row.get('j_gene'))
    j_full = _j_germline(species, row.get('j_gene'))
    cdr3, at118 = [], None
    for i in inside:
        domain[i] = 'V'
        num = int(re.match(r'\d+', imgt[i]).group())
        if num <= 104:
            segment[i] = 'V'
        elif num >= 118:
            segment[i] = 'J'
            if imgt[i] == '118':
                at118 = i
        else:
            cdr3.append(i)
            segment[i] = 'N/D'
    # the V gene's leader before it. ANARCI's species call can be off (a human chain called
    # mouse), and the match is checked residue by residue, so both are tried.
    leaders = [_allele(_stitchr_leaders(sp), row.get('v_gene'))
               for sp in dict.fromkeys([species, 'human', 'mouse'])]
    _mark_leader(seq, inside[0], leaders, segment)
    for i in cdr3:                                   # V runs on into the CDR3 ...
        if not _matches(v_germ, imgt[i], seq[i]):
            break
        segment[i] = 'V'
    if j_full is not None and at118 is not None:     # ... and J back from 118
        aa, k118 = j_full
        for back, i in enumerate(reversed([i for i in cdr3 if i < at118]), start=1):
            if segment[i] == 'V' or back > k118 or aa[k118 - back] != seq[i]:
                break
            segment[i] = 'J'
    else:
        for i in reversed(cdr3):
            if segment[i] == 'V' or not _matches(j_germ, imgt[i], seq[i]):
                break
            segment[i] = 'J'

    # the constant domain, in what follows the V domain
    start = inside[-1] + 1
    genes = C_GENES.get(chain)
    if genes and n - start >= _MIN_C_LENGTH:
        c = number_cdomain(seq[start:], genes=genes)
        for k, lab in enumerate(c['imgt']):
            if lab != '-':
                i = start + k
                domain[i], segment[i] = 'C', 'C'
                imgt[i], region[i] = offset_label(lab), cdomain_region(lab)
    return domain, segment, imgt, region
