"""MHC chains: their gene, IMGT numbering and groove helices.

:func:`annotate_mhc` is the MHC counterpart of :func:`tcrkit.anarci`; :func:`best_mhc`,
:func:`mhc_span` and :func:`mhc_chain` are the gene matching behind it, for one sequence.

Numbering follows IMGT's G-DOMAIN scheme, numbered on the templates of
https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html
(Lefranc et al., *Dev. Comp. Immunol.* 29, 917-938, 2005).
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
import pandas as pd

from ._util import ERROR_COL, collapse_columns, collect, segments as _segments, squeeze as _squeeze
#: Membrane regions carried over from the matched allele, and the span column each is
#: recorded in - the same names :func:`tcrkit.mhc_sequence` cuts on.
from .pmhc import TOPOLOGY_LETTERS
from .pmhc import TOPOLOGY as _TOPOLOGY_COLS
from .gdomain import _MIN_SCORE, GDOMAINS, HELIX_METHODS, number_gdomains, tpm_bench_helix

__all__ = ['annotate_mhc', 'gene_domains', 'best_mhc', 'mhc_chain', 'mhc_span',
           'mhc_species', 'mhc_database', 'B2M_SEQUENCE', 'HELIX_METHODS', 'COLLAPSIBLE']

_COLUMNS = ['gene', 'allele', 'allele_identity', 'allele_ties', 'allele_source', 'mhc_class', 'score', 'span', 'domains', 'domain_spans',
            'helix', 'regions', *_TOPOLOGY_COLS, ERROR_COL]
_RESIDUE_COLUMNS = ['position', 'aa', 'gene', 'allele', 'domain', 'imgt', 'region',
                    'helix', 'topology']

_TOPOLOGY = tuple(_TOPOLOGY_COLS.items())


# ----------------------------------------------- matching a sequence to its gene

#: Human B2M, the database's B2M entry.
B2M_SEQUENCE = ('MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFT'
                'PTEKDEYACRVNHVTLSQPKIVKWDRDM')

#: An MHC gene is called for sequences of at least this many residues whose best match
#: scores more than :data:`_MHC_MIN_SCORE` per residue.
_MHC_MIN_LENGTH = 30
_MHC_MIN_SCORE = 0.07


def _local_aligner():
    from Bio import Align
    return Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=-1,
                                 open_gap_score=-2, extend_gap_score=-1)


@lru_cache(maxsize=1)
def mhc_database() -> pd.DataFrame:
    """The MHC sequences genes are matched on: HLA class I, HLA-E, class II, mouse and
    B2M, then the other species' chains (TCR3d's: macaque, pig, chicken, ...).

    Columns ``allele``, ``sequence``, ``gene`` and ``class`` (1 or 2; 1 for B2M). The
    order matters: of equally good matches the first one wins.
    """
    from .pmhc import load_mhc_table, region_of

    db = load_mhc_table()
    db = db[db['in_detect']].sort_values('detect_rank')
    # match on the mature chain: a deposited structure starts where the signal peptide
    # was cleaved, so the leader would only dilute the score
    return pd.DataFrame({
        'allele': db['allele'].values,
        'sequence': [region_of(r, 'mature')[0] for _, r in db.iterrows()],
        'gene': db['gene'].values,
        'class': db['detect_class'].astype(int).values,
    })


def mhc_chain(gene: str | None) -> str | None:
    """Which chain of an MHC molecule `gene` codes for: ``'heavy'`` (class I heavy
    chain), ``'alpha'`` / ``'beta'`` (class II), ``'B2M'``, or None for a gene that
    is not in :func:`mhc_database`.

    Class comes from the database; for class II the letter before the gene number
    says which chain - DRA, DQA1, H2 AA and EA are alpha, DRB1, DQB1, AB and EB beta.
    """
    if gene == 'B2M':
        return 'B2M'
    classes = _gene_classes()
    if gene not in classes:
        return None
    if classes[gene] == 1:
        return 'heavy'
    return 'alpha' if re.sub(r'\d+$', '', gene).endswith('A') else 'beta'


def mhc_species(allele: str | None) -> str | None:
    """The species of an allele of :func:`mhc_database` (``'human'``, ``'mouse'``,
    ``'pig'``, ``'chicken'``, ...), or None for an allele that is not in it.

    >>> mhc_species('A0201'), mhc_species('not-an-allele')
    ('human', None)

    The species of a chain's closest allele: ``mhc_species(best_mhc(sequence)[1])``.
    An allele of the MHC Motif Atlas (what :func:`annotate_mhc` reports by default) is
    looked up there too, once its table is cached.
    """
    found = _allele_species().get(allele)
    if found is None and isinstance(allele, str):
        from .motifatlas import species_of
        found = species_of(allele)
    return found


@lru_cache(maxsize=1)
def _allele_species() -> dict[str, str]:
    from .pmhc import load_mhc_table
    table = load_mhc_table()
    return dict(zip(table['allele'], table['species']))


@lru_cache(maxsize=1)
def _gene_classes() -> dict[str, int]:
    db = mhc_database()
    return dict(zip(db['gene'], db['class']))


def best_mhc(sequence: str, db: pd.DataFrame | None = None
             ) -> tuple[str | None, str | None, float]:
    """Best-matching MHC gene, allele, and score per residue (``None, None, 0`` if none).

    The score is the local alignment score divided by the longer sequence's length; of
    equally good matches the first in :func:`mhc_database` wins. Results against the
    bundled database (`db` None) are cached by sequence.
    """
    if db is None:
        return _best_bundled(sequence)
    return _best_match(sequence, db)


@lru_cache(maxsize=4096)
def _best_bundled(sequence: str) -> tuple[str | None, str | None, float]:
    return _best_match(sequence, mhc_database())


#: How many of the closest genes go through to the second pass. Two is not enough; the
#: answer stops changing at three.
_GENES_CONSIDERED = 4


def _scan(sequence: str, db: pd.DataFrame, rows) -> tuple[int | None, float]:
    """The row of `db` (by position) that `sequence` matches best, and its score."""
    aligner = _local_aligner()
    seqs = db['sequence']
    best_i, best = None, 0.0
    for i in rows:
        score = aligner.score(sequence, seqs.iat[i]) / max(len(sequence), len(seqs.iat[i]))
        if score > best:
            best_i, best = i, score
    return best_i, best


def _best_match(sequence: str, db: pd.DataFrame) -> tuple[str | None, str | None, float]:
    """Gene, allele and score per residue of the closest database entry.

    Two passes: one representative allele per gene narrows the field to the
    :data:`_GENES_CONSIDERED` closest genes, then only those genes' alleles are scored.
    A few dozen alignments instead of one per allele, for the same answer - carrying
    fewer than three genes forward does change it, because a gene like DRB1 spans
    serotypes too different for one representative to speak for.
    """
    first: dict = {}                             # gene -> first row, the representative
    for i, gene in enumerate(db['gene']):
        first.setdefault(gene, i)

    aligner = _local_aligner()
    seqs = db['sequence']
    ranked = sorted(
        ((aligner.score(sequence, seqs.iat[i]) / max(len(sequence), len(seqs.iat[i])), g)
         for g, i in first.items()), reverse=True)
    if not ranked or ranked[0][0] <= 0:
        return None, None, 0.0

    keep = {g for _, g in ranked[:_GENES_CONSIDERED]}
    best_i, best = _scan(sequence, db, [i for i, g in enumerate(db['gene']) if g in keep])
    if best_i is None:
        return None, None, 0.0
    return db['gene'].iat[best_i], db['allele'].iat[best_i], best


def mhc_span(sequence: str, db: pd.DataFrame | None = None) -> tuple[int, int] | None:
    """Where the MHC part lies in `sequence`: ``(start, end)``, 0-based, end exclusive.

    The part covered by the local alignment to its best-matching MHC sequence; residues
    outside it (tags, linkers, a fused peptide) are not part of the MHC. None if no MHC
    sequence matches.
    """
    gene, allele, _ = best_mhc(sequence, db)
    if gene is None:
        return None
    db = mhc_database() if db is None else db
    ref = db.loc[db['allele'] == allele, 'sequence'].iloc[0]
    coords = _local_aligner().align(sequence, ref)[0].coordinates
    return int(coords[0][0]), int(coords[0][-1])


# ------------------------------------------------------------------- genes

def gene_domains(gene: str | None) -> tuple[str, ...]:
    """The G domains a chain of `gene` has: ``('G-ALPHA1', 'G-ALPHA2')`` for a class I
    heavy chain, ``('G-ALPHA',)`` / ``('G-BETA',)`` for class II alpha / beta, ``()``
    for B2M, and all four (:data:`~tcrkit.gdomain.GDOMAINS`) for None or an unknown gene.
    """
    kind = mhc_chain(gene)
    if kind is None:
        return GDOMAINS
    return {'heavy': ('G-ALPHA1', 'G-ALPHA2'), 'alpha': ('G-ALPHA',),
            'beta': ('G-BETA',), 'B2M': ()}[kind]


def _ranges(positions) -> str:
    """1-based positions -> ``'50-86, 140-176'``."""
    positions = sorted(positions)
    if not positions:
        return ''
    runs, start, prev = [], positions[0], positions[0]
    for p in positions[1:]:
        if p != prev + 1:
            runs.append((start, prev))
            start = p
        prev = p
    runs.append((start, prev))
    return ', '.join(f'{a}-{b}' if a != b else f'{a}' for a, b in runs)


def _topology_of(seq: str, allele: str | None):
    """Carry the matched allele's membrane topology onto `seq`.

    The reference spans are DeepTMHMM's, measured on the allele's own sequence and
    carried over here by alignment - so a chain that is not exactly a bundled allele,
    or is only a fragment of one, still gets its regions marked.

    Returns:
        ``(prediction, per-residue labels, {region: '1-283' ...})``.
    """
    from .gdomain import _map_ranges
    from .pmhc import _spans, load_mhc_table

    labels = np.full(len(seq), '-', dtype=object)
    if allele is None:
        return None, labels, {}
    table = load_mhc_table().set_index('allele')
    if allele not in table.index:
        return None, labels, {}
    row = table.loc[allele]

    flat, owner = [], []
    for name, col in _TOPOLOGY:
        for span in _spans(row.get(col)):
            flat.append(span)
            owner.append(name)
    if not flat:
        return row.get('prediction'), labels, {}

    spans: dict[str, list[int]] = {}
    mapped = _map_ranges(row['sequence'], seq, flat)
    for name, (ref_start, ref_end), (start, end) in zip(owner, flat, mapped):
        if start is None or end is None:
            continue
        start = max(1, start)
        # a region the sequence does not reach collapses onto its last residue; that is
        # an alignment artefact, not a one-residue transmembrane helix
        if ref_end - ref_start >= 4 and end - start <= 1:
            continue
        labels[start - 1:end] = name
        spans.setdefault(name, []).extend(range(start, end + 1))
    return row.get('prediction'), labels, {k: _ranges(v) for k, v in spans.items()}


def _match_allele(seq: str, ref: str) -> tuple[int, int, float]:
    """Where `seq` aligns to the allele sequence `ref`, and how identical it is there.

    Returns ``(start, end, identity)``: the 0-based, end-exclusive span of `seq` covered
    by the local alignment (as :func:`mhc_span`), and the fraction of the
    alignment's columns, gaps included, where the two residues are the same.
    """
    aln = _local_aligner().align(seq, ref)[0]
    same = aligned = 0
    for (s0, s1), (r0, _) in zip(*aln.aligned):
        aligned += s1 - s0
        same += sum(a == b for a, b in zip(seq[s0:s1], ref[r0:r0 + s1 - s0]))
    coords = aln.coordinates
    start, end = int(coords[0][0]), int(coords[0][-1])
    columns = (end - start) + int(coords[1][-1] - coords[1][0]) - aligned
    return start, end, same / columns if columns else 0.0


def _blank(seq: str) -> pd.DataFrame:
    """A per-residue table with nothing numbered."""
    return pd.DataFrame({'aa': list(seq), 'domain': '-', 'imgt': '-', 'region': '',
                         'helix': False, 'topology': '-', 'template': ''})


# ------------------------------------------------------------- one sequence

#: The membrane regions a residue before the first numbered one may be in for
#: :func:`_number_the_leader` to number it.
_LEADER_TOPOLOGY = ('signal_peptide', 'extracellular')


def _number_the_leader(num: pd.DataFrame, labels) -> pd.DataFrame:
    """Number the residues before the first IMGT position backwards from it.

    The usual convention for leader peptides: counting back from IMGT 1, the residue
    before it is ``-1`` (there is no 0); counting back from a negative IMGT position
    (class II G-ALPHA starts at ``-2``) it continues ``-3``, ``-4`` ... So the whole
    N-terminus, the signal peptide and any mature residues IMGT leaves out (DRB1 has
    four before IMGT 1), is one run of negative numbers on the same axis;
    ``topology`` still says which of them are the signal peptide.

    Only applied when the first numbered residue is IMGT ``1`` or a negative position,
    and everything before it is signal peptide or extracellular. Nothing after a
    domain (the membrane anchor, the tail) is touched.
    """
    numbered = [i for i, v in enumerate(num['imgt']) if v != '-']
    if not numbered or numbered[0] == 0:
        return num
    first = numbered[0]
    label = num['imgt'].iloc[first]
    if not re.fullmatch(r'-?\d+', label) or int(label) > 1:
        return num
    if any(labels[i] not in _LEADER_TOPOLOGY for i in range(first)):
        return num
    top = int(label) if int(label) < 0 else 0            # the count skips 0
    imgt = list(num['imgt'])
    for i in range(first):
        imgt[i] = str(top - (first - i))
    return num.assign(imgt=imgt)


def _number_clike(seq: str, num: pd.DataFrame, kind: str | None) -> pd.DataFrame:
    """`num` with the C-LIKE domain after its last numbered residue added: class I
    alpha3, class II alpha2 / beta2, or B2M's domain, on IMGT's C-LIKE scheme plus
    :data:`~tcrkit.cdomain.C_OFFSET` (``'1023'``)."""
    from .cdomain import MHC_CLIKE_GENES, number_cdomain, offset_label

    if kind is None:
        genes = [g for k, gs in MHC_CLIKE_GENES.items() if k != 'B2M' for g in gs]
    else:
        genes = MHC_CLIKE_GENES[kind]
    numbered = np.flatnonzero((num['imgt'] != '-').to_numpy())
    start = int(numbered[-1]) + 1 if len(numbered) and kind != 'B2M' else 0
    if len(seq) - start < 40:
        return num
    c = number_cdomain(seq[start:], genes=genes)
    inside = np.flatnonzero((c['imgt'] != '-').to_numpy())
    rows = inside + start
    num = num.copy()
    num.loc[rows, 'domain'] = 'C-LIKE'
    num.loc[rows, 'imgt'] = [offset_label(x) for x in c['imgt'].iloc[inside]]
    num.loc[rows, 'region'] = c['region'].iloc[inside].to_numpy()
    num.loc[rows, 'template'] = c['template'].iloc[inside].to_numpy()
    return num


def _number_the_b2m_leader(num: pd.DataFrame, labels) -> pd.DataFrame:
    """B2M's leader counted back from its first C-LIKE position (``-1`` before 1)."""
    numbered = [i for i, v in enumerate(num['imgt']) if v != '-']
    if not numbered:
        return num
    first = numbered[0]
    plain = num.assign(imgt=['-'] * first + ['1'] + ['-'] * (len(num) - first - 1))
    lead = _number_the_leader(plain, labels)['imgt'].iloc[:first].tolist()
    return num.assign(imgt=lead + list(num['imgt'].iloc[first:]))


#: Where :func:`annotate_mhc` takes the ``allele`` from.
ALLELE_SOURCES = ('motifatlas', 'uniprot')


def _annotate_one(seq: str, helix: str, mark_regions: bool = False,
                  c_like: bool = False, allele_source: str = 'motifatlas'
                  ) -> tuple[dict, pd.DataFrame]:
    """``(row of the wide table, per-residue table)`` for one cleaned sequence."""
    if len(seq) < _MHC_MIN_LENGTH:
        return ({ERROR_COL: f'too short for an MHC chain ({len(seq)} residues, needs '
                            f'{_MHC_MIN_LENGTH})'}, _blank(seq))
    rec: dict = {}
    gene, allele, score = best_mhc(seq)
    rec['score'] = round(score, 3)
    if gene is not None and score > _MHC_MIN_SCORE:
        db = mhc_database()
        ref = db.loc[db['allele'] == allele].iloc[0]
        start, end, identity = _match_allele(seq, ref['sequence'])
        cls = ref['class']
        rec.update(gene=gene, allele=allele, allele_identity=round(identity, 3),
                   allele_source='uniprot',
                   mhc_class='-' if gene == 'B2M' else {1: 'I', 2: 'II'}[cls],
                   span=f'{start + 1}-{end}')
        if allele_source == 'motifatlas' and gene != 'B2M':
            from .motifatlas import allele_for
            hit = allele_for(seq, gene, mhc_species(allele))
            if hit is not None:
                rec.update(allele=hit[0], allele_identity=hit[1], allele_ties=hit[2],
                           allele_source='motifatlas')
    else:
        gene, allele = None, None

    # the membrane regions, carried over from the bundled allele that matched
    _, labels, spans = _topology_of(seq, allele)
    rec['regions'] = ', '.join(name for name, _ in _TOPOLOGY if spans.get(name)) or None
    if mark_regions:
        rec['topology_seq'] = ''.join(TOPOLOGY_LETTERS.get(x, '-') for x in labels)
    rec.update({name: spans.get(name) for name, _ in _TOPOLOGY})

    domains = gene_domains(gene)
    if not domains:                                  # B2M: no groove domain
        num = _blank(seq).assign(topology=labels)
        if c_like:                                   # but a C-LIKE one
            num = _number_the_b2m_leader(_number_clike(seq, num, 'B2M'), labels)
            inside = np.flatnonzero((num['domain'] == 'C-LIKE').to_numpy()) + 1
            if len(inside):
                rec.update(domains='C-LIKE', domain_spans=_ranges(inside))
        return {**rec, ERROR_COL: None}, num

    num = number_gdomains(seq, domains)
    found = num[num['domain'] != '-']
    if found.empty:
        best = max(num.attrs['scores'].values())
        what = f'{gene} chain' if gene else 'sequence (no MHC gene matches either)'
        rec[ERROR_COL] = (f"no IMGT {' / '.join(domains)} domain found in this {what} "
                          f'(best score {best:.2f}, needs {_MIN_SCORE})')
        return rec, _number_the_leader(num.assign(topology=labels), labels)
    if helix == 'tpm-bench':
        num['helix'] = tpm_bench_helix(seq)
    num = _number_the_leader(num.assign(topology=labels), labels)
    if c_like:
        num = _number_clike(seq, num, mhc_chain(gene))
        found = num[num['domain'] != '-']
    order = list(dict.fromkeys(found['domain']))
    rec.update(domains=', '.join(order),
               domain_spans=', '.join(_ranges(found.index[found['domain'] == d] + 1)
                                      for d in order),
               helix=_ranges(np.flatnonzero(num['helix'].to_numpy()) + 1),
               **{ERROR_COL: None})
    return rec, num


# -------------------------------------------------------------------- segments

#: The per-residue columns :func:`annotate_mhc` can `collapse` on.
COLLAPSIBLE = ('domain', 'region', 'helix', 'topology')


# -------------------------------------------------------------- the capability

def annotate_mhc(sequences, helix: str = 'imgt', per_residue: bool = False,
                 mark_regions: bool = False, col_mapping=None,
                 squeeze: bool = False, collapse=None, domains: bool = False,
                 show=False, allele_source: str = 'motifatlas') -> pd.DataFrame:
    """Which MHC gene each sequence is, its IMGT numbering and its groove helices.

    Each sequence is matched to the bundled alleles for its gene, then to that gene's
    MHC Motif Atlas alleles for its closest allele, and its G domains are numbered on
    IMGT's G-DOMAIN scheme:
    https://www.imgt.org/IMGTrepertoireMH/Proteins/protein/G-DOMAIN/Gdomains.html

    >>> out = annotate_mhc('GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEG'
    ...                    'PEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGY'
    ...                    'HQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRR'
    ...                    'YLENGKETLQRT')        # HLA-A*02:01, alpha1 + alpha2
    >>> out[['gene', 'allele', 'domains', 'helix']].iloc[0].tolist()
    ['A', 'A0201', 'G-ALPHA1, G-ALPHA2', '50-86, 140-176']

    Args:
        sequences: one amino-acid sequence or any collection (list, dict, DataFrame...).
        helix: ``'imgt'`` (default, fixed IMGT positions) or ``'tpm-bench'`` (MHC Motif
            Atlas ranges, http://mhcmotifatlas.org).
        per_residue: one row per residue instead of one per sequence.
        collapse: one row per run of residues sharing these per-residue column(s), any of
            :data:`COLLAPSIBLE`; ``True`` means ``['domain', 'region']``.
        domains: also number the membrane-proximal C-LIKE domain - class I alpha3,
            class II alpha2 / beta2, and B2M's only domain - on IMGT's C-LIKE scheme
            (:mod:`tcrkit.cdomain`) plus 1000, so its labels never repeat the G
            domain's: C-LIKE position 23 is ``'1023'``, ``45.1`` is ``'1045.1'``.
            Those residues get ``domain`` ``'C-LIKE'`` and their strand / loop in
            ``region`` (``'A'``, ``'AB'``, ... ``'G'``), and ``'C-LIKE'`` is listed in
            ``domains`` / ``domain_spans``. B2M's leader is numbered back from it.
        show: print each sequence with its annotation under it
            (:func:`tcrkit.print_residues`, in colour) and return None, instead of a table.
            ``True``, or a dict of ``print_residues`` options, e.g.
            ``{'width': None}``, ``{'color': False}``.
        mark_regions: add ``topology_seq``, one letter per residue (``S``/``O``/``M``/
            ``I``/``-``), as in :func:`tcrkit.mhc_sequence`.
        col_mapping: where the sequence lives in a frame, e.g. ``{'sequence': 'mhc_seq'}``.
        squeeze: return the plain result instead of a frame.
        allele_source: where ``allele`` comes from. ``'motifatlas'`` (default): the
            closest of the gene's MHC Motif Atlas alleles (13,966 class I, 865 class II;
            :func:`tcrkit.motifatlas.closest_allele`), downloaded on first use. Class I
            is compared over the groove, all the Motif Atlas has. Falls back to the
            bundled reference, with a warning, when the Motif Atlas cannot be had, and
            for genes it lacks (B2M, most non-human). ``'uniprot'``: the closest
            bundled allele, as before 0.1.11. The gene, numbering, helices and
            membrane regions do not depend on it.

    Returns:
        Per sequence: ``gene``, ``allele`` (the closest allele - the gene is reliable,
        the allele approximate), ``allele_identity`` (1.0 = identical where aligned),
        ``allele_ties`` (how many other Motif Atlas alleles are as close; they share
        that stretch of sequence), ``allele_source`` (``'motifatlas'`` or
        ``'uniprot'``),
        ``mhc_class``, ``score``, ``span``, ``domains``, ``domain_spans``, ``helix``
        (1-based ranges), ``regions`` plus one range column per membrane region, and
        ``error``.

        Per residue: ``position``, ``aa``, ``gene``, ``allele``, ``domain``, ``imgt``,
        ``region`` (IMGT strand A-D, turn, or ``'helix'`` = IMGT 50-92), ``helix``
        (bool) and ``topology``. The ``helix`` column is the reference for the groove
        helices: the residues TCR docking angles are measured on, narrower than
        ``region == 'helix'``.

        With `collapse`: the column value(s), ``start`` / ``end``, ``length``,
        ``imgt_start`` / ``imgt_end`` and ``sequence`` per segment.
    """
    if helix not in HELIX_METHODS:
        raise ValueError(f'helix must be one of {HELIX_METHODS}, got {helix!r}')
    if allele_source not in ALLELE_SOURCES:
        raise ValueError(f'allele_source must be one of {ALLELE_SOURCES}, '
                         f'got {allele_source!r}')
    collapse = collapse_columns(collapse, COLLAPSIBLE, default=['domain', 'region'])
    if collapse is not None:
        per_residue = True
    inp = collect(('sequence',), (sequences,), col_mapping=col_mapping)
    recs, parts, templates = [], [], {}
    for idx, seq in zip(inp.index, inp['sequence']):
        seq = '' if seq is None or (isinstance(seq, float) and pd.isna(seq)) else str(seq)
        seq = re.sub(r'\s+', '', seq).upper()
        if not seq or not re.fullmatch(r'[A-Z]+', seq):
            recs.append({ERROR_COL: 'not an amino-acid sequence'})
            continue
        rec, num = _annotate_one(seq, helix, mark_regions, domains, allele_source)
        recs.append(rec)
        if per_residue or show:
            part = num.drop(columns='template')
            part.insert(0, 'position', np.arange(1, len(num) + 1))
            part.insert(2, 'gene', rec.get('gene'))
            part.insert(3, 'allele', rec.get('allele'))
            part.index = [idx] * len(part)
            parts.append(part)
            templates[idx] = '; '.join(t for t in dict.fromkeys(num['template']) if t)
    if show:
        from .display import show as _show
        return _show(pd.concat(parts), show) if parts else None
    if per_residue:
        out = (pd.concat(parts) if parts else pd.DataFrame(columns=_RESIDUE_COLUMNS))
        if collapse is not None:
            out = _segments(out, collapse)
        out.attrs['templates'] = templates
        return _squeeze(out) if squeeze else out
    columns = list(_COLUMNS)
    if mark_regions:
        columns.insert(columns.index('regions') + 1, 'topology_seq')
    out = pd.DataFrame(recs, index=inp.index).reindex(columns=columns)
    out.insert(0, 'sequence', list(inp['sequence']))
    return _squeeze(out, echo=('sequence',)) if squeeze else out
