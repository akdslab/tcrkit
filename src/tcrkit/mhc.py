"""MHC chains: their gene, IMGT numbering and groove helices.

One capability, :func:`annotate_mhc`, the MHC counterpart of :func:`tcrkit.anarci`::

    annotate_mhc(seq)                          # one row per sequence: gene, domains, helix
    annotate_mhc(seq, per_residue=True)        # one row per residue: IMGT label, helix
    annotate_mhc({'A2': seq_a2, 'DRB1': seq_drb1})     # any collection, as everywhere

and the gene matching behind it, for one sequence::

    best_mhc(seq)                    # ('A', 'A0201', 0.53): gene, allele, score per residue
    mhc_span('GGGGSGGGGS' + seq)     # (10, 192): where the MHC part is, without the tag
    mhc_chain('DRB1')                # 'beta': which chain of the molecule a gene codes for

Each sequence is first matched to the bundled HLA and B2M sequences
(:func:`best_mhc`) for its gene and closest allele.

The gene call is reliable; the allele is approximate. The bundled database holds 184
alleles (human class I and II, HLA-E/G, a few mouse H2 and B2M), a small fraction of
the tens of thousands of known HLA alleles, so ``allele`` is the *closest of these*, not
a typing result. ``allele_identity`` says how close: 1.0 means your sequence is identical
to that allele over the part that aligns (other alleles may still be identical there
too, for example ones that differ only outside the groove), and anything lower means
the true allele is not bundled. The gene says which
groove domains (G domains) the chain has: G-ALPHA1 and G-ALPHA2 for a class I heavy
chain, G-ALPHA for a class II alpha chain, G-BETA for a class II beta chain, none for
B2M. Those are numbered on IMGT's G-DOMAIN scheme (:func:`tcrkit.gdomain.number_gdomains`):
strands A-D (positions 1-49) and the helix (50-92B), with insertions such as ``'61A'``.
An IMGT position means the same place in every allele. A sequence that matches no gene
is still numbered, trying all four domains.

The ``helix`` column marks the groove-helix residues that TCR docking angles are measured
on, by one of :data:`HELIX_METHODS`:

- ``'imgt'`` (default): fixed IMGT positions (:data:`tcrkit.gdomain.HELIX`), so the same
  residues in every allele: HLA-A2 50-86 and 140-176, HLA-DR alpha 46-78, beta 54-64
  and 67-91. They are part of IMGT's helix region (50-92B), not all of it.
- ``'tpm-bench'``: the helix ranges of the closest MHC Motif Atlas gene
  (http://mhcmotifatlas.org), carried over by sequence alignment, as the tpm-bench
  benchmark does. The same residues as ``'imgt'`` on class I; on class II the ends
  differ by a few residues, gene by gene.
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
           'mhc_database', 'B2M_SEQUENCE', 'HELIX_METHODS', 'COLLAPSIBLE']

_COLUMNS = ['gene', 'allele', 'allele_identity', 'mhc_class', 'score', 'span', 'domains', 'domain_spans',
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
    """The MHC sequences genes are matched on: HLA class I, HLA-E, class II, then B2M.

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

def _number_the_leader(num: pd.DataFrame, labels) -> pd.DataFrame:
    """Number signal-peptide residues backwards from IMGT 1: ``-1`` is the last one.

    The usual convention for leader peptides, and for these chains the frames
    line up - the G domain starts at the mature N-terminus, so the residue before IMGT
    1 is -1. Only applied when that actually holds: the first numbered residue must be
    IMGT ``1`` and everything before it the signal peptide. Nothing else unnumbered
    (the C-like domain, the membrane anchor, the tail) is touched - those belong to
    other IMGT schemes, and a label here would mean two frames in one column.
    """
    numbered = [i for i, v in enumerate(num['imgt']) if v != '-']
    if not numbered:
        return num
    first = numbered[0]
    if num['imgt'].iloc[first] != '1' or first == 0:
        return num
    if any(labels[i] != 'signal_peptide' for i in range(first)):
        return num
    imgt = list(num['imgt'])
    for i in range(first):
        imgt[i] = str(i - first)                 # -first ... -1
    return num.assign(imgt=imgt)


def _annotate_one(seq: str, helix: str, mark_regions: bool = False
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
                   mhc_class='-' if gene == 'B2M' else {1: 'I', 2: 'II'}[cls],
                   span=f'{start + 1}-{end}')
    else:
        gene = None

    # the membrane regions, carried over from whichever allele matched
    _, labels, spans = _topology_of(seq, rec.get('allele'))
    rec['regions'] = ', '.join(name for name, _ in _TOPOLOGY if spans.get(name)) or None
    if mark_regions:
        rec['topology_seq'] = ''.join(TOPOLOGY_LETTERS.get(x, '-') for x in labels)
    rec.update({name: spans.get(name) for name, _ in _TOPOLOGY})

    domains = gene_domains(gene)
    if not domains:                                  # B2M: no groove domain to number
        return {**rec, ERROR_COL: None}, _blank(seq).assign(topology=labels)

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
    order = list(dict.fromkeys(found['domain']))
    rec.update(domains=', '.join(order),
               domain_spans=', '.join(_ranges(found.index[found['domain'] == d] + 1)
                                      for d in order),
               helix=_ranges(np.flatnonzero(num['helix'].to_numpy()) + 1),
               **{ERROR_COL: None})
    return rec, _number_the_leader(num.assign(topology=labels), labels)


# -------------------------------------------------------------------- segments

#: The per-residue columns :func:`annotate_mhc` can `collapse` on.
COLLAPSIBLE = ('domain', 'region', 'helix', 'topology')


# -------------------------------------------------------------- the capability

def annotate_mhc(sequences, helix: str = 'imgt', per_residue: bool = False,
                 mark_regions: bool = False, col_mapping=None,
                 squeeze: bool = False, collapse=None) -> pd.DataFrame:
    """Which MHC gene each sequence is, its IMGT numbering and its groove helices.

    >>> out = annotate_mhc('GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEG'
    ...                    'PEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGY'
    ...                    'HQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRR'
    ...                    'YLENGKETLQRT')        # HLA-A*02:01, alpha1 + alpha2
    >>> out[['gene', 'allele', 'domains', 'helix']].iloc[0].tolist()
    ['A', 'A0201', 'G-ALPHA1, G-ALPHA2', '50-86, 140-176']

    Args:
        sequences: one amino-acid sequence or any collection of them - a list, Series,
            dict of name -> sequence, list of records, or a DataFrame.
        helix: how the groove helices are found, ``'imgt'`` (default) or
            ``'tpm-bench'`` (see :data:`HELIX_METHODS`).
        per_residue: one row per residue instead of one per sequence.
        collapse: one row per *segment* instead: a run of consecutive residues that
            share the value of this per-residue column, e.g. ``'topology'``, or of
            several, e.g. ``['domain', 'region']``. Any of
            :data:`COLLAPSIBLE`; ``True`` means ``['domain', 'region']``. Implies the
            per-residue numbering.
        mark_regions: add ``topology_seq``, a string as long as the sequence with one
            letter per residue - ``S`` signal peptide, ``O`` extracellular, ``M``
            transmembrane, ``I`` cytoplasmic, ``-`` unknown. The same alphabet and
            column :func:`tcrkit.mhc_sequence` uses, estimated here from the allele
            this sequence matched.
        col_mapping: where the sequence lives when the input is a frame or records,
            e.g. ``{'sequence': 'mhc_seq'}``.
        squeeze: return the plain result instead of a frame.

    Returns:
        Indexed like the input. ``sequence`` echoes it, then

        - ``gene`` (``'A'``, ``'DRB1'``, ``'B2M'``, ...), ``allele`` (the closest
          bundled allele, e.g. ``'A0201'``), ``allele_identity`` (the fraction of
          identical residues between your sequence and that allele over the part that
          aligns; 1.0 is an exact match there), ``mhc_class`` (``'I'``, ``'II'``, or
          ``'-'`` for B2M), ``score`` (alignment score per residue) and ``span``
          (1-based range of the MHC part);
        - ``domains``: the G domains found, in sequence order; ``domain_spans``: their
          1-based residue ranges;
        - ``helix``: the 1-based ranges of the groove-helix residues;
        - ``regions``: which topology regions this sequence contains, in order, then
          one column per region - ``signal_peptide``,
          ``extracellular``, ``transmembrane``, ``cytoplasmic`` - giving that
          region's 1-based range in *this* sequence, carried over from the matched
          allele by alignment. A chain that is not exactly a bundled allele, or only
          a fragment of one, is marked from the closest allele;
        - ``error``: why nothing was found, e.g. for a TCR chain or a peptide.

        With ``per_residue=True``: one row per residue, indexed by the input's index
        repeated: ``position`` (1-based), ``aa``, ``gene``, ``allele``, ``domain``,
        ``imgt`` (the IMGT label, ``'-'`` outside the G domains), ``region`` (IMGT's
        strand ``'A'`` ... ``'D'`` or ``'helix'``, 50-92B), ``helix`` (True / False,
        by `helix`) and
        ``topology`` (``'signal'``, ``'extracellular'``, ``'transmembrane'``,
        ``'cytoplasmic'``, or ``'-'``).
        ``attrs['templates']`` holds the IMGT entries each sequence was numbered on.

        With `collapse`: one row per segment, in sequence order, indexed by the input's
        index repeated: the `collapse` column(s), then ``start`` / ``end`` (1-based,
        inclusive) and ``length``, ``imgt_start`` / ``imgt_end`` (the first and last
        IMGT label in the segment, None if it has none), and its ``sequence``.
    """
    if helix not in HELIX_METHODS:
        raise ValueError(f'helix must be one of {HELIX_METHODS}, got {helix!r}')
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
        rec, num = _annotate_one(seq, helix, mark_regions)
        recs.append(rec)
        if per_residue:
            part = num.drop(columns='template')
            part.insert(0, 'position', np.arange(1, len(num) + 1))
            part.insert(2, 'gene', rec.get('gene'))
            part.insert(3, 'allele', rec.get('allele'))
            part.index = [idx] * len(part)
            parts.append(part)
            templates[idx] = '; '.join(t for t in dict.fromkeys(num['template']) if t)
    if per_residue:
        out = (pd.concat(parts) if parts else pd.DataFrame(columns=_RESIDUE_COLUMNS))
        out.attrs['templates'] = templates
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
