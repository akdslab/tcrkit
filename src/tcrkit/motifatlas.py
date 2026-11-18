"""The MHC Motif Atlas: MHC sequences, their ligands, and what the ligands say.

Nothing is shipped with tcrkit. The four data files are downloaded from
http://mhcmotifatlas.org the first time they are needed, then read from a local cache
(``$TCRKIT_CACHE/motifatlas``, default ``~/.cache/tcrkit``)::

    tcrkit.mhc_sequence('B*51:250', source='motifatlas', region='groove')
    tcrkit.annotate_mhc(seq)['allele']          # the closest Motif Atlas allele
    tcrkit.peptide_anchors('LLFGYPVYV', 'HLA-A*02:01')
    tcrkit.peptide_motif('DRB1*01:01')          # amino-acid frequencies per position

What is in it:

- **MHC sequences.** 13,966 class I alleles (human A/B/C/E/G, 8 mouse H2), each only
  the 181-residue alpha1 + alpha2 groove (from mature residue 2); 865 class II chains,
  full length with their signal peptide (human DR/DQ/DP, cattle, mouse, chicken).
- **Ligands.** About 483,000 class I peptides (8-14 residues) of 135 alleles and
  600,000 class II peptides of 95 alleles or alpha/beta pairs, each class II one with
  its 9-residue binding core. The class II ``__reverse`` sets (peptides bound in the
  reverse orientation) are left out.

Allele names are normalized once, when the table is first built (``allele``, as
:func:`tcrkit.normalize_mhc` writes it); the Motif Atlas's own spelling is kept in
``motifatlas_name``. The class II chains get their membrane topology from tcrkit's own
reference by alignment, as DecoderTCR's do, so a region means the same stretch of the
chain as for the bundled sequences. Class I has only the groove, so only
``region='groove'`` is available for it.

Set ``TCRKIT_NO_DOWNLOAD=1`` to never download; whatever needs the data then falls back
or says so in ``error``.

Data: Tadros DM, Eggenschwiler S, Racle J, Gfeller D. The MHC Motif Atlas: a database of
MHC binding specificities and ligands. Nucleic Acids Research 51(D1), D428-D437 (2023).
https://doi.org/10.1093/nar/gkac965
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import time
import urllib.request
import warnings
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from ._util import ERROR_COL, collect, squeeze as _squeeze

__all__ = ['peptide_anchors', 'peptide_motif', 'plot_motif', 'load_motifatlas_table',
           'load_motifatlas_peptides', 'fetch_motifatlas', 'closest_allele', 'CITATION',
           'MIN_LIGANDS', 'ANCHOR_FRACTION']

CITATION = ('Tadros DM, Eggenschwiler S, Racle J, Gfeller D. The MHC Motif Atlas: a '
            'database of MHC binding specificities and ligands. Nucleic Acids Research '
            '51(D1), D428-D437 (2023). https://doi.org/10.1093/nar/gkac965')

BASE_URL = 'http://mhcmotifatlas.org/data/'
#: The files, by the name they are cached under.
FILES = {
    'sequences_I.txt': 'classI/MHC_I_sequences.txt',
    'sequences_II.txt': 'classII/MHC_II_sequences.txt',
    'peptides_I.txt': 'classI/MS/Peptides/all_peptides.txt',
    'peptides_II.txt': 'classII/MS/Peptides/all_peptides.txt',
}

#: An allele (and, for class I, peptide length) needs at least this many ligands for a
#: motif; fewer, and :func:`peptide_anchors` reports an error instead of guessing.
MIN_LIGANDS = 50

#: A position is an anchor when its information content is at least this fraction of
#: the motif's highest.
ANCHOR_FRACTION = 0.5

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

_SPECIES = {'Bos sp.': 'cattle', 'Bos taurus': 'cattle', 'Gallus gallus': 'chicken'}
_TABLE_VERSION = 2                     # bump when the processed table's layout changes


# ===================================================================== the files

def cache_dir() -> Path:
    root = Path(os.environ.get('TCRKIT_CACHE', Path.home() / '.cache' / 'tcrkit'))
    return root / 'motifatlas'


def downloads_allowed() -> bool:
    return os.environ.get('TCRKIT_NO_DOWNLOAD', '').strip().lower() in ('', '0', 'false', 'no')


def fetch_motifatlas(force: bool = False, files=None) -> dict[str, Path]:
    """The local copies of the Motif Atlas files, downloading the missing ones.

    Returns ``{name: path}``. Each download is recorded in ``manifest.json`` next to the
    files (URL, date, size, SHA-256): the Motif Atlas has no versions, so that is what
    says which copy you used. ``force=True`` downloads again.
    """
    folder = cache_dir()
    names = list(FILES) if files is None else list(files)
    out, manifest_path = {}, folder / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for name in names:
        path = folder / name
        if path.exists() and not force:
            out[name] = path
            continue
        if not downloads_allowed():
            raise RuntimeError(f'MHC Motif Atlas {name} is not in {folder}, and '
                               'TCRKIT_NO_DOWNLOAD is set.')
        folder.mkdir(parents=True, exist_ok=True)
        url = BASE_URL + FILES[name]
        with urllib.request.urlopen(url, timeout=300) as resp:
            data = resp.read()
        if not data.startswith(b'Allele'):
            raise RuntimeError(f'{url} did not return a Motif Atlas table.')
        tmp = path.with_suffix('.part')
        tmp.write_bytes(data)
        tmp.replace(path)
        manifest[name] = {'url': url, 'downloaded': time.strftime('%Y-%m-%d'),
                          'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
        manifest_path.write_text(json.dumps(manifest, indent=1))
        out[name] = path
    return out


# ===================================================================== allele names

_SPLIT = re.compile(r'^(?:(?P<prefix>BoLA|Gaga)_)?(?P<gene>[A-Z]+\d*)_(?P<a>\d+)_(?P<b>\d+)$')


def _ligand_name(part: str) -> str:
    """A ligand file's ``'DRB1_01_01'`` / ``'BoLA_DRB3_001_01'`` -> ``'DRB1*01:01'`` /
    ``'BoLA-DRB3*001:01'``; anything else unchanged (``'A0201'``, ``'H2_Aa_b'``)."""
    m = _SPLIT.match(part)
    if not m:
        return part
    prefix = f"{m['prefix']}-" if m['prefix'] else ''
    return f"{prefix}{m['gene']}*{m['a']}:{m['b']}"


def _normalized(names) -> dict:
    """``{name: (allele, gene, mhc_class, species)}``, through :func:`normalize_mhc`."""
    from .normalize import normalize_mhc

    names = list(dict.fromkeys(names))
    norm = normalize_mhc(names, errors='ignore').reindex(
        columns=['allele', 'gene', 'mhc_class', 'species'])
    out = {}
    for name, a, g, c, s in zip(names, norm['allele'], norm['gene'], norm['mhc_class'],
                                norm['species']):
        ok = isinstance(a, str)
        out[name] = (a if ok else name, g if ok else None, c if ok else None,
                     _SPECIES.get(s, s) if isinstance(s, str) else None)
    return out


def _ligand_key(name: str, names: dict) -> str:
    """The allele key of a ligand set: ``'DRB1*0101'``, ``'DPA1*0103__DPB1*0101'``."""
    return '__'.join(names[_ligand_name(p)][0] for p in name.split('__'))


# ===================================================================== the sequences

def _spans(value: str | None) -> str | None:
    """annotate_mhc's ``'25-307'`` / ``'1-14, 20-30'`` -> ``'[(25, 307)]'``, or None."""
    if not isinstance(value, str) or not value.strip():
        return None
    pairs = []
    for part in value.split(','):
        a, _, b = part.strip().partition('-')
        pairs.append((int(a), int(b or a)))
    return str(pairs)


def _build_table(paths: dict) -> pd.DataFrame:
    """The processed sequence table, from the raw files: names normalized, topology and
    groove spans carried over. Slow (it annotates the class II chains); run once."""
    from .mhc import annotate_mhc

    s1 = pd.read_csv(paths['sequences_I.txt'], sep='\t').assign(mhc_group='I')
    s2 = pd.read_csv(paths['sequences_II.txt'], sep='\t').assign(mhc_group='II')
    raw = pd.concat([s1, s2], ignore_index=True)
    raw.columns = ['motifatlas_name', 'sequence', 'mhc_group']
    # 14 cattle chains carry alignment gaps ('-'): the sequence is what is left
    raw['sequence'] = raw['sequence'].str.strip().str.upper().str.replace('-', '', regex=False)
    names = _normalized(raw['motifatlas_name'])
    raw['allele'] = [names[n][0] for n in raw['motifatlas_name']]
    raw['gene'] = [names[n][1] for n in raw['motifatlas_name']]
    raw['mhc_class'] = [names[n][2] or g for n, g in zip(raw['motifatlas_name'],
                                                       raw['mhc_group'])]
    raw['species'] = [names[n][3] for n in raw['motifatlas_name']]
    raw = raw.drop_duplicates('allele', keep='first').reset_index(drop=True)
    raw['length'] = raw['sequence'].str.len()
    raw['source'] = 'motifatlas'

    # which alleles have ligands (preferred when several alleles are equally close)
    ligand = set()
    for name in ('peptides_I.txt', 'peptides_II.txt'):
        if name in paths:
            sets = pd.read_csv(paths[name], sep='\t', usecols=['Allele'])['Allele'].unique()
            parts = [_ligand_name(p) for s in sets for p in s.split('__') if p != 'reverse']
            ligand |= {a for a, *_ in _normalized(parts).values()}
    raw['has_ligands'] = raw['allele'].isin(ligand)

    # class I: the groove alone, all of it extracellular; class II: full chains, their
    # topology and groove carried over from the bundled reference
    one = raw['mhc_group'] == 'I'
    raw['covers'] = np.where(one, 'groove', 'full')
    for col in ('signal_spans', 'outside_spans', 'tmhelix_spans', 'inside_spans',
                'groove_spans', 'topology_from'):
        raw[col] = None
    raw.loc[one, 'outside_spans'] = [str([(1, n)]) for n in raw.loc[one, 'length']]
    raw.loc[one, 'groove_spans'] = raw.loc[one, 'outside_spans']
    two = raw.index[~one]
    if len(two):
        ann = annotate_mhc(dict(zip(raw.loc[two, 'allele'], raw.loc[two, 'sequence'])),
                           allele_source='uniprot')
        for region, col in (('signal_peptide', 'signal_spans'),
                            ('extracellular', 'outside_spans'),
                            ('transmembrane', 'tmhelix_spans'),
                            ('cytoplasmic', 'inside_spans')):
            raw.loc[two, col] = [_spans(v) for v in ann[region]]
        raw.loc[two, 'groove_spans'] = [_groove(v) for v in ann['domain_spans']]
        raw.loc[two, 'topology_from'] = ['carried over from ' + a if isinstance(a, str)
                                         else None for a in ann['allele']]
    raw['mature_start'] = [ast_first_end(v) for v in raw['signal_spans']]
    raw['in_detect'] = False
    cols = ['allele', 'motifatlas_name', 'gene', 'mhc_class', 'species', 'sequence',
            'length', 'source', 'covers', 'has_ligands', 'topology_from', 'signal_spans',
            'outside_spans', 'tmhelix_spans', 'inside_spans', 'groove_spans',
            'mature_start', 'in_detect']
    return raw[cols]


def _groove(domain_spans) -> str | None:
    """annotate_mhc's G-domain spans ``'30-122'`` -> the groove, first to last residue."""
    if not isinstance(domain_spans, str) or not domain_spans.strip():
        return None
    ends = [int(x) for x in re.findall(r'\d+', domain_spans)]
    return str([(min(ends), max(ends))])


def ast_first_end(value) -> int:
    """The signal peptide's length from its span list, 0 if it has none."""
    m = re.match(r'\[\(1, (\d+)\)', value) if isinstance(value, str) else None
    return int(m.group(1)) if m else 0


@lru_cache(maxsize=1)
def load_motifatlas_table(force_download: bool = False) -> pd.DataFrame:
    """The Motif Atlas MHC sequences as a tcrkit MHC table, one row per allele.

    Columns: ``allele`` (normalized, as :func:`tcrkit.normalize_mhc` writes it),
    ``motifatlas_name`` (the Motif Atlas's spelling), ``gene``, ``mhc_class``,
    ``species``, ``sequence``, ``length``, ``source`` (``'motifatlas'``), ``covers``
    (``'groove'`` for class I, which holds only alpha1 + alpha2; ``'full'`` for class
    II), ``has_ligands``, the topology span columns of :func:`tcrkit.load_mhc_table`
    plus ``groove_spans``, ``topology_from`` and ``mature_start``. Built once from the
    downloaded files (normalizing every name and annotating the class II chains takes a
    minute) and cached as ``motifatlas/table.csv``. ``attrs`` has the citation.
    """
    paths = fetch_motifatlas(force=force_download)
    cached = cache_dir() / 'table.csv'
    stamp = cache_dir() / 'table.version'
    fresh = (cached.exists() and stamp.exists() and not force_download
             and stamp.read_text().strip() == str(_TABLE_VERSION)
             and all(cached.stat().st_mtime >= p.stat().st_mtime for p in paths.values()))
    if fresh:
        table = pd.read_csv(cached, dtype={'mhc_class': str})
    else:
        table = _build_table(paths)
        table.to_csv(cached, index=False)
        stamp.write_text(str(_TABLE_VERSION))
    table.attrs['citation'] = CITATION
    return table


@lru_cache(maxsize=1)
def available() -> bool:
    """Whether the Motif Atlas table can be had: cached, or downloadable now."""
    try:
        load_motifatlas_table()
        return True
    except Exception:                                     # offline, or no permission
        return False


# ===================================================================== closest allele

@lru_cache(maxsize=1)
def _by_gene() -> dict:
    """``{(gene, species): (alleles, sequences, has_ligands, k-mer sets)}``."""
    table = load_motifatlas_table()
    out = {}
    for (gene, species), g in table.groupby(['gene', 'species'], dropna=False):
        order = g.sort_values(['has_ligands', 'allele'], ascending=[False, True])
        seqs = order['sequence'].tolist()
        out[(gene, species)] = (order['allele'].tolist(), seqs,
                                order['has_ligands'].tolist(), [_kmers(s) for s in seqs])
    return out


def _kmers(seq: str, k: int = 4) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


#: How many of the k-mer-closest alleles are aligned when nothing matches exactly.
_SHORTLIST = 10


def closest_allele(sequence: str, gene: str, species: str | None = None):
    """The Motif Atlas allele of `gene` closest to `sequence`.

    Only that gene's alleles are searched (the gene comes from
    :func:`tcrkit.mhc.best_mhc`). An allele whose sequence `sequence` contains, or that
    contains `sequence`, is identical to it; otherwise the :data:`_SHORTLIST` alleles
    sharing the most 4-mers are aligned and the most identical wins. Of equally close
    alleles one with ligands comes first, then the lowest name.

    Returns:
        ``(allele, identity, ties)`` - ``ties`` is how many other alleles are equally
        close - or None when the Motif Atlas has no allele of that gene and species.
    """
    from .mhc import _match_allele

    index = _by_gene()
    entry = index.get((gene, species))
    if entry is None and species is None:
        entry = next((v for (g, _), v in index.items() if g == gene), None)
    if entry is None:
        return None
    alleles, seqs, _, kmers = entry
    exact = [i for i, s in enumerate(seqs) if s in sequence or sequence in s]
    if exact:
        return alleles[exact[0]], 1.0, len(exact) - 1
    query = _kmers(sequence)
    shared = np.fromiter((len(query & k) for k in kmers), dtype=int, count=len(kmers))
    short = np.argsort(-shared, kind='stable')[:_SHORTLIST]
    scored = [(round(_match_allele(sequence, seqs[i])[2], 3), -int(r), int(i))
              for r, i in enumerate(short)]
    best = max(scored)
    ties = sum(1 for s in scored if s[0] == best[0]) - 1
    return alleles[best[2]], best[0], ties


# ===================================================================== ligands

@lru_cache(maxsize=1)
def load_motifatlas_peptides() -> pd.DataFrame:
    """The Motif Atlas ligands: ``allele`` (the normalized key, ``'A0201'``,
    ``'DRB1*0101'``, ``'DPA1*0103__DPB1*0101'``), ``motifatlas_name``, ``mhc_class``
    (``'I'`` / ``'II'``), ``peptide``, ``length`` and, for class II, ``core``. The
    ``__reverse`` class II sets are left out."""
    paths = fetch_motifatlas(files=['peptides_I.txt', 'peptides_II.txt'])
    p1 = pd.read_csv(paths['peptides_I.txt'], sep='\t').assign(mhc_class='I', Core=None)
    p2 = pd.read_csv(paths['peptides_II.txt'], sep='\t').assign(mhc_class='II')
    pep = pd.concat([p1, p2], ignore_index=True)
    pep.columns = [c.lower() for c in pep.columns]
    pep = pep[~pep['allele'].str.endswith('__reverse')]
    sets = pep['allele'].unique()
    names = _normalized({_ligand_name(p) for s in sets for p in s.split('__')})
    key = {s: _ligand_key(s, names) for s in sets}
    pep = pep.rename(columns={'allele': 'motifatlas_name'})
    pep.insert(0, 'allele', pep['motifatlas_name'].map(key))
    pep['length'] = pep['peptide'].str.len()
    return pep[['allele', 'motifatlas_name', 'mhc_class', 'peptide', 'length',
                'core']].reset_index(drop=True)


@lru_cache(maxsize=2)
def _background(mhc_class: str) -> np.ndarray:
    """Amino-acid frequencies over every ligand of a class: the motifs' background."""
    pep = load_motifatlas_peptides()
    text = ''.join(pep.loc[pep['mhc_class'] == mhc_class, 'peptide'])
    counts = np.array([text.count(a) for a in AMINO_ACIDS], dtype=float) + 1
    return counts / counts.sum()


def _resolve_set(allele: str) -> tuple[str | None, str | None]:
    """The ligand set an allele name means: ``(key, error)``."""
    from .normalize import normalize_mhc

    sets = load_motifatlas_peptides()['allele'].unique()
    parts = [p for p in re.split(r'__|/', str(allele)) if p]
    norm = normalize_mhc(parts, errors='ignore').reindex(columns=['allele'])['allele'].tolist()
    if any(not isinstance(a, str) for a in norm):
        return None, f'{allele!r} is not an allele name'
    key = '__'.join(a for a in norm if not a.startswith('DRA'))   # DR sets omit DRA
    if key in sets:
        return key, None
    within = [s for s in sets if all(a in s.split('__') for a in key.split('__'))]
    if len(within) == 1:
        return within[0], None
    if within:
        return None, (f'{allele!r} is in {len(within)} Motif Atlas pairs '
                      f"({', '.join(within[:4])}{', ...' if len(within) > 4 else ''}); "
                      'name the pair')
    return None, f'no Motif Atlas ligands for {key}'


def _motif(key: str, length: int | None):
    """``(p, ic, n, mhc_class)``: the background-normalized frequencies (positions x
    amino acids), information content per position, and how many ligands made them."""
    pep = load_motifatlas_peptides()
    rows = pep[pep['allele'] == key]
    cls = rows['mhc_class'].iat[0]
    if cls == 'II':
        seqs = rows['core'].dropna()
        seqs = seqs[seqs.str.len() == 9]
    else:
        seqs = rows.loc[rows['length'] == length, 'peptide']
    n = len(seqs)
    if n == 0:
        return None, None, 0, cls
    width = len(seqs.iat[0])
    counts = np.zeros((width, 20))
    lookup = {a: i for i, a in enumerate(AMINO_ACIDS)}
    for s in seqs:
        for pos, a in enumerate(s):
            j = lookup.get(a)
            if j is not None:
                counts[pos, j] += 1
    bg = _background(cls)
    freq = (counts + bg) / (counts.sum(axis=1, keepdims=True) + 1)   # one pseudo-ligand
    p = freq / bg
    p = p / p.sum(axis=1, keepdims=True)
    ic = math.log2(20) + (p * np.log2(p)).sum(axis=1)
    return p, ic, n, cls


def peptide_motif(allele: str, length: int = 9) -> pd.DataFrame:
    """The binding motif of an allele, from its Motif Atlas ligands.

    One row per position (1-based), one column per amino acid: the frequency at that
    position divided by its frequency in all ligands of the class, renormalized to 1 -
    the Motif Atlas's own normalization, with the class's ligands standing in for the
    proteome. ``ic`` is the information content of each position, in bits (0 = no
    preference, log2(20) = 4.32 = one amino acid only). Class I motifs depend on the
    peptide `length`; class II ones are of the 9-residue binding core, any length.
    ``attrs`` has ``n_ligands``, ``mhc_class`` and ``allele``.
    """
    key, err = _resolve_set(allele)
    if err:
        raise ValueError(err)
    p, ic, n, cls = _motif(key, length)
    if p is None:
        raise ValueError(f'no {length}-mer ligands for {key}')
    out = pd.DataFrame(p, columns=list(AMINO_ACIDS), index=range(1, len(p) + 1))
    out.index.name = 'position'
    out['ic'] = ic
    out.attrs.update(n_ligands=n, mhc_class=cls, allele=key)
    return out


def plot_motif(allele: str, length: int = 9, ax=None):
    """A sequence logo of :func:`peptide_motif`: letter heights = frequency x
    information content, as on the Motif Atlas. Needs logomaker
    (``pip install 'tcrkit[logos]'``). Returns the matplotlib Axes."""
    try:
        import logomaker
    except ImportError as e:
        raise ImportError("plot_motif needs logomaker: pip install 'tcrkit[logos]'") from e
    motif = peptide_motif(allele, length)
    heights = motif[list(AMINO_ACIDS)].mul(motif['ic'], axis=0)
    logo = logomaker.Logo(heights, ax=ax, color_scheme='chemistry')
    kind = 'core' if motif.attrs['mhc_class'] == 'II' else f'{length}-mer'
    logo.ax.set_title(f"{motif.attrs['allele']} {kind} ({motif.attrs['n_ligands']} ligands)")
    logo.ax.set_ylabel('bits')
    return logo.ax


def _best_core(peptide: str, p: np.ndarray, bg: np.ndarray) -> int:
    """The 0-based start of the 9-mer of `peptide` that fits the core motif best."""
    lookup = {a: i for i, a in enumerate(AMINO_ACIDS)}
    odds = np.log2(np.clip(p, 1e-9, None))
    best, start = -np.inf, 0
    for s in range(len(peptide) - p.shape[0] + 1):
        score = sum(odds[k, lookup[a]] for k, a in enumerate(peptide[s:s + p.shape[0]])
                    if a in lookup)
        if score > best:
            best, start = score, s
    return start


def peptide_anchors(peptides, alleles=None, col_mapping=None,
                    squeeze: bool = False) -> pd.DataFrame:
    """The probable anchor positions of peptides on an allele, from its ligands.

    For class I, the allele's ligands of the same length give a motif; for class II, the
    binding cores of its ligands do, and the peptide's core is the 9-mer that fits it
    best (or the recorded one, if the peptide is itself a Motif Atlas ligand). Anchors
    are the motif positions whose information content is at least
    :data:`ANCHOR_FRACTION` of the highest - typically P2 and the last position of a
    class I peptide, P1, P4, P6 and P9 of a class II core. A statistical guess from the
    ligands, not a structure: it needs :data:`MIN_LIGANDS` ligands, else ``error``.

    >>> peptide_anchors('LLFGYPVYV', 'HLA-A*02:01')['anchors'].iat[0]   # doctest: +SKIP
    'P2, P9'

    Args:
        peptides, alleles: a peptide and its allele (any spelling; a class II pair as
            ``'DQA1*01:02__DQB1*06:02'``), or two parallel collections, or a frame /
            records with ``peptide`` and ``allele`` passed as `peptides`.
        col_mapping: e.g. ``{'peptide': 'epitope', 'allele': 'mhc'}``.
        squeeze: return the plain result instead of a frame.

    Returns:
        Indexed like the input: ``peptide``, ``allele``, then ``motif_allele`` (the
        ligand set used), ``mhc_class``, ``n_ligands``, ``core`` / ``core_start``
        (class II; 1-based in the peptide), ``anchors`` (motif positions, ``'P2, P9'``),
        ``anchor_positions`` (1-based in the peptide), ``anchor_residues`` (the
        peptide's residues there), ``anchor_ic`` (their information content, bits) and
        ``error``.
    """
    inp = collect(('peptide', 'allele'), (peptides, alleles), col_mapping=col_mapping)
    cols = ['motif_allele', 'mhc_class', 'n_ligands', 'core', 'core_start', 'anchors',
            'anchor_positions', 'anchor_residues', 'anchor_ic', ERROR_COL]
    pep_table = None
    rows = []
    for idx, pep, allele in zip(inp.index, inp['peptide'], inp['allele']):
        row = dict.fromkeys(cols)
        pep = re.sub(r'\s+', '', pep).upper() if isinstance(pep, str) else ''
        try:
            if not re.fullmatch(r'[A-Z]+', pep):
                raise ValueError('not a peptide sequence')
            key, err = _resolve_set(allele)
            if err:
                raise ValueError(err)
            if pep_table is None:
                pep_table = load_motifatlas_peptides()
            p, ic, n, cls = _motif(key, len(pep))
            row.update(motif_allele=key, mhc_class=cls, n_ligands=n)
            if n < MIN_LIGANDS:
                what = 'cores' if cls == 'II' else f'{len(pep)}-mer ligands'
                raise ValueError(f'{n} {what} for {key}, needs {MIN_LIGANDS}')
            offset = 0
            if cls == 'II':
                if len(pep) < 9:
                    raise ValueError('a class II peptide needs at least 9 residues')
                known = pep_table[(pep_table['allele'] == key) & (pep_table['peptide'] == pep)]
                core = known['core'].dropna()
                offset = (pep.find(core.iat[0]) if len(core) and core.iat[0] in pep
                          else _best_core(pep, p, _background('II')))
                row.update(core=pep[offset:offset + 9], core_start=offset + 1)
            top = ic.max()
            where = [k for k in range(len(ic)) if ic[k] >= ANCHOR_FRACTION * top]
            row.update(anchors=', '.join(f'P{k + 1}' for k in where),
                       anchor_positions=', '.join(str(offset + k + 1) for k in where),
                       anchor_residues=', '.join(pep[offset + k] for k in where),
                       anchor_ic=', '.join(f'{ic[k]:.2f}' for k in where))
        except (ValueError, RuntimeError) as e:
            row[ERROR_COL] = str(e)
        rows.append(row)
    out = pd.concat([inp[['peptide', 'allele']],
                     pd.DataFrame(rows, index=inp.index, columns=cols)], axis=1)
    return _squeeze(out, ('peptide', 'allele')) if squeeze else out


# ===================================================================== annotate_mhc

_warned = False


def allele_for(sequence: str, gene: str, species: str | None):
    """What :func:`tcrkit.annotate_mhc` reports by default: :func:`closest_allele`, or
    None (with one warning per session) when the Motif Atlas cannot be had."""
    global _warned
    if not available():
        if not _warned:
            warnings.warn('MHC Motif Atlas not available (offline, or TCRKIT_NO_DOWNLOAD '
                          "set): alleles come from the bundled reference "
                          "(allele_source='uniprot').", stacklevel=3)
            _warned = True
        return None
    return closest_allele(sequence, gene, species)


def species_of(allele: str) -> str | None:
    """The species of a Motif Atlas allele, if the table is already cached."""
    if not (cache_dir() / 'table.csv').exists():
        return None
    try:
        table = load_motifatlas_table()
    except Exception:
        return None
    hit = table.loc[table['allele'] == allele, 'species']
    return hit.iat[0] if len(hit) else None
