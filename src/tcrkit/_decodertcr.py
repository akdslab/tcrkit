"""DecoderTCR's HLA reference: an opt-in second source of MHC sequences.

Nothing is shipped with tcrkit: the file is downloaded from DecoderTCR's GitHub
repository the first time it is asked for, at a pinned commit and checked against a
pinned SHA-256, then read from a local cache::

    tcrkit.mhc_sequence('HLA-A*02:01', source='decodertcr')    # downloads once
    tcrkit.load_mhc_table(source='decodertcr')                  # the whole table

DecoderTCR stores *pairs*; this splits them into single chains, so each chain can be
asked for on its own (``'DRB1*04:01'``) or as a pair (``paired=True``, or an explicit
``'DRA*01:01__DRB1*04:01'``). The pairs each chain came from are in ``pairs``.

DecoderTCR's sequences are mature chains with no topology. Each one is annotated with
:func:`tcrkit.annotate_mhc` against tcrkit's own (UniProt) reference, which carries the
matched allele's signal peptide / extracellular / transmembrane / cytoplasmic spans
onto it by alignment - so ``region='extracellular'`` cuts the same stretch of the chain
as it does for the bundled sequence, and the IMGT numbering is the same.

DecoderTCR (https://github.com/Biohub/DecoderTCR) is MIT-licensed, Copyright 2026 Chan
Zuckerberg Biohub. Its ``hla_reference.json`` maps 88 HLA pairs (``'A0201__B2M'``,
``'DRA0101__DRB10401'``) to the mature sequences of both chains, ``HLA_a`` / ``HLA_b``:
87 distinct chains. They are fetched by the user rather than redistributed here.

The cache is ``$TCRKIT_CACHE/decodertcr/<commit>/`` (default ``~/.cache/tcrkit``).
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from pathlib import Path

import pandas as pd

#: The DecoderTCR commit the file is taken from (main on 2026-08-04), and the file's
#: SHA-256 at that commit. A different file means the pin was edited, not that GitHub
#: changed underneath.
COMMIT = '3e3d9889d26d79635f674940ddab039c4dc6f9f9'
SHA256 = '9d425bd350b91bd9a595af6ed4137268251e0ca810efd7a959f4faae5d2b59d7'
PATH = 'src/DecoderTCR/reconstruct/data/hla_reference.json'
URL = f'https://raw.githubusercontent.com/Biohub/DecoderTCR/{COMMIT}/{PATH}'
CITATION = ('Lai B, Englund M, Bharanikumar R, Nocedal I, Davariashtiyani A, Perera J, Khan AA. DecoderTCR: Compositional Pretraining and Entropy-Guided Decoding for TCR-pMHC Interactions. ICML (2026). https://openreview.net/pdf?id=yzes8qBM70'
            f' Data: https://github.com/Biohub/DecoderTCR/tree/{COMMIT} (MIT licence).')


def cache_dir() -> Path:
    root = Path(os.environ.get('TCRKIT_CACHE', Path.home() / '.cache' / 'tcrkit'))
    return root / 'decodertcr' / COMMIT


def fetch_decodertcr(force: bool = False) -> Path:
    """The local copy of ``hla_reference.json``, downloading it the first time."""
    path = cache_dir() / 'hla_reference.json'
    if path.exists() and not force:
        if hashlib.sha256(path.read_bytes()).hexdigest() == SHA256:
            return path
        path.unlink()                                  # corrupt or partial: fetch again
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(URL, timeout=120) as resp:
        data = resp.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SHA256:
        raise RuntimeError(f'DecoderTCR hla_reference.json at {COMMIT[:8]} has SHA-256 '
                           f'{digest}, expected {SHA256}; not using it.')
    tmp = path.with_suffix('.part')
    tmp.write_bytes(data)
    tmp.replace(path)
    return path


def load_decodertcr_table(force_download: bool = False) -> pd.DataFrame:
    """DecoderTCR's HLA chains as a tcrkit MHC table, one row per chain.

    Columns: ``allele`` (tcrkit's key, as :func:`tcrkit.normalize_mhc` writes it -
    ``A0201``, ``DRB1*0401``; ``B2M_human`` for B2M), ``sequence`` (the mature chain,
    as DecoderTCR gives it), ``length``, ``gene``, ``mhc_class``, ``species``,
    ``mhc_prefix``, ``source`` (``'decodertcr'``), ``decodertcr_name`` (DecoderTCR's own
    name for the chain), ``pairs`` (the DecoderTCR keys the chain occurs in), the topology span columns of :func:`tcrkit.load_mhc_table`
    carried over from the matched bundled allele (``topology_from`` names it), and
    ``in_detect`` False - gene matching keeps using the bundled reference. ``attrs``
    carries the commit and the citation.
    """
    from .mhc import annotate_mhc
    from .normalize import normalize_mhc

    raw = json.loads(fetch_decodertcr(force_download).read_text())
    chains: dict[str, dict] = {}
    for key, pair in raw.items():
        for name, seq in zip(key.split('__'), (pair['HLA_a'], pair['HLA_b'])):
            entry = chains.setdefault(name, {'name': name, 'sequence': seq, 'pairs': []})
            if entry['sequence'] != seq:
                raise ValueError(f'DecoderTCR gives {name} two different sequences')
            entry['pairs'].append(key)
    table = pd.DataFrame(list(chains.values()))
    norm = normalize_mhc(list(table['name']), errors='coerce')
    table['allele'] = [('B2M_human' if n == 'B2M' else a)
                       for n, a in zip(table['name'], norm['allele'])]
    table['gene'] = ['B2M' if n == 'B2M' else g for n, g in zip(table['name'], norm['gene'])]
    table['mhc_class'] = ['-' if n == 'B2M' else c
                          for n, c in zip(table['name'], norm['mhc_class'])]
    table['species'] = 'human'
    table['mhc_prefix'] = 'HLA'
    table['length'] = table['sequence'].str.len()
    table['source'] = 'decodertcr'
    table['pairs'] = table['pairs'].map(', '.join)

    # topology from tcrkit's own reference, carried over by alignment
    ann = annotate_mhc(dict(zip(table['allele'], table['sequence'])))
    for region, col in (('signal_peptide', 'signal_spans'), ('extracellular', 'outside_spans'),
                        ('transmembrane', 'tmhelix_spans'), ('cytoplasmic', 'inside_spans')):
        table[col] = [_as_spans(v) for v in ann[region]]
    table['mature_start'] = [s[0][1] if s and s[0][0] == 1 else 0
                             for s in map(_parse_spans, table['signal_spans'])]
    # DeepTMHMM's classes: SP+TM, TM, SP (signal peptide only), GLOB (neither)
    table['prediction'] = [
        {(True, True): 'SP+TM', (False, True): 'TM', (True, False): 'SP',
         (False, False): 'GLOB'}[(isinstance(sig, str), isinstance(tm, str))]
        for sig, tm in zip(table['signal_spans'], table['tmhelix_spans'])]
    table['topology_from'] = ['carried over from ' + a if isinstance(a, str) else None
                              for a in ann['allele']]
    table['in_detect'] = False
    table = table.rename(columns={'name': 'decodertcr_name'})
    table = table[['allele', 'gene', 'mhc_class', 'species', 'mhc_prefix', 'sequence',
                   'length', 'source', 'decodertcr_name', 'pairs', 'prediction', 'topology_from',
                   'signal_spans', 'outside_spans', 'tmhelix_spans', 'inside_spans',
                   'mature_start', 'in_detect']]
    table.attrs.update(commit=COMMIT, citation=CITATION)
    return table


def _as_spans(ranges):
    """annotate_mhc's ``'25-307'`` / ``'1-14, 20-30'`` -> ``'[(25, 307)]'``, or None."""
    if not isinstance(ranges, str) or not ranges.strip():
        return None
    return str([tuple(int(x) for x in r.split('-')) for r in ranges.split(', ')])


def _parse_spans(value):
    import ast
    return ast.literal_eval(value) if isinstance(value, str) else []
