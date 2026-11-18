"""Download everything tcrkit fetches on first use, in one go.

    tcrkit.download_data()                       # Stitchr (human, mouse), Motif Atlas, DecoderTCR
    tcrkit.download_data(['motifatlas'])         # just one
    tcrkit download                              # the same, from the command line

Handy before going offline (e.g. onto compute nodes without internet). Each source is
downloaded only if it is not there yet, unless ``force=True``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

__all__ = ['download_data', 'DATA_SOURCES']

#: What :func:`download_data` can fetch.
DATA_SOURCES = ('stitchr', 'motifatlas', 'decodertcr')


def _stitchr(species: str, force: bool) -> dict:
    """IMGT germline data for Stitchr, via its own ``stitchrdl``."""
    from Stitchr import stitchrfunctions as fxn

    folder = Path(fxn.data_dir) / species.upper()
    if (folder / 'TRB.fasta').exists() and not force:
        return {'status': 'cached', 'path': str(folder)}
    exe = (shutil.which('stitchrdl')
           or shutil.which('stitchrdl', path=str(Path(sys.executable).parent)))
    if not exe:
        raise FileNotFoundError('stitchrdl not found; install stitchr in this environment')
    proc = subprocess.run([exe, '-s', species], capture_output=True, text=True)
    if proc.returncode != 0 or not (folder / 'TRB.fasta').exists():
        raise RuntimeError(f'stitchrdl -s {species} failed: '
                           f'{(proc.stderr or proc.stdout).strip()[-500:]}')
    return {'status': 'downloaded', 'path': str(folder)}


def _motifatlas(force: bool) -> dict:
    from . import motifatlas

    had = all((motifatlas.cache_dir() / name).exists() for name in motifatlas.FILES)
    motifatlas.fetch_motifatlas(force=force)
    motifatlas.load_motifatlas_table.cache_clear()
    motifatlas.load_motifatlas_table(force_download=False)  # normalize the names once
    return {'status': 'cached' if had and not force else 'downloaded',
            'path': str(motifatlas.cache_dir())}


def _decodertcr(force: bool) -> dict:
    from . import _decodertcr

    had = (_decodertcr.cache_dir() / 'hla_reference.json').exists()
    path = _decodertcr.fetch_decodertcr(force=force)
    return {'status': 'cached' if had and not force else 'downloaded',
            'path': str(path.parent)}


def download_data(sources=None, species=('human', 'mouse'), force: bool = False
                  ) -> pd.DataFrame:
    """Download the data tcrkit otherwise fetches on first use.

    Args:
        sources: any of :data:`DATA_SOURCES` - ``'stitchr'`` (IMGT germlines for
            :func:`tcrkit.stitch`, through ``stitchrdl``), ``'motifatlas'`` (MHC Motif
            Atlas sequences and ligands, for :func:`tcrkit.annotate_mhc`'s alleles,
            ``source='motifatlas'`` and :func:`tcrkit.peptide_anchors`) and
            ``'decodertcr'`` (DecoderTCR's HLA reference). Default: all three.
        species: the Stitchr species to download.
        force: download again even if already there.

    Returns:
        One row per item: ``source``, ``species`` (Stitchr only), ``status``
        (``'downloaded'``, ``'cached'`` or ``'failed'``), ``path`` and ``error``. A
        failure is reported, not raised, so the other sources still download.
    """
    sources = DATA_SOURCES if sources is None else (
        (sources,) if isinstance(sources, str) else tuple(sources))
    unknown = set(sources) - set(DATA_SOURCES)
    if unknown:
        raise ValueError(f'sources must be among {DATA_SOURCES}, got {sorted(unknown)}')
    species = (species,) if isinstance(species, str) else tuple(species)
    jobs = []
    for source in sources:
        if source == 'stitchr':
            jobs += [(source, sp, lambda sp=sp: _stitchr(sp, force)) for sp in species]
        elif source == 'motifatlas':
            jobs.append((source, None, lambda: _motifatlas(force)))
        else:
            jobs.append((source, None, lambda: _decodertcr(force)))
    rows = []
    for source, sp, run in jobs:
        try:
            rows.append({'source': source, 'species': sp, **run(), 'error': None})
        except Exception as e:                  # report it, keep going
            rows.append({'source': source, 'species': sp, 'status': 'failed',
                         'path': None, 'error': f'{type(e).__name__}: {e}'})
    return pd.DataFrame(rows, columns=['source', 'species', 'status', 'path', 'error'])
