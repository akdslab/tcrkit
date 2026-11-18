"""Test-wide settings, for the tests and the doctests: no downloads, and an empty cache.

The suite must not depend on the network or on what happens to be in ``~/.cache/tcrkit``,
so every test runs with ``TCRKIT_NO_DOWNLOAD=1`` and a fresh ``TCRKIT_CACHE``. Without the
MHC Motif Atlas, :func:`tcrkit.annotate_mhc` falls back to the bundled alleles, as it does
offline. Tests of downloaded sources build their own small stand-ins
(``tests/test_motifatlas.py``, ``tests/test_mhc_sources.py``).
"""

import os
import tempfile
import warnings

os.environ['TCRKIT_NO_DOWNLOAD'] = '1'
os.environ['TCRKIT_CACHE'] = tempfile.mkdtemp(prefix='tcrkit_test_cache_')
warnings.filterwarnings('ignore', message='MHC Motif Atlas not available')
