"""tcrkit - shared TCR-processing toolkit for the group.

Five capabilities, all re-exported here:

    stitch          build a full TCR nt/aa sequence from V/J calls + CDR3   (Stitchr)
    anarci          annotate V/J/CDR regions of a protein sequence          (ANARCI)
    normalize_tcr   standardize TCR gene/allele names to IMGT               (tidytcells)
    normalize_mhc   standardize HLA/MHC allele names                        (mhcgnomes)
    mhc_sequence    look up an MHC allele's protein sequence                (reference table)
    add_cdr3_junctions  put the IMGT junction residues back on a CDR3    (reference table)

They live in five modules: ``stitch``, ``anarci``, ``normalize`` (allele *names*),
``pmhc`` (an allele to its *sequence*) and ``cdr3_junction``.

:mod:`tcrkit.tcr_pmhc` assembles whole complexes from the lot;
:mod:`tcrkit.metrics` adds repertoire-level measures on top - :func:`diversity`
(richness, entropy, clonality), :func:`morisita_horn` (overlap),
:func:`generation_probability` (OLGA) and :func:`tcrdist` (tcrdist3). The last two need
the optional ``metrics`` extra.

:func:`build_tcr_pmhc` runs all of them end to end over a table of TCR pairs and their
MHC restriction - normalize, complete the junctions, stitch, look up the MHC sequences,
annotate with ANARCI, and check the CDR3 survived the round trip.

One function each. Every one takes a single value *or* any collection - a list, dict,
dict of lists, list of records, Series or DataFrame - detects which you gave, and returns
a DataFrame with your input echoed back beside the result. Nothing assumes a particular
column naming convention: pass ``col_mapping=`` to say where your fields live.

    >>> import tcrkit
    >>> tcrkit.normalize_tcr('TCRBV20S1', species='human')['standardized'][0]
    'TRBV20-1'
    >>> tcrkit.normalize_tcr('TCRBV20S1', species='human', squeeze=True)
    'TRBV20-1'

Pass ``squeeze=True`` when you want the bare value instead of a frame.

The capabilities compose as separate steps. A CDR3 reported without its IMGT junction
residues will not stitch, so repair it first and then stitch::

    df['cdr3'] = tcrkit.add_cdr3_junctions(df, col_mapping={'cdr3': 'cdr3_aa',
                                                      'j': 'j_call'})['cdr3_fixed']
    out = tcrkit.stitch(df, chain='beta')

The same five are available from the command line::

    tcrkit stitch --v TRBV19 --j TRBJ2-7 --cdr3 CASSIRSSYEQYF
    tcrkit --help

Note: importing a name here rebinds it over the submodule of the same name, so
``tcrkit.stitch`` is the *function*. To reach a module's other members use
``from tcrkit.stitch import stitch_pair`` (or ``import tcrkit.stitch as stitch_mod``).
"""

from __future__ import annotations

__version__ = '0.1.0'

from .anarci import CDR_IMGT_LOCS, anarci, read_sequences
from .cdr3_junction import (
    add_cdr3_junctions,
    build_junction_table,
    junction_of,
    load_junction_table,
)
from .metrics import diversity, generation_probability, morisita_horn, tcrdist
from .normalize import genes_match, mhc_chain_of, normalize_mhc, normalize_tcr, to_gene
from .pmhc import bundled_mhc_table_path, load_mhc_table, mhc_sequence
from .tcr_pmhc import build_tcr_pmhc
from .stitch import DEFAULT_CONSTANT, stitch, stitch_pair

__all__ = [
    '__version__',
    # 1. stitch
    'stitch', 'stitch_pair', 'DEFAULT_CONSTANT',
    # 2. anarci
    'anarci', 'read_sequences', 'CDR_IMGT_LOCS',
    # 3. normalize_tcr
    'normalize_tcr', 'to_gene', 'genes_match',
    # 4. normalize_mhc
    'normalize_mhc', 'mhc_chain_of',
    # pMHC sequences
    'mhc_sequence', 'load_mhc_table', 'bundled_mhc_table_path',
    # 5. cdr3 junctions
    'add_cdr3_junctions', 'junction_of', 'load_junction_table',
    'build_junction_table',
    # the lot, end to end
    'build_tcr_pmhc',
    # repertoire metrics (the last two need the 'metrics' extra)
    'diversity', 'morisita_horn', 'generation_probability', 'tcrdist',
]
