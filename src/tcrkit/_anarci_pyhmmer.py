"""In-process, binary-free replacement for ANARCI's hmmscan subprocess.

This module makes ANARCI number sequences using ``pyhmmer`` (HMMER bound as a
C-extension, executed in-process) instead of shelling out to the ``hmmscan``
binary. It is applied by monkeypatching ``anarci.anarci.run_hmmer`` so that the
public ``anarci.run_anarci`` / ``anarci.anarci`` entry points keep working
unchanged -- no ``hmmscan`` on PATH and no conda required.

Provenance
----------
The pyhmmer approach (``run_pyhmmer`` / ``_find_pyhmmer_hsps`` and the pyhmmer
field access in ``_parse_hmmer_query`` / ``_hmm_alignment_to_states``) is derived
from oxpig/ANARCI PR #38 by Daniel Prihoda ("use pyhmmer to drop the hmmscan
binary"), which is NOT in the PyPI ``anarci`` release. ANARCI is BSD-3-Clause
(Oxford Protein Informatics Group); this derivative keeps that license.

The numbering logic itself is ported faithfully (byte-for-byte) from the
INSTALLED ``anarci.anarci`` module's ``_parse_hmmer_query`` /
``_hmm_alignment_to_states`` / ``_domains_are_same`` so that the produced state
vectors -- and therefore every downstream metric -- are identical to the
hmmscan path. Only the per-domain *field access* is swapped from the Biopython
``Hmmer3TextParser`` HSP objects to pyhmmer ``Domain`` / ``Alignment`` objects:

    Biopython HSP field          ->  pyhmmer Domain / Alignment field
    -------------------------------------------------------------------
    hsp.hit_id                   ->  hsp.alignment.hmm_name        (e.g. 'human_H')
    hsp.hit_description          ->  hsp.alignment.hmm_accession
    hsp.evalue                   ->  hsp.i_evalue
    hsp.bitscore                 ->  hsp.score
    hsp.bias                     ->  hsp.bias
    hsp.query_start              ->  hsp.alignment.target_from - 1
    hsp.query_end                ->  hsp.alignment.target_to
    hsp.hit_start (hmm)          ->  hsp.alignment.hmm_from - 1
    hsp.hit_end   (hmm)          ->  hsp.alignment.hmm_to
    hsp.aln_annotation['RF']     ->  derived from hsp.alignment.hmm_sequence
    hsp.aln_annotation['PP']     ->  derived from hsp.alignment.target_sequence

The ``RF``/``PP`` annotation lines that the Biopython parser exposed are not
available from pyhmmer, but the per-column match/insert/delete classification
they encode is recovered exactly from the standard HMMER alignment encoding:
an insert column has ``'.'`` in the model (``hmm_sequence``) line, and a delete
column has ``'-'`` in the target (``target_sequence``) line. We rebuild
equivalent ``reference_string`` ('x' for match, '.' otherwise) and
``state_string`` ('.' for delete, '0' otherwise) marker strings so the rest of
the ported state-vector loop is unchanged.

NB: in pyhmmer >= 0.12 ``hmm_name`` / ``hmm_accession`` / ``Hit.accession`` are
``str``; in older releases they are ``bytes``. ``_text`` normalises both.
"""
from __future__ import annotations

import os
import sys

import pyhmmer

# Pull the constants / helpers ANARCI's own run_hmmer relies on straight from the
# installed package so this shim tracks whatever ANARCI is installed.
from anarci.anarci import (
    HMM_path,
    all_reference_states,
    get_hmm_length,
)


def _text(value):
    """Return a ``str`` from a pyhmmer name/accession that may be ``str`` or ``bytes``.

    pyhmmer >= 0.12 returns ``str`` for ``hmm_name``/``hmm_accession``/``Hit.accession``;
    earlier releases (and oxpig/ANARCI PR #38) return ``bytes``. Handle both.
    """
    if isinstance(value, bytes):
        return value.decode()
    return value


# ---------------------------------------------------------------------------
# Domain dedup + state-vector logic. Ported byte-identically from the installed
# anarci.anarci; only the HSP field access is changed (Biopython HSP -> pyhmmer
# Domain/Alignment) per the table in the module docstring.
# ---------------------------------------------------------------------------
class _HspView:
    """Biopython-HSP-shaped read-only view of a pyhmmer ``Domain``.

    Only the three fields the dedup check touches, mapped per the table in the
    module docstring. pyhmmer ``Domain`` objects are immutable C extensions, so
    the adapter has to wrap rather than annotate them.
    """
    __slots__ = ('_dom',)

    def __init__(self, dom):
        self._dom = dom

    @property
    def hit_id(self):
        return _text(self._dom.alignment.hmm_name)

    @property
    def query_start(self):
        return self._dom.alignment.target_from - 1

    @property
    def query_end(self):
        return self._dom.alignment.target_to


def _domains_are_same(dom1, dom2):
    """
    Check to see if the domains are overlapping.
    @param dom1:
    @param dom2:

    @return: True or False

    Rather than reimplementing the overlap test, delegate to the installed
    ``anarci.anarci._domains_are_same`` through ``_HspView``. The lookup is done
    at call time and on the module attribute, so a monkeypatch of it -- e.g.
    ``tcrkit.anarci._split_chain_types``, which wraps the check to also require a
    matching chain type -- applies on the pyhmmer path exactly as it does on the
    hmmscan path. A local copy of the body would silently ignore such patches.
    """
    _aa = sys.modules['anarci.anarci']
    return _aa._domains_are_same(_HspView(dom1), _HspView(dom2))


def _hmm_alignment_to_states(hsp, n, seq_length, order):
    """
    Take a hit hsp and turn the alignment into a state vector with sequence indices

    Ported from anarci.anarci._hmm_alignment_to_states. The ``reference_string``
    (RF) and ``state_string`` (PP) marker lines the Biopython parser provided are
    reconstructed from the pyhmmer alignment: an insert column is ``'.'`` in the
    HMM (model) line and a delete column is ``'-'`` in the target line. Everything
    after that derivation is identical to the installed implementation.

    The installed ANARCI tagged each Biopython HSP with ``hsp.order`` and read it
    back here; pyhmmer ``Domain`` objects are immutable C extensions and cannot
    carry an extra attribute, so ``order`` is passed explicitly (same value -- the
    domain's index in sequence order -- as ``anarci.anarci`` assigned). This is the
    only behavioural deviation from the installed body and it is value-identical.
    """
    ali = hsp.alignment

    # pyhmmer alignment strings. In the standard HMMER alignment encoding:
    #   - the model line ('hmm_sequence') has '.' at insert-state columns
    #   - the target line ('target_sequence') has '-' at delete-state columns
    hmm_seq = ali.hmm_sequence.upper()
    target_seq = ali.target_sequence

    # Reconstruct the Biopython-style RF / PP annotation lines so the original
    # numbering loop below is reproduced byte-for-byte:
    #   reference_string[i] == 'x'  <=>  match column (else insert)
    #   state_string[i]     == '.'  <=>  delete column
    reference_string = ''.join('.' if h == '.' else 'x' for h in hmm_seq)
    state_string = ''.join('.' if t == '-' else '0' for t in target_seq)

    assert len(reference_string) == len(state_string), "Aligned reference and state strings had different lengths. Don't know how to handle"

    # Extract the start an end points of the hmm states and the sequence
    # These are python indices i.e list[ start:end ] and therefore start will be one less than in the text file
    _hmm_start = ali.hmm_from - 1
    _hmm_end = ali.hmm_to

    _seq_start = ali.target_from - 1
    _seq_end = ali.target_to

    # Extact the full length of the HMM hit
    species, ctype = _text(ali.hmm_name).split('_')
    _hmm_length = get_hmm_length(species, ctype)

    # Handle cases where there are n terminal modifications.
    # In most cases the user is going to want these included in the numbered domain even though they are not 'antibody like' and
    # not matched to the germline. Only allow up to a maximum of 5 unmatched states at the start of the domain
    # Adds a bug here if there is a very short linker between a scfv domains with a modified n-term second domain
    # Thus this is only done for the first identified domain ( hence order attribute on hsp )
    if order == 0 and _hmm_start and _hmm_start < 5:
        n_extend = _hmm_start
        if _hmm_start > _seq_start:
            n_extend = min(_seq_start, _hmm_start - _seq_start)
        state_string = '8' * n_extend + state_string
        reference_string = 'x' * n_extend + reference_string
        _seq_start = _seq_start - n_extend
        _hmm_start = _hmm_start - n_extend

    # Handle cases where the alignment should be extended to the end of the j-element
    # This occurs when there a c-terminal modifications of the variable domain that are significantly different to germline
    # Extension is only made when half of framework 4 has been recognised and there is only one domain recognised.
    if n == 1 and _seq_end < seq_length and (123 < _hmm_end < _hmm_length):  # Extend forwards
        n_extend = min(_hmm_length - _hmm_end, seq_length - _seq_end)
        state_string = state_string + '8' * n_extend
        reference_string = reference_string + 'x' * n_extend
        _seq_end = _seq_end + n_extend
        _hmm_end = _hmm_end + n_extend

    # Generate lists for the states and the sequence indices that are included in this alignment
    hmm_states = all_reference_states[_hmm_start: _hmm_end]
    sequence_indices = list(range(_seq_start, _seq_end))
    h, s = 0, 0  # initialise the current index in the hmm and the sequence

    state_vector = []
    # iterate over the state string (or the reference string)
    for i in range(len(state_string)):
        if reference_string[i] == "x":  # match state
            state_type = "m"
        else:  # insert state
            state_type = "i"

        if state_string[i] == ".":  # overloading if deleted relative to reference. delete_state
            state_type = "d"
            sequence_index = None
        else:
            sequence_index = sequence_indices[s]
        # Store the alignment as the state identifier (uncorrected IMGT annotation) and the index of the sequence

        state_vector.append(((hmm_states[h], state_type), sequence_index))

        # Updates to the indices
        if state_type == "m":
            h += 1
            s += 1
        elif state_type == "i":
            s += 1
        else:  # delete state
            h += 1

    return state_vector


def _parse_hmmer_query(hsps, seq_len, hmmer_species=None):
    """
    Parse the pyhmmer domain hits for a single sequence.

    Ported from anarci.anarci._parse_hmmer_query. The original received a
    Biopython ``query`` object and read ``query.hsps`` / ``query.seq_len``; this
    version receives the already-collected list of pyhmmer ``Domain`` objects and
    the sequence length directly (the hits/length having been gathered in
    ``_find_pyhmmer_hsps`` / ``run_pyhmmer``). Bit-score thresholding happens at
    hit-collection time (matching hmmscan, where only domains above threshold are
    reported), so the ``>= bit_score_threshold`` guards here are unconditional.

    @param hsps: list of pyhmmer Domain objects for one sequence (most significant first irrelevant)
    @param seq_len: length of the query sequence

    The function will identify multiple domains if they have been found and provide the details for the best alignment for each domain.
    This allows the ability to identify single chain fvs and engineered antibody sequences as well as the capability in the future for identifying constant domains.
    """
    hit_table = [['id', 'description', 'evalue', 'bitscore', 'bias',
                  'query_start', 'query_end']]

    # Find the best hit for each domain in the sequence.

    top_descriptions, domains, state_vectors = [], [], []

    if hsps:  # We have some hits
        # If we have specified a species, check to see we have hits for that species
        # Otherwise revert back to using any species
        if hmmer_species:
            hit_correct_species = []
            for hsp in hsps:
                for species in hmmer_species:
                    if _text(hsp.alignment.hmm_name).startswith(species):
                        hit_correct_species.append(hsp)

            if hit_correct_species:
                hsp_list = hit_correct_species
            else:
                print("Limiting hmmer search to species %s was requested but hits did not achieve a high enough bitscore. Reverting to using any species" % (hmmer_species))
                hsp_list = hsps
        else:
            hsp_list = hsps

        for hsp in sorted(hsp_list, key=lambda x: x.i_evalue):  # Iterate over the matches of the domains in order of their e-value (most significant first)
            new = True
            for i in range(len(domains)):  # Check to see if we already have seen the domain
                if _domains_are_same(domains[i], hsp):
                    new = False
                    break
            ali = hsp.alignment
            hit_table.append([_text(ali.hmm_name), _text(ali.hmm_accession), hsp.i_evalue, hsp.score, hsp.bias, ali.target_from - 1, ali.target_to])
            if new:  # It is a new domain and this is the best hit. Add it for further processing.
                domains.append(hsp)
                top_descriptions.append(dict(list(zip(hit_table[0], hit_table[-1]))))  # Add the last added to the descriptions list.

        # Reorder the domains according to the order they appear in the sequence.
        ordering = sorted(list(range(len(domains))), key=lambda x: domains[x].alignment.target_from)
        domains = [domains[_] for _ in ordering]
        top_descriptions = [top_descriptions[_] for _ in ordering]

    ndomains = len(domains)
    for i in range(ndomains):  # If any significant hits were identified parse and align them to the reference state.
        # (installed ANARCI set domains[i].order = i here; pyhmmer Domains are
        #  immutable so `order` is threaded through as the explicit `i` argument.)
        species, chain = top_descriptions[i]["id"].split("_")
        state_vectors.append(_hmm_alignment_to_states(domains[i], ndomains, seq_len, i))  # Alignment to the reference states.
        top_descriptions[i]["species"] = species  # Reparse
        top_descriptions[i]["chain_type"] = chain
        top_descriptions[i]["query_start"] = state_vectors[-1][0][-1]  # Make sure the query_start agree if it was changed

    return hit_table, state_vectors, top_descriptions


# ---------------------------------------------------------------------------
# pyhmmer hit collection (derived from oxpig/ANARCI PR #38). hmmsearch is used
# because pyhmmer has no hmmscan; passing Z=len(hmms) makes the reported
# e-values match what hmmscan would report against the same database.
# ---------------------------------------------------------------------------
def _find_pyhmmer_hsps(sequence_list, hmms, ncpu=None, bit_score_threshold=80):
    alphabet = hmms[0].alphabet
    Z = len(hmms)  # CRITICAL: database size so e-values match hmmscan

    sequences = [
        pyhmmer.easel.TextSequence(sequence=seq, accession=str(i).encode()).digitize(alphabet)
        for i, (sid, seq) in enumerate(sequence_list)
    ]

    # hmmsearch (hmmscan unavailable in pyhmmer) returns a list of hits per HMM;
    # Z=len(hmms) normalises the e-values to match a per-sequence hmmscan run.
    hits_per_hmm = pyhmmer.hmmsearch(
        hmms,
        sequences,
        cpus=0 if ncpu is None else ncpu,
        Z=Z,
    )

    hsps_batch = [[] for _ in range(len(sequence_list))]

    # Distribute the domain hits back to their query sequence (transpose
    # hmmsearch's per-HMM layout to ANARCI's per-sequence layout) and apply the
    # bit-score threshold up front, exactly as hmmscan reporting would.
    for hmm_hits in hits_per_hmm:
        for hit in hmm_hits:
            for hsp in hit.domains:
                if hsp.score >= bit_score_threshold:
                    hsps_batch[int(_text(hit.accession))].append(hsp)

    return hsps_batch


def run_pyhmmer(sequence_list, hmm_database="ALL", ncpu=None, bit_score_threshold=80, hmmer_species=None):
    """
    Run the sequences in sequence list against a precompiled hmm_database, in-process.

    Drop-in replacement for ``anarci.anarci.run_hmmer`` (same return shape: a list
    with one ``(hit_table, state_vectors, top_descriptions)`` tuple per input
    sequence) that uses pyhmmer instead of the hmmscan subprocess.

    @param sequence_list: a list of (name, sequence) tuples. Both are strings
    @param hmm_database: The hmm database to use. Currently, all hmms are in the ALL database.
    @param ncpu: The number of cpu's to allow hmmer to use.
    @param bit_score_threshold: bit-score threshold for a hit to be considered.
    @param hmmer_species: optional list of species to restrict hits to.
    """
    assert hmm_database in ["ALL"], "Unknown HMM database %s" % hmm_database
    hmm_full_path = os.path.join(HMM_path, '{}.hmm'.format(hmm_database))
    if not os.path.exists(hmm_full_path):
        print("ANARCI HMM database not found:", hmm_full_path, file=sys.stderr)
        raise FileNotFoundError(hmm_full_path)

    with pyhmmer.plan7.HMMFile(hmm_full_path) as hmm_file:
        hmms = list(hmm_file)

    hsps_per_sequence = _find_pyhmmer_hsps(
        sequence_list, hmms, ncpu=ncpu, bit_score_threshold=bit_score_threshold
    )
    assert len(sequence_list) == len(hsps_per_sequence), \
        'Unexpected mismatch in returned sequence hits: {} != {}'.format(len(sequence_list), len(hsps_per_sequence))

    results = []
    for (sid, seq), hsps in zip(sequence_list, hsps_per_sequence):
        results.append(_parse_hmmer_query(hsps=hsps, seq_len=len(seq), hmmer_species=hmmer_species))
    return results


_PATCHED = False


def patch_anarci():
    """Monkeypatch ``anarci.anarci.run_hmmer`` to the in-process pyhmmer version.

    Idempotent. After this, ``anarci.run_anarci`` / ``anarci.anarci`` number
    sequences with no ``hmmscan`` subprocess. ``hmmerpath`` / ``ncpu`` arguments
    that ANARCI threads through to ``run_hmmer`` are accepted (``hmmerpath`` is
    simply ignored -- there is no binary to locate).
    """
    global _PATCHED
    if _PATCHED:
        return
    # NB: `anarci/__init__.py` does `from .anarci import anarci`, so the attribute
    # `anarci.anarci` is the *function*, shadowing the submodule. `import anarci.anarci
    # as x` would therefore bind the function, not the module. Pull the real module
    # object out of sys.modules so we patch the right `run_hmmer`.
    import anarci  # ensure the package (and its .anarci submodule) is imported
    _aa = sys.modules["anarci.anarci"]

    def _run_hmmer(sequence_list, hmm_database="ALL", hmmerpath="", ncpu=None,
                   bit_score_threshold=80, hmmer_species=None):
        # hmmerpath is accepted for signature compatibility and ignored: pyhmmer
        # runs in-process, so there is no hmmscan binary to locate.
        return run_pyhmmer(
            sequence_list,
            hmm_database=hmm_database,
            ncpu=ncpu,
            bit_score_threshold=bit_score_threshold,
            hmmer_species=hmmer_species,
        )

    _aa.run_hmmer = _run_hmmer
    # anarci.anarci imports run_hmmer into its own module namespace and calls it
    # unqualified inside `anarci()` and `check_for_j()`, so rebinding the module
    # attribute above is sufficient. Expose marker for tests/debugging.
    _aa._run_hmmer_is_pyhmmer = True
    _PATCHED = True
