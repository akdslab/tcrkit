"""Smoke tests: one per capability, on trivial hardcoded input.

A fast sanity check that each wrapper is wired to its backend and still produces the
values it always has - not a validation suite. Run with::

    uv run pytest                 # or: pytest

The expected values are not invented: they were captured from the original
``immune-tcr-data`` implementations (Stitchr 1.1.3.1 / tidytcells 2.2.0) and confirmed
to match tcrkit's output, so a regression here means behaviour drifted from the source
these were ported from.
"""

import pandas as pd
import pytest

import tcrkit
from tcrkit.pmhc import bundled_mhc_table_path, load_mhc_table

BETA_AA = ('MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGTTAKGEVPDGYNVSR'
           'LKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL')
ALPHA_AA = ('AQKITQTQPGMFVQEKEAVTLDCTYDTSDPSYGLFWYKQPSSGEMIFLIYQGSYDQQNATEGRYSLNFQ'
            'KARKSANLVISASQLGDSAMYFCAMSQLNSGNTPLVFGKGTRLSVIA')


# ---------------------------------------------------------------------------
# the shared contract: one function per capability, any input shape, a frame out
# ---------------------------------------------------------------------------

def test_all_five_capabilities_are_exported_at_the_top_level():
    for name in ('stitch', 'anarci', 'normalize_tcr', 'normalize_mhc',
                 'add_cdr3_junctions', 'mhc_sequence'):
        assert callable(getattr(tcrkit, name)), name


def test_no_arity_variants_remain():
    """The whole point of the unified API: no *_many / *_frame siblings."""
    for gone in ('anarci_many', 'normalize_tcr_many', 'normalize_tcr_frame',
                 'normalize_mhc_many', 'fix_junctions', 'fix_junction_frame',
                 'stitch_frame', 'stitch_batch', 'stitch_batch_single',
                 'fix_junction', 'load_table', 'build_table'):
        assert not hasattr(tcrkit, gone), f'{gone} should no longer be public'


def test_legacy_compat_modules_are_gone():
    for gone in ('util_stitchr', 'util_thimble', 'util_anarci', 'util_allele',
                 'util_cdr3'):
        assert not hasattr(tcrkit, gone), f'{gone} should have been removed'


def test_the_two_nomenclature_capabilities_share_one_module():
    from tcrkit import normalize

    for name in ('normalize_tcr', 'normalize_mhc', 'to_gene', 'genes_match',
                 'mhc_chain_of'):
        assert hasattr(normalize, name), name


@pytest.mark.parametrize('capability, single, many', [
    ('normalize_tcr', ('TCRBV20S1',), (['TCRBV20S1', 'TRAV23/6'],)),
    ('normalize_mhc', ('A*0201',), (['A*0201', 'H2-Kb'],)),
    ('add_cdr3_junctions', ('ASSYSGNTEAF', 'TRBJ1-1'),
     (['ASSYSGNTEAF', 'AVMDSNYQLI'], ['TRBJ1-1', 'TRAJ33'])),
])
def test_every_capability_returns_a_frame_for_one_and_for_many(capability, single, many):
    func = getattr(tcrkit, capability)
    one, lots = func(*single), func(*many)
    assert isinstance(one, pd.DataFrame) and len(one) == 1
    assert isinstance(lots, pd.DataFrame) and len(lots) == 2


def test_input_shapes_all_reach_the_same_answer():
    """A value, a list, a dict of columns, records and a frame must agree."""
    expected = ['CASSYSGNTEAFF', 'CAVMDSNYQLIW']
    cdr3s, js = ['ASSYSGNTEAF', 'AVMDSNYQLI'], ['TRBJ1-1', 'TRAJ33']

    shapes = [
        tcrkit.add_cdr3_junctions(cdr3s, js),                                   # parallel
        tcrkit.add_cdr3_junctions({'cdr3': cdr3s, 'j': js}),                    # dict of cols
        tcrkit.add_cdr3_junctions([{'cdr3': c, 'j': j} for c, j in zip(cdr3s, js)]),  # records
        tcrkit.add_cdr3_junctions(pd.DataFrame({'cdr3': cdr3s, 'j': js})),      # frame
        tcrkit.add_cdr3_junctions(pd.DataFrame({'x': cdr3s, 'y': js}),          # frame + mapping
                            col_mapping={'cdr3': 'x', 'j': 'y'}),
    ]
    for got in shapes:
        assert list(got['cdr3_fixed']) == expected


def test_index_is_preserved_so_results_assign_straight_back():
    df = pd.DataFrame({'cdr3': ['ASSYSGNTEAF', 'AVMDSNYQLI'],
                       'j': ['TRBJ1-1', 'TRAJ33']}, index=['x', 'y'])
    out = tcrkit.add_cdr3_junctions(df)
    assert list(out.index) == ['x', 'y']
    df['fixed'] = out['cdr3_fixed']          # aligns on the index, no .values needed
    assert list(df['fixed']) == ['CASSYSGNTEAFF', 'CAVMDSNYQLIW']


def test_squeeze_gives_the_bare_value():
    assert tcrkit.normalize_tcr('TCRBV20S1', species='human', squeeze=True) == 'TRBV20-1'
    assert list(tcrkit.normalize_tcr(['TCRBV20S1', 'TRAV23/6'], species='human',
                                     squeeze=True)) == ['TRBV20-1', 'TRAV23/DV6']
    assert tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1',
                               squeeze=True) == 'CASSYSGNTEAFF'
    assert isinstance(tcrkit.normalize_mhc('A*0201', squeeze=True), dict)


# ---------------------------------------------------------------------------
# 1. stitch
# ---------------------------------------------------------------------------

def test_stitch_beta():
    out = tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF',
                        chain='beta', species='HUMAN')
    r = out.iloc[0]
    assert r['error'] == '-'
    assert len(r['sequence']) == 924
    assert r['sequence_aa'].startswith('MSNQVLCCVVLCFLGANTVD')
    # the CDR3 that went in must be readable in the translation that came out
    assert 'CASSIRSSYEQYF' in r['sequence_aa']


def test_stitch_alpha():
    r = tcrkit.stitch('TRAV12-2', 'TRAJ23', 'CAVNTGGGNKLTF',
                      chain='alpha', species='HUMAN').iloc[0]
    assert r['error'] == '-'
    assert len(r['sequence']) == 819
    assert r['sequence_aa'].startswith('MKSLRVLLVILWLQLSWVWS')


def test_stitch_does_not_repair_the_cdr3_itself():
    """Repairing the junction is fix_junction's job, run as a step before stitching."""
    out = tcrkit.stitch('TRBV19', 'TRBJ2-7', 'ASSIRSSYEQY', chain='beta')
    r = out.iloc[0]
    assert r['error'] != '-'
    assert r['sequence'] == '-'
    # stitch stitches exactly what it was given - no repair columns, no retry bookkeeping
    assert list(out.columns) == ['v', 'j', 'cdr3', 'sequence', 'sequence_aa',
                                 'error', 'stitchr_warnings']
    for gone in ('cdr3_fix', 'cdr3_changed', 'junction_fix', 'species_fix',
                 'n_attempts', 'attempt_log'):
        assert gone not in out.columns, f'{gone} should no longer be produced'


def test_stitch_takes_no_fix_argument():
    """A leftover fix= must fail loudly, not be silently ignored."""
    with pytest.raises(TypeError, match='add_cdr3_junctions'):
        tcrkit.stitch('TRBV19', 'TRBJ2-7', 'ASSIRSSYEQY', chain='beta', fix=True)
    with pytest.raises(TypeError, match='unexpected keyword'):
        tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', typo=1)


def test_fix_junction_then_stitch_is_the_two_step_flow():
    """The documented pipeline: repair the CDR3 column, then stitch it."""
    df = pd.DataFrame({'v': ['TRBV19', 'TRAV1-2'],
                       'j': ['TRBJ2-7', 'TRAJ33'],
                       'cdr3': ['ASSIRSSYEQY', 'AVMDSNYQLI'],
                       'chain': ['beta', 'alpha']})

    # step 1 - the J gene decides the 118 residue (F for TRBJ2-7, W for TRAJ33)
    df['cdr3'] = tcrkit.add_cdr3_junctions(df)['cdr3_fixed']
    assert list(df['cdr3']) == ['CASSIRSSYEQYF', 'CAVMDSNYQLIW']

    # step 2
    out = tcrkit.stitch(df, progress_bar=False)
    assert (out['error'] == '-').all()
    # identical to stitching those completed CDR3s directly
    assert out['sequence'][0] == tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF',
                                               chain='beta')['sequence'][0]


def test_stitch_run_anarci_reads_the_stitched_sequence_back():
    """The stitched sequence must number back to the CDR3 that went into it."""
    out = tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta',
                        run_anarci=True, progress_bar=False)
    r = out.iloc[0]
    assert r['error'] == '-'
    assert r['chain_type_anarci'] == 'B'
    assert r['CDR3_aa_anarci'] in r['cdr3']           # ANARCI reports 105-117
    assert bool(r['cdr3_match_anarci']) is True


def test_constant_auto_uses_stitchrs_j_gene_aware_rule():
    """DEFAULT_CONSTANT is flat; Stitchr picks TRBC2*01 for a TRBJ2 rearrangement."""
    def beta(v, j, cdr3, **kw):
        return tcrkit.stitch(v, j, cdr3, chain='beta', **kw)['sequence_aa'][0]

    # TRBJ2 - the two policies disagree
    mine = beta('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF')
    auto = beta('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', constant='auto')
    assert mine != auto
    assert auto == beta('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', constant='TRBC2*01')
    assert mine == beta('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', constant='TRBC1*01')

    # TRBJ1, and alpha - they agree
    assert (beta('TRBV15', 'TRBJ1-2', 'CATSDLKVHSGNYGYTF', constant='auto')
            == beta('TRBV15', 'TRBJ1-2', 'CATSDLKVHSGNYGYTF'))


def test_constant_spellings():
    def beta(**kw):
        return tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF',
                             chain='beta', **kw)['sequence_aa'][0]

    assert beta(constant='constant') == beta()          # synonym for the default
    assert beta(constant='default') == beta()
    assert beta(constant={'beta': 'TRBC2*01'}) == beta(constant='TRBC2*01')


def test_constant_auto_does_not_undo_no_leader():
    """Stitchr's autofill fills the leader from V too; no_leader must still win."""
    with_leader = tcrkit.stitch('TRAV12-2', 'TRAJ23', 'CAVNTGGGNKLTF', chain='alpha',
                                constant='auto')['sequence_aa'][0]
    without = tcrkit.stitch('TRAV12-2', 'TRAJ23', 'CAVNTGGGNKLTF', chain='alpha',
                            constant='auto', no_leader=True)['sequence_aa'][0]
    assert len(without) < len(with_leader)


def test_stitch_missing_input_is_reported_not_raised():
    r = tcrkit.stitch(None, 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta').iloc[0]
    assert r['sequence'] == '-'
    assert r['error'] == 'missing v/j/cdr3 input'


def test_stitch_frame_without_a_col_mapping_is_a_clear_error():
    """tcrkit must not guess column names - that was the source project's convention."""
    df = pd.DataFrame([{'v_call': 'TRBV19', 'j_call': 'TRBJ2-7',
                        'cdr3_aa': 'CASSIRSSYEQYF'}])
    with pytest.raises(ValueError, match='col_mapping'):
        tcrkit.stitch(df, chain='beta')


def test_stitch_takes_chain_and_species_per_row():
    df = pd.DataFrame([
        {'v': 'TRBV19', 'j': 'TRBJ2-7', 'cdr3': 'CASSIRSSYEQYF', 'chain': 'beta'},
        {'v': 'TRAV12-2', 'j': 'TRAJ23', 'cdr3': 'CAVNTGGGNKLTF', 'chain': 'alpha'},
    ])
    out = tcrkit.stitch(df, progress_bar=False)
    assert (out['error'] == '-').all()
    assert list(out['sequence'].str.len()) == [924, 819]


def test_stitch_reports_a_bad_row_without_aborting_the_batch():
    out = tcrkit.stitch([('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF'),
                         ('NOTAGENE', 'TRBJ2-7', 'CASSIRSSYEQYF')],
                        chain='beta', progress_bar=False)
    assert out['error'][0] == '-'
    assert out['error'][1] != '-'


# ---------------------------------------------------------------------------
# 2. anarci
# ---------------------------------------------------------------------------

def test_anarci_numbers_a_beta_chain():
    r = tcrkit.anarci(BETA_AA).iloc[0]
    assert r['error'] == '-'
    assert r['chain_type'] == 'B'
    assert r['species'] == 'human'
    assert r['v_gene'] == 'TRBV6-2*01'
    assert r['j_gene'] == 'TRBJ2-3*01'
    # CDRs are IMGT 27-38 / 56-65 / 105-117, so CDR3 excludes the C104 and F118
    assert r['CDR1_aa'] == 'MNHDY'
    assert r['CDR2_aa'] == 'SVGEGT'
    assert r['CDR3_aa'] == 'ASSFTDTQY'
    # positions index back into the sequence that was passed in
    assert BETA_AA[r['CDR3_aa_start'] - 1:r['CDR3_aa_end']] == r['CDR3_aa']


def test_anarci_allow_prevents_an_alpha_being_called_delta():
    """Without allow=, ANARCI can match the overlapping delta locus and report a TRDJ."""
    loose = tcrkit.anarci(ALPHA_AA).iloc[0]
    tight = tcrkit.anarci(ALPHA_AA, allow={'A', 'B'}).iloc[0]

    assert loose['chain_type'] == 'D' and loose['j_gene'].startswith('TRDJ')
    assert tight['chain_type'] == 'A' and tight['j_gene'] == 'TRAJ29*01'


def test_anarci_allow_none_does_not_restrict():
    """allow=None must mean 'no restriction', not crash ANARCI."""
    r = tcrkit.anarci(BETA_AA, allow=None).iloc[0]
    assert r['error'] == '-'
    assert r['chain_type'] == 'B'


def test_anarci_works_without_an_hmmscan_binary():
    """The point of the bundled pyhmmer backend: no hmmscan needed."""
    import shutil

    from tcrkit.anarci import backend

    tcrkit.anarci(BETA_AA)                    # force the backend decision
    if shutil.which('hmmscan') is None:
        assert backend() == 'pyhmmer'
    assert backend() in ('pyhmmer', 'hmmscan')


def test_anarci_keeps_the_ids_you_gave_it():
    out = tcrkit.anarci({'b': BETA_AA, 'a': ALPHA_AA}, allow={'A', 'B'},
                        progress_bar=False)
    assert list(out.index) == ['b', 'a']
    assert list(out['chain_type']) == ['B', 'A']


def test_anarci_returns_one_row_per_input_including_failures():
    out = tcrkit.anarci([BETA_AA, 'NOTAREALSEQUENCE'], progress_bar=False)
    assert len(out) == 2
    assert out.iloc[0]['error'] == '-'
    assert out.iloc[1]['error'] != '-'      # reported, not raised
    assert out.iloc[1]['CDR3_aa'] == '-'


def test_anarci_imgt_index_maps_positions_to_residues():
    out = tcrkit.anarci(BETA_AA, allow={'A', 'B'}, include_imgt_index=True,
                        progress_bar=False)
    imgt = out['imgt_dict'].iloc[0]
    assert imgt['104'] == 'C' and imgt['118'] == 'F'
    assert ''.join(imgt[str(p)] for p in range(105, 118) if str(p) in imgt) == \
        out['CDR3_aa'].iloc[0]


def test_anarci_separate_index_aligns_chains_on_shared_positions():
    wide = tcrkit.anarci({'a': ALPHA_AA, 'b': BETA_AA, 'junk': 'NOTASEQUENCE'},
                         allow={'A', 'B'}, separate_index=True, progress_bar=False)
    assert list(wide.index) == ['_position', '_suffix', 'a', 'b', 'junk']
    assert wide.loc['_position', '104'] == 104          # positions are ints to slice on
    assert (wide.loc[['a', 'b'], '104'] == 'C').all()   # same column, both chains
    assert (wide.loc[['a', 'b'], '118'] == 'F').all()
    positions = wide.loc['_position'] != -1
    assert (wide.loc['junk', positions] == '-').all()   # a failure gaps out, not raises
    assert wide.loc['junk', 'error'] != '-'


# ---------------------------------------------------------------------------
# 3. normalize_tcr
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('symbol, species, expected', [
    ('TCRBV20S1', 'human', 'TRBV20-1'),
    ('TRAV23/6', 'human', 'TRAV23/DV6'),        # IMGT TRAV/DV dual naming
    ('TCRAV13S1', 'human', 'TRAV22'),
    ('TRBV20-1*01', 'human', 'TRBV20-1*01'),
    ('TRAV14D-3/DV8', 'mouse', 'TRAV14D-3/DV8'),
    ('TRBV6-5', 'any', 'TRBV6-5'),              # species detected
])
def test_normalize_tcr(symbol, species, expected):
    assert tcrkit.normalize_tcr(symbol, species=species)['standardized'][0] == expected


def test_normalize_tcr_reports_the_species_it_resolved_under():
    out = tcrkit.normalize_tcr(['TRBV6-5', 'TRAV14D-3/DV8'])     # species='any'
    assert list(out['species']) == ['human', 'mouse']


def test_normalize_tcr_unresolvable_is_empty_not_an_error():
    out = tcrkit.normalize_tcr('junk-not-a-gene', species='human', log_failures=False)
    assert pd.isna(out['standardized'][0])


def test_species_may_be_a_scalar_a_sequence_or_a_column():
    """The same precedence rule everywhere: sequence > scalar > input column > default."""
    out = tcrkit.normalize_tcr(['TRBV20-1'] * 3, species=['human', 'mouse', 'any'],
                               log_failures=False)
    assert list(out['species']) == ['human', 'mouse', 'human']
    # an explicit scalar beats a species column the frame happens to carry
    df = pd.DataFrame({'symbol': ['TRBV6-5'], 'species': ['mouse']})
    assert tcrkit.normalize_tcr(df, species='human')['species'][0] == 'human'
    with pytest.raises(ValueError, match='entries but the input has'):
        tcrkit.normalize_tcr(['a', 'b'], species=['human'])


def test_an_unrecognised_species_raises_only_when_you_passed_it():
    """A typo in the keyword is a bug; the same junk in the caller's column is data."""
    with pytest.raises(ValueError, match='Unrecognised species'):
        tcrkit.normalize_tcr('TRBV20-1', species='zebra')
    with pytest.raises(ValueError, match='Unrecognised species'):
        tcrkit.normalize_tcr(['TRBV20-1'] * 2, species=['human', 'zebra'])
    with pytest.raises(ValueError, match='Unrecognised species'):
        tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', species='zebra')

    # but a free-text organism column falls back to 'any' instead of blowing up
    df = pd.DataFrame({'symbol': ['TRBV20-1'], 'species': ['zebra']})
    assert tcrkit.normalize_tcr(df)['standardized'][0] == 'TRBV20-1'


def test_a_valid_but_wrong_species_silently_changes_the_gene():
    """Documented trap: tidytcells coerces rather than refusing, so check this yourself."""
    out = tcrkit.normalize_tcr(['TRBV20-1', 'TRBV20-1'], species=['human', 'mouse'],
                               log_failures=False)
    assert out['standardized'][0] == 'TRBV20-1'
    assert out['standardized'][1] == 'TRBV20'      # a different gene, and no error


def test_stitch_chain_may_be_a_sequence():
    out = tcrkit.stitch(['TRBV19', 'TRAV12-2'], ['TRBJ2-7', 'TRAJ23'],
                        ['CASSIRSSYEQYF', 'CAVNTGGGNKLTF'],
                        chain=['beta', 'alpha'], progress_bar=False)
    assert (out['error'] == '-').all()
    assert list(out['sequence'].str.len()) == [924, 819]


def test_normalize_tcr_per_row_species_from_a_column():
    df = pd.DataFrame({'symbol': ['TRBV6-5', 'TRAV14D-3/DV8'],
                       'species': ['human', 'mouse']})
    out = tcrkit.normalize_tcr(df)
    assert list(out['standardized']) == ['TRBV6-5', 'TRAV14D-3/DV8']


def test_to_gene_drops_the_allele():
    assert tcrkit.to_gene('TRBV20-1*01', species='human')['standardized'][0] == 'TRBV20-1'


def test_genes_match_is_a_plain_predicate():
    assert tcrkit.genes_match('TRAV23/DV6', 'TRAV23/6')
    assert tcrkit.genes_match('TRAV10', 'TRAV10*01')
    assert not tcrkit.genes_match('TRAV10', 'TRBV19')


# ---------------------------------------------------------------------------
# 4. normalize_mhc
# ---------------------------------------------------------------------------

def test_normalize_mhc_class_i():
    r = tcrkit.normalize_mhc('A*0201').iloc[0]
    assert r['allele_input'] == 'A*0201'
    assert r['allele'] == 'A0201'
    assert r['gene'] == 'A'
    assert r['mhc_class'] == 'I'
    assert r['species'] == 'human'


def test_normalize_mhc_echoes_the_input_first():
    out = tcrkit.normalize_mhc(['A*0201', 'H2-Kb'])
    assert out.columns[0] == 'allele_input'
    assert list(out['allele_input']) == ['A*0201', 'H2-Kb']


def test_normalize_mhc_mouse():
    r = tcrkit.normalize_mhc('H2-Kb').iloc[0]
    assert r['species'] == 'mouse'
    assert r['mhc_prefix'] == 'H2'


def test_normalize_mhc_pair_is_ordered_by_gene_not_by_input():
    """'beta__alpha' comes back canonically, with the swap recorded."""
    r = tcrkit.normalize_mhc('DRB5*0101__DRA*0101').iloc[0]
    assert r['mhc_a_allele'] == 'DRA0101'
    assert r['mhc_b_allele'] == 'DRB5*0101'
    assert bool(r['chain_order_swapped']) is True
    assert r['mhc_class'] == 'II'


def test_normalize_mhc_paired_implies_the_paired_layout():
    """paired=True must give a/b columns for every class, not just class II."""
    for allele in ('A*0201', 'H2-Kb', 'DRB1*04:01'):
        r = tcrkit.normalize_mhc(allele, paired=True, squeeze=True)
        assert 'mhc_a_allele' in r and 'mhc_b_allele' in r, allele
        assert 'allele' not in r, f'{allele} fell back to the flat layout'
    # a class-I heavy chain is paired with its invariant B2M light chain
    assert tcrkit.normalize_mhc('A*0201', paired=True,
                                squeeze=True)['mhc_b_allele'] == 'B2M_human'
    assert tcrkit.normalize_mhc('H2-Kb', paired=True,
                                squeeze=True)['mhc_b_allele'] == 'B2M_mouse'


def test_mhc_chain_of():
    assert tcrkit.mhc_chain_of('HLA-A*02:01') == 'alpha'   # class-I heavy
    assert tcrkit.mhc_chain_of('HLA-DRB1*04:01') == 'beta'
    assert tcrkit.mhc_chain_of('B2M') == 'beta'            # class-I light


def test_normalize_mhc_coerce_records_failures():
    out = tcrkit.normalize_mhc(['A*0201', 'not-an-allele'], errors='coerce')
    assert len(out) == 2
    assert out.iloc[0]['allele'] == 'A0201'
    assert pd.isna(out.iloc[1]['allele'])
    assert isinstance(out.iloc[1]['error'], str)


def test_normalize_mhc_raise_is_the_default():
    with pytest.raises(Exception):
        tcrkit.normalize_mhc(['A*0201', 'not-an-allele'])


# ---------------------------------------------------------------------------
# mhc_sequence
# ---------------------------------------------------------------------------
# The MHC sequence reference is derived from IMGT/HLA and is deliberately not in version
# control, so a fresh clone cannot run these. Everything else in the suite can.

needs_mhc_ref = pytest.mark.skipif(
    not bundled_mhc_table_path().exists(),
    reason='mhc_sequences.csv is not in version control - see tcrkit.pmhc.load_mhc_table')


def test_missing_mhc_reference_raises_with_instructions(tmp_path):
    with pytest.raises(FileNotFoundError, match='not distributed with the source'):
        load_mhc_table(tmp_path / 'absent.csv')


@needs_mhc_ref
def test_mhc_sequence_lookup():
    r = tcrkit.mhc_sequence('HLA-A*02:01').iloc[0]
    assert r['allele'] == 'A0201'                 # normalized to the reference's key
    assert r['sequence'].startswith('GSHSMRYFFTSVSRPGRGEP')
    assert r['length'] == len(r['sequence'])
    assert r['mhc_class'] == 'I' and r['species'] == 'human'
    assert pd.isna(r['error'])


@needs_mhc_ref
def test_mhc_sequence_accepts_any_spelling_of_the_same_allele():
    out = tcrkit.mhc_sequence(['HLA-A*02:01', 'A*0201', 'A0201'])
    assert out['allele'].nunique() == 1
    assert out['sequence'].nunique() == 1


@needs_mhc_ref
def test_mhc_sequence_miss_is_reported_not_raised():
    out = tcrkit.mhc_sequence(['A*0201', 'A2']).iloc[1]     # A2 is a serotype
    assert pd.isna(out['sequence'])
    assert 'not in the reference' in out['error']


@needs_mhc_ref
def test_mhc_sequence_paired_resolves_both_chains():
    out = tcrkit.mhc_sequence(['A*0201', 'DRB1*04:01'], paired=True)
    cls_i, cls_ii = out.iloc[0], out.iloc[1]
    # class I: heavy chain + the invariant B2M light chain, both found
    assert cls_i['mhc_a_allele'] == 'A0201' and cls_i['mhc_b_allele'] == 'B2M_human'
    assert isinstance(cls_i['mhc_a_sequence'], str)
    assert isinstance(cls_i['mhc_b_sequence'], str)
    # class II: the alpha/beta pair
    assert cls_ii['mhc_a_allele'] == 'DRA0101'
    assert cls_ii['mhc_b_allele'] == 'DRB1*0401'
    assert isinstance(cls_ii['mhc_b_sequence'], str)


@needs_mhc_ref
def test_mhc_sequence_finds_a_species_tagged_b2m_on_its_own():
    """Regression: the class-I light chain looked up as a single allele.

    normalize_mhc emits 'B2M_human' / 'B2M_mouse' for a pair, and the reference is keyed
    that way - but normalizing one of those labels on its own goes through mhcgnomes and
    comes back as a bare 'B2M' / 'H2-B2M' with the species dropped, which used to miss.
    Real pMHC tables carry mhc.b as its own column, so this is the common case.
    """
    out = tcrkit.mhc_sequence(['B2M_human', 'B2M_mouse']).set_index('allele_input')
    for label in ('B2M_human', 'B2M_mouse'):
        assert isinstance(out.loc[label, 'sequence'], str), label
        assert pd.isna(out.loc[label, 'error']), label
    # the two species have different B2M sequences
    assert out.loc['B2M_human', 'sequence'] != out.loc['B2M_mouse', 'sequence']


@needs_mhc_ref
def test_mhc_sequence_full_is_longer_than_the_domain():
    short = tcrkit.mhc_sequence('A*0201')['sequence'][0]
    full = tcrkit.mhc_sequence('A*0201', full=True)['sequence'][0]
    assert full.startswith(short) and len(full) > len(short)


@needs_mhc_ref
def test_mhc_reference_has_the_expected_shape():
    tab = load_mhc_table()
    assert len(tab) == 251
    assert {'allele', 'sequence', 'full_sequence', 'gene', 'mhc_class',
            'species'} <= set(tab.columns)
    assert tab['allele'].is_unique


# ---------------------------------------------------------------------------
# 5. add_cdr3_junctions
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('cdr3, j, expected', [
    ('ASSYSGNTEAF', 'TRBJ1-1', 'CASSYSGNTEAFF'),      # gains both flanks
    ('AVMDSNYQLI', 'TRAJ33*01', 'CAVMDSNYQLIW'),      # allele suffix tolerated
    ('CASSYSGNTEAFF', 'TRBJ1-1', 'CASSYSGNTEAFF'),    # already flanked -> untouched
    ('ASSIRSSYEQY', 'TRBJ2-7', 'CASSIRSSYEQYF'),
    ('AVRD', 'TRAJ35', 'CAVRDC'),                     # TRAJ35 uses C/C, not C/F
])
def test_add_cdr3_junctions(cdr3, j, expected):
    out = tcrkit.add_cdr3_junctions(cdr3, j).iloc[0]
    assert out['cdr3_fixed'] == expected
    assert out['changed'] == (expected != cdr3)


def test_add_cdr3_junctions_unknown_gene_is_empty():
    out = tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'NOTAGENE').iloc[0]
    assert pd.isna(out['cdr3_fixed'])
    assert bool(out['changed']) is False


def test_add_cdr3_junctions_keep_unknown_falls_back_to_the_original():
    out = tcrkit.add_cdr3_junctions(['ASSYSGNTEAF', 'ASSYSGNTEAF'],
                              ['TRBJ1-1', 'NOTAGENE'], keep_unknown=True)
    assert list(out['cdr3_fixed']) == ['CASSYSGNTEAFF', 'ASSYSGNTEAF']


def test_optional_fields_are_not_echoed_when_unused():
    out = tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1')
    assert 'v' not in out.columns and 'species' not in out.columns
    assert 'v' in tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1', 'TRBV19').columns


def test_add_cdr3_junctions_reports_the_lookup_it_used():
    """The table's reasoning is visible, not hidden."""
    out = tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1').iloc[0]
    assert out['junction'] == 'C/F'
    assert out['purity'] == pytest.approx(1.0)
    assert out['n'] == 211265
    assert out['level'] == 'J'


def test_junction_of():
    out = tcrkit.junction_of(['TRBJ1-1', 'TRAJ33', 'TRAJ35']).set_index('j')
    assert list(out['junction']) == ['C/F', 'C/W', 'C/C']
    assert out.loc['TRAJ35', 'purity'] == pytest.approx(0.999932, abs=1e-6)
    assert set(out['level']) == {'J'}


def test_junction_of_unknown_gene():
    out = tcrkit.junction_of('NOTAGENE').iloc[0]
    assert out['junction'] is None or pd.isna(out['junction'])
    assert out['n'] == 0


def test_bundled_lookup_table_is_installed():
    """The lookup table must ship with the package, not be read from a project path."""
    from tcrkit.cdr3_junction import bundled_junction_table_path, load_junction_table

    assert bundled_junction_table_path().exists()
    tab = load_junction_table()
    assert len(tab) > 1000
    assert {'level', 'species', 'v_gene', 'j_gene', 'junction',
            'purity', 'n', 'locus', 'alphabet'} <= set(tab.columns)


# ---------------------------------------------------------------------------
# packaging
# ---------------------------------------------------------------------------

def test_stitch_pair_stitches_both_chains_in_one_call():
    """The one function that is not per-chain: a whole receptor from one row."""
    df = pd.DataFrame([{'a_v': 'TRAV12-2', 'a_j': 'TRAJ23', 'a_cdr3': 'CAVNTGGGNKLTF',
                        'b_v': 'TRBV19', 'b_j': 'TRBJ2-7', 'b_cdr3': 'CASSIRSSYEQYF'}])
    out = tcrkit.stitch_pair(
        df, col_mapping={'alpha_v': 'a_v', 'alpha_j': 'a_j', 'alpha_cdr3': 'a_cdr3',
                         'beta_v': 'b_v', 'beta_j': 'b_j', 'beta_cdr3': 'b_cdr3'},
        species='human')
    assert out['error'][0] == '-'
    for chain in ('alpha', 'beta'):
        assert len(out[f'{chain}_sequence_aa'][0]) > 100


def test_read_sequences_reads_plain_and_fasta(tmp_path):
    plain = tmp_path / 'seqs.txt'
    plain.write_text(f'{BETA_AA}\n{ALPHA_AA}\n')
    assert tcrkit.read_sequences(plain) == [('id_1', BETA_AA), ('id_2', ALPHA_AA)]

    fasta = tmp_path / 'seqs.fa'
    fasta.write_text(f'>beta\n{BETA_AA}\n>alpha\n{ALPHA_AA}\n')
    assert tcrkit.read_sequences(fasta) == [('beta', BETA_AA), ('alpha', ALPHA_AA)]
    # and the ids carry through into anarci
    assert list(tcrkit.anarci(tcrkit.read_sequences(fasta),
                              allow={'A', 'B'}, progress_bar=False).index) == ['beta',
                                                                                'alpha']

    with pytest.raises(FileNotFoundError):
        tcrkit.read_sequences(tmp_path / 'nope.txt')


def test_cdr_imgt_locs_matches_what_anarci_reports():
    """The documented IMGT spans must be the ones the CDR columns actually use."""
    assert tcrkit.CDR_IMGT_LOCS == {'cdr1': (27, 38), 'cdr2': (56, 65),
                                    'cdr3': (105, 117)}
    r = tcrkit.anarci(BETA_AA, allow={'A', 'B'}).iloc[0]
    # CDR3 is 105-117, so it excludes the conserved C104 / F118
    assert BETA_AA[r['CDR3_aa_start'] - 1:r['CDR3_aa_end']] == r['CDR3_aa']
    assert not r['CDR3_aa'].startswith('C')


def test_build_junction_table_learns_from_your_own_reference():
    ref = pd.DataFrame({'cdr3': ['CASSYSGNTEAFF'] * 5 + ['CAVMDSNYQLIW'] * 5,
                        'j': ['TRBJ1-1'] * 5 + ['TRAJ33'] * 5,
                        'sp': ['human'] * 10})
    tab = tcrkit.build_junction_table(ref, 'cdr3', 'j', species_col='sp')
    assert len(tab) > 0
    assert {'level', 'species', 'j_gene', 'junction'} <= set(tab.columns)
    # a table built from only these two J genes still completes them correctly
    assert tcrkit.add_cdr3_junctions('ASSYSGNTEAF', 'TRBJ1-1', table=tab,
                                     squeeze=True) == 'CASSYSGNTEAFF'
    # and knows nothing about a J gene it never saw
    assert pd.isna(tcrkit.add_cdr3_junctions('ASSIRSSYEQY', 'TRBJ2-7', table=tab,
                                             squeeze=True))


def test_cli_entry_point_builds():
    from tcrkit.cli import build_parser

    parser = build_parser()
    for cmd in ('stitch', 'anarci', 'normalize-tcr', 'normalize-mhc',
                'add-cdr3-junctions'):
        assert parser.parse_args([cmd]).command == cmd


@pytest.mark.parametrize('argv', [
    ['stitch', '--v', 'TRBV19', '--j', 'TRBJ2-7', '--cdr3', 'CASSIRSSYEQYF'],
    ['anarci', '--seq', BETA_AA, '--allow', 'A', 'B', '-q'],
    ['normalize-tcr', 'TCRBV20S1', '--species', 'human'],
    ['normalize-mhc', 'A*0201'],
    ['add-cdr3-junctions', '--cdr3', 'ASSYSGNTEAF', '--j', 'TRBJ1-1'],
])
def test_cli_subcommands_actually_run(argv, capsys):
    """Dispatch each subcommand, not just parse it - imports only fail at dispatch."""
    from tcrkit.cli import main

    assert main(argv) == 0
    assert capsys.readouterr().out.strip()


# ---------------------------------------------------------------------------
# build_tcr_pmhc - all five capabilities end to end
# ---------------------------------------------------------------------------

PMHC_RAW = pd.DataFrame([
    {'epitope': 'GILGFVFTL',
     'alpha_v': 'TRAV1-2', 'alpha_j': 'TRAJ33', 'alpha_cdr3': 'AVMDSNYQLI',
     'beta_v': 'TCRBV20S1', 'beta_j': 'TRBJ2-7', 'beta_cdr3': 'ASSLGQAYEQY',
     'mhc_a': 'HLA-A*02:01'},
    {'epitope': 'BROKEN',
     'alpha_v': 'TRAV1-2', 'alpha_j': 'TRAJ33', 'alpha_cdr3': 'AVMDSNYQLI',
     'beta_v': 'NOTAGENE', 'beta_j': 'TRBJ2-7', 'beta_cdr3': 'ASSLGQAYEQY',
     'mhc_a': 'A2'},
])


@needs_mhc_ref
def test_build_tcr_pmhc_runs_every_step():
    out = tcrkit.build_tcr_pmhc(PMHC_RAW, progress_bar=False)
    good = out.iloc[0]

    assert good['epitope'] == 'GILGFVFTL'            # extra columns carried through
    assert good['beta_v'] == 'TRBV20-1'              # normalized
    assert good['beta_cdr3'] == 'CASSLGQAYEQYF'      # junction completed
    assert good['alpha_cdr3'] == 'CAVMDSNYQLIW'      # ... using the J gene, so W not F
    assert len(good['beta_sequence_aa']) > 250       # stitched
    assert good['mhc_a_allele'] == 'A0201'           # MHC resolved, both chains
    assert good['mhc_b_allele'] == 'B2M_human'
    assert isinstance(good['mhc_a_sequence'], str)
    assert good['beta_anarci_v'] == 'TRBV20-1*01'    # ANARCI agrees with the V call
    assert bool(good['beta_cdr3_ok']) and bool(good['alpha_cdr3_ok'])
    assert good['error'] == '-' and bool(good['ok'])


@needs_mhc_ref
def test_build_tcr_pmhc_reports_per_row_instead_of_raising():
    out = tcrkit.build_tcr_pmhc(PMHC_RAW, run_anarci=False, progress_bar=False)
    bad = out.iloc[1]
    assert not bool(bad['ok'])
    assert 'beta V call not recognised' in bad['error']
    assert 'did not stitch' in bad['error']
    assert "'A2' not in the reference" in bad['error']
    # the good row is still built
    assert bool(out.iloc[0]['ok'])


@needs_mhc_ref
def test_build_tcr_pmhc_single_chain_and_col_mapping():
    df = pd.DataFrame([{'v_call': 'TRBV19', 'j_call': 'TRBJ2-7',
                        'junction_aa': 'ASSIRSSYEQY', 'mhc.a': 'DRA*01:01',
                        'mhc.b': 'DRB1*04:01', 'study': 's1'}])
    out = tcrkit.build_tcr_pmhc(
        df, col_mapping={'beta_v': 'v_call', 'beta_j': 'j_call',
                         'beta_cdr3': 'junction_aa', 'mhc_a': 'mhc.a',
                         'mhc_b': 'mhc.b'}, progress_bar=False)
    assert out['study'][0] == 's1'
    assert out['beta_cdr3'][0] == 'CASSIRSSYEQYF'
    assert 'alpha_sequence' not in out.columns        # no alpha columns -> skipped
    assert out['mhc_a_allele'][0] == 'DRA0101'        # explicit class-II pair
    assert out['mhc_b_allele'][0] == 'DRB1*0401'
    assert out['mhc_class'][0] == 'II'
    assert bool(out['ok'][0])


MY_MHC_TABLE = pd.DataFrame({
    'allele': ['HLA-A*02:01', 'B2M_human'],              # deliberately not normalized
    'sequence': ['GSHSMRYFFTSVSRPGRGEP', 'IQRTPKIQVYSRHPAENGK'],
})


def test_mhc_sequence_takes_your_own_table():
    out = tcrkit.mhc_sequence(['A*0201', 'DRB1*04:01'], table=MY_MHC_TABLE)
    assert out['sequence'][0] == 'GSHSMRYFFTSVSRPGRGEP'   # found despite the spelling
    assert pd.isna(out['sequence'][1])                    # absent from my table
    assert 'not in the reference' in out['error'][1]


def test_mhc_table_needs_allele_and_sequence_columns():
    with pytest.raises(ValueError, match="needs \\['sequence'\\]"):
        tcrkit.mhc_sequence('A*0201', table=pd.DataFrame({'allele': ['A0201']}))


def test_build_tcr_pmhc_survives_a_missing_mhc_reference(monkeypatch):
    def gone(*args, **kwargs):
        raise FileNotFoundError('no reference here.')

    monkeypatch.setattr('tcrkit.pmhc.load_mhc_table', gone)
    with pytest.warns(UserWarning, match='MHC sequences left empty'):
        out = tcrkit.build_tcr_pmhc(PMHC_RAW.iloc[[0]], run_anarci=False,
                                    progress_bar=False)
    assert out['ok'].all()                                # the TCR half still built
    assert out['mhc_a_allele'][0] == 'A0201'              # allele still normalized
    assert out['mhc_a_sequence'].isna().all()             # just no sequences


def test_build_tcr_pmhc_takes_your_own_mhc_table():
    out = tcrkit.build_tcr_pmhc(PMHC_RAW.iloc[[0]], run_anarci=False, progress_bar=False,
                                mhc_table=MY_MHC_TABLE)
    assert out['mhc_a_sequence'][0] == 'GSHSMRYFFTSVSRPGRGEP'
    assert out['mhc_b_allele'][0] == 'B2M_human'          # partner inferred for class I


def test_build_tcr_pmhc_needs_at_least_one_chain():
    with pytest.raises(ValueError, match='No chain found'):
        tcrkit.build_tcr_pmhc(pd.DataFrame([{'epitope': 'GILGFVFTL'}]))
