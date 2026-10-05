"""Smoke tests for the repertoire metrics.

Every metric except `morisita_horn` wraps an external package, so those tests are
skipped when the ``metrics`` extra is not installed and the suite still passes in a
plain install. `morisita_horn` is computed in tcrkit and always runs.
"""

import importlib.util
import math
import warnings

import pandas as pd
import pytest

import tcrkit

has_skbio = pytest.mark.skipif(importlib.util.find_spec('skbio') is None,
                               reason="needs the 'metrics' extra (scikit-bio)")
has_olga = pytest.mark.skipif(importlib.util.find_spec('olga') is None,
                              reason="needs the 'metrics' extra (olga)")
has_tcrdist = pytest.mark.skipif(importlib.util.find_spec('tcrdist') is None,
                                 reason="needs the 'metrics' extra (tcrdist3)")


# ---------------------------------------------------------------------------
# diversity
# ---------------------------------------------------------------------------

@has_skbio
def test_diversity_counts_raw_clone_labels():
    out = tcrkit.diversity(['A', 'A', 'A', 'B', 'C']).iloc[0]
    assert out['n_clones'] == 3
    assert out['n_cells'] == 5
    assert out['top_clone_fraction'] == pytest.approx(0.6)


@has_skbio
def test_diversity_accepts_abundances_directly():
    labels = tcrkit.diversity(['A', 'A', 'A', 'B', 'C']).iloc[0]
    counts = tcrkit.diversity([3, 1, 1]).iloc[0]
    for metric in ('n_clones', 'n_cells', 'shannon', 'clonality'):
        assert labels[metric] == pytest.approx(counts[metric])


@has_skbio
def test_entropy_is_in_bits_and_maximal_when_even():
    """Four equally abundant clones carry exactly log2(4) = 2 bits."""
    even = tcrkit.diversity([25, 25, 25, 25]).iloc[0]
    assert even['shannon'] == pytest.approx(2.0)
    assert even['pielou_evenness'] == pytest.approx(1.0)
    assert even['clonality'] == pytest.approx(0.0)


@has_skbio
def test_clonality_is_one_for_a_single_clone():
    out = tcrkit.diversity([100]).iloc[0]
    assert out['n_clones'] == 1
    assert out['clonality'] == pytest.approx(1.0)
    assert out['top_clone_fraction'] == pytest.approx(1.0)


@has_skbio
def test_clonality_orders_repertoires_as_expected():
    """A repertoire dominated by one clone must score higher than an even one."""
    df = pd.DataFrame({'sample': ['dominated'] * 3 + ['even'] * 3,
                       'cdr3': list('ABC') * 2,
                       'count': [100, 1, 1, 34, 33, 33]})
    out = tcrkit.diversity(df, by='sample',
                           col_mapping={'clone': 'cdr3'}).set_index('repertoire')
    assert out.loc['dominated', 'clonality'] > out.loc['even', 'clonality']
    assert out.loc['even', 'clonality'] == pytest.approx(0.0, abs=1e-3)
    assert list(out.index) == ['dominated', 'even']


@has_skbio
def test_diversity_simpson_and_dominance_are_complements():
    """scikit-bio's naming: simpson = 1 - sum(p^2), dominance = sum(p^2)."""
    out = tcrkit.diversity([50, 30, 20]).iloc[0]
    assert out['simpson'] + out['dominance'] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# morisita_horn
# ---------------------------------------------------------------------------

def test_morisita_horn_bounds():
    mh = lambda a, b: tcrkit.morisita_horn(a, b, squeeze=True)
    assert mh({'A': 10, 'B': 5}, {'A': 10, 'B': 5}) == pytest.approx(1.0)
    assert mh({'A': 10, 'B': 5}, {'C': 10, 'D': 5}) == 0.0


def test_morisita_horn_always_returns_a_frame():
    """The uniform contract: every capability hands back a DataFrame."""
    two = tcrkit.morisita_horn({'A': 10}, {'A': 10})
    assert isinstance(two, pd.DataFrame) and two.shape == (2, 2)
    assert isinstance(tcrkit.morisita_horn({'A': 10}, {'A': 10}, long=True),
                      pd.DataFrame)


def test_morisita_horn_ignores_sequencing_depth():
    """The point of the index: the same composition at 10x depth still scores 1."""
    shallow = {'A': 10, 'B': 5, 'C': 1}
    deep = {'A': 100, 'B': 50, 'C': 10}
    assert tcrkit.morisita_horn(shallow, deep, squeeze=True) == pytest.approx(1.0)


def test_morisita_horn_pairwise_matrix():
    df = pd.DataFrame({'sample': ['s1'] * 3 + ['s2'] * 3,
                       'cdr3': list('ABC') * 2,
                       'count': [100, 1, 1, 34, 33, 33]})
    m = tcrkit.morisita_horn(df, by='sample', col_mapping={'clone': 'cdr3'})
    assert list(m.index) == ['s1', 's2']
    assert m.loc['s1', 's1'] == pytest.approx(1.0)
    assert m.loc['s1', 's2'] == pytest.approx(m.loc['s2', 's1'])   # symmetric
    assert 0 < m.loc['s1', 's2'] < 1


def test_morisita_horn_needs_a_second_repertoire():
    with pytest.raises(ValueError, match='other='):
        tcrkit.morisita_horn({'A': 1, 'B': 2})


# ---------------------------------------------------------------------------
# generation_probability (OLGA)
# ---------------------------------------------------------------------------

@has_olga
def test_generation_probability():
    out = tcrkit.generation_probability('CASSIRSSYEQYF').iloc[0]
    assert out['pgen'] == pytest.approx(1.804633841190399e-07, rel=1e-6)
    assert out['log10_pgen'] == pytest.approx(math.log10(out['pgen']))
    assert pd.isna(out['error'])


@has_olga
def test_generation_probability_conditioning_on_genes_lowers_it():
    """Conditioning restricts the recombinations that could have made it.

    V and J act independently in OLGA, so a lone V must count too - it used to be
    silently discarded unless both were given.
    """
    pgen = lambda **kw: tcrkit.generation_probability('CASSIRSSYEQYF', **kw)['pgen'][0]
    loose, v_only, both = pgen(), pgen(v='TRBV19'), pgen(v='TRBV19', j='TRBJ2-7')
    assert 0 < both < loose
    assert v_only < loose                       # a V on its own still narrows it
    assert v_only == pytest.approx(both, rel=1e-3)


@has_olga
def test_genes_are_only_echoed_when_given():
    """An input field the caller left out should not come back as a column of None."""
    assert 'v' not in tcrkit.generation_probability('CASSIRSSYEQYF').columns
    assert 'j' not in tcrkit.generation_probability('CASSIRSSYEQYF').columns
    with_j = tcrkit.generation_probability('CASSIRSSYEQYF', j='TRBJ2-7')
    assert 'j' in with_j.columns and 'v' not in with_j.columns


@has_olga
def test_generation_probability_alpha_chain():
    out = tcrkit.generation_probability('CAVMDSNYQLIW', chain='alpha').iloc[0]
    assert out['pgen'] > 0


@has_olga
def test_generation_probability_bad_cdr3_is_reported_not_printed(capsys):
    out = tcrkit.generation_probability(['CASSIRSSYEQYF', 'NOTACDR3']).iloc[1]
    assert out['pgen'] == 0.0
    assert out['error'] == 'not an amino-acid sequence'
    assert capsys.readouterr().out == ''        # OLGA's own complaint stays muted


@has_olga
def test_generation_probability_unknown_model():
    with pytest.raises(ValueError, match='No OLGA model'):
        tcrkit.generation_probability('CASSIRSSYEQYF', chain='delta')


# ---------------------------------------------------------------------------
# tcrdist (tcrdist3)
# ---------------------------------------------------------------------------

TCRS = pd.DataFrame({
    'beta_cdr3': ['CASSIRSSYEQYF', 'CASSLGQAYEQYF', 'CASSYSIRGSRGEQYF'],
    'beta_v': ['TRBV19', 'TRBV20-1', 'TRBV6-5'],
    'beta_j': ['TRBJ2-7', 'TRBJ2-7', 'TRBJ2-7'],
}, index=['t1', 't2', 't3'])


@has_tcrdist
def test_tcrdist_matrix():
    m = tcrkit.tcrdist(TCRS)
    assert m.shape == (3, 3)
    assert list(m.index) == ['t1', 't2', 't3']     # the caller's index is kept
    assert (m.values.diagonal() == 0).all()        # a TCR is zero distance from itself
    assert (m.values == m.values.T).all()          # symmetric
    assert m.loc['t1', 't2'] > 0


@has_tcrdist
def test_tcrdist_emits_no_warnings():
    """tcrdist3 warns about its own stale db_file check on every stock call."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        tcrkit.tcrdist(TCRS)
    noise = [str(w.message) for w in caught
             if 'db_file must be' in str(w.message)
             or 'IProgress not found' in str(w.message)]
    assert not noise, noise


@has_tcrdist
def test_tcrdist_still_lets_real_warnings_through():
    """Only the two known-noise messages are filtered, not warnings in general."""
    from tcrkit.metrics import _quiet_tcrdist

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        with _quiet_tcrdist():
            warnings.warn('something you should actually see', UserWarning)
    assert [str(w.message) for w in caught] == ['something you should actually see']


@has_tcrdist
def test_tcrdist_long_form():
    long = tcrkit.tcrdist(TCRS, long=True)
    assert list(long.columns) == ['i', 'j', 'distance']
    assert len(long) == 6                          # 3x3 minus the diagonal
    assert (long.i != long.j).all()


@has_tcrdist
def test_tcrdist_missing_column_says_which():
    with pytest.raises(ValueError, match='beta_cdr3'):
        tcrkit.tcrdist(TCRS.drop(columns=['beta_cdr3']))
