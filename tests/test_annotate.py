"""tcrkit.annotate: what a sequence is, without saying first."""

import pandas as pd
import pytest

import tcrkit

BETA = ('MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGTTAKGEVPDGYNVSRLKKQNFLLGLES'
        'AAPSQTSVYFCASSFTDTQYFGPGTRLTVL')
ALPHA = ('AQKITQTQPGMFVQEKEAVTLDCTYDTSDPSYGLFWYKQPSSGEMIFLIYQGSYDQQNATEGRYSLNFQKARKSANLVIS'
         'ASQLGDSAMYFCAMSQLNSGNTPLVFGKGTRLSVIA')


@pytest.fixture(scope='module')
def mixed():
    t = tcrkit.load_mhc_table().set_index('allele')
    seqs = {k: t.at[k, 'sequence'] for k in ('A0201', 'DRA0101', 'DRB1*0101', 'B2M_human',
                                             'H2-Kb')}
    return {**seqs, 'beta': BETA, 'alpha': ALPHA, 'pep': 'GILGFVFTL',
            'junk': 'ACDEFGHIKLMNPQRSTVWY' * 4, 'empty': ''}


def test_every_chain_type_is_recognised(mixed):
    out = tcrkit.annotate(mixed, allow={'A', 'B'}, progress_bar=False)
    assert out['chain_type'].to_dict() == {
        'A0201': 'MHC_I', 'DRA0101': 'MHC_II_alpha', 'DRB1*0101': 'MHC_II_beta',
        'B2M_human': 'B2M', 'H2-Kb': 'MHC_I', 'beta': 'TCR_beta', 'alpha': 'TCR_alpha',
        'pep': 'peptide', 'junk': 'unknown', 'empty': 'unknown'}
    assert out['error'].isna().sum() == 8
    assert out.loc['junk', 'error'] and out.loc['empty', 'error'] == 'empty sequence'


def test_each_row_carries_the_annotation_of_its_kind(mixed):
    out = tcrkit.annotate(mixed, allow={'A', 'B'}, progress_bar=False)
    b, a2 = out.loc['beta'], out.loc['A0201']
    assert b['v_gene'].startswith('TRBV') and b['CDR3_aa'] == 'ASSFTDTQY'
    assert pd.isna(b['gene']) and pd.isna(b['helix'])
    assert (a2['gene'], a2['allele'], a2['mhc_class'], a2['species']) == ('A', 'A0201', 'I', 'human')
    assert a2['helix'] and pd.isna(a2['v_gene'])
    assert out.loc['H2-Kb', 'species'] == 'mouse'


def test_matches_the_specific_functions():
    tcr = tcrkit.anarci(BETA, allow={'A', 'B'}, progress_bar=False).iloc[0]
    got = tcrkit.annotate(BETA, allow={'A', 'B'}, progress_bar=False).iloc[0]
    assert (got['v_gene'], got['j_gene'], got['CDR3_aa']) == \
        (tcr['v_gene'], tcr['j_gene'], tcr['CDR3_aa'])


def test_per_residue_for_any_chain(mixed):
    res = tcrkit.annotate(mixed, per_residue=True, allow={'A', 'B'}, progress_bar=False)
    first = res.groupby(level=0, sort=False)['chain_type'].first()
    assert first['beta'] == 'TCR_beta' and first['A0201'] == 'MHC_I'
    assert first['pep'] == 'peptide' and 'empty' not in first
    for name in ('beta', 'A0201', 'pep'):
        assert len(res.loc[[name]]) == len(mixed[name])
    assert {'position', 'aa', 'chain_type', 'domain', 'imgt', 'region'} <= set(res.columns)


def test_squeeze_and_col_mapping():
    r = tcrkit.annotate('GILGFVFTL', squeeze=True, progress_bar=False)
    assert r['chain_type'] == 'peptide'
    df = pd.DataFrame({'seq': [BETA, 'GILGFVFTL']})
    out = tcrkit.annotate(df, col_mapping={'sequence': 'seq'}, progress_bar=False)
    assert out['chain_type'].tolist() == ['TCR_beta', 'peptide']


def test_cli_annotate(capsys):
    from tcrkit.cli import main

    assert main(['annotate', '--seq', BETA, 'GILGFVFTL', '-q']) == 0
    out = capsys.readouterr().out
    assert 'TCR_beta' in out and 'peptide' in out
