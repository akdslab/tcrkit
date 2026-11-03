"""An MHC chain's gene, IMGT numbering and groove helices (tcrkit.annotate_mhc),
and the gene matching behind it (tcrkit.mhc.best_mhc, mhc_span, mhc_chain)."""

import pandas as pd
import pytest

import tcrkit
from tcrkit import mhc
from tcrkit.cli import main

# HLA-A*02:01 alpha1 + alpha2, as in PDB 1AO7 chain A
HLA_A2 = ('GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLG'
          'TLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLR'
          'AYLEGTCVEWLRRYLENGKETLQRT')
# PDB 4MS8 chain A (H2-Ld): residues missing where alpha1 meets alpha2
H2_LD_GAPPED = ('MGPHSMRYYETATSPRYTSVGYVDDKEFVRFDSDAENPRYEPQVPWMEQEGPEYWERITQIAKGQEQWFRVNLRTLLGY'
                'HTLQWMYGCDVGSDGRLLRGYEQFAYDGCDYIALNEDLRTWTAADMAAQITRRKWEQAGAAEYYRAYLEGECVEWLHR'
                'YLKNG')
# PDB 7RDV chain B (I-E beta): diverges from every IMGT entry at the end of the helix
IE_BETA = ('ERHFVVQFKGECYYTNGTQRIRLVTRYIYNREEYVRYDSDVGEYRAVTELGRPDAEYWNSQPEILERTRAEVDTACRHNY'
           'EGPETSTSLRRLEQPNVAISLSRTEALNHHNTLVCSVTDFYPAKIKVRWFRNGQEETVGVSSTQLIRNGDWTFQVLVML'
           'EMTPHQGEVYTCHVEHPSLKSPITVEW')
DRA = ('IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEIMTKRSN'
       'YTPITN')
B2M = ('IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYA'
       'CRVNHVTLSQPKIVKWDRDM')
BETA = ('MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGTTAKGEVPDGYNVSR'
        'LKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL')


@pytest.fixture(scope='module')
def wide():
    return tcrkit.annotate_mhc({'A2': HLA_A2, '4MS8': H2_LD_GAPPED, 'IE': IE_BETA,
                                'DRA': DRA, 'b2m': B2M, 'tcr': BETA, 'pep': 'LLFGYPVYV'})


def test_gene_domains():
    assert mhc.gene_domains('A') == ('G-ALPHA1', 'G-ALPHA2')
    assert mhc.gene_domains('DQA1') == ('G-ALPHA',)
    assert mhc.gene_domains('DRB1') == ('G-BETA',)
    assert mhc.gene_domains('B2M') == ()
    assert mhc.gene_domains(None) == mhc.gene_domains('H2-K') == tcrkit.gdomain.GDOMAINS
    # mouse: I-A and I-E alpha / beta, H2-K heavy chain
    assert mhc.gene_domains('AA') == mhc.gene_domains('EA') == ('G-ALPHA',)
    assert mhc.gene_domains('AB') == mhc.gene_domains('EB') == ('G-BETA',)
    assert mhc.gene_domains('K') == ('G-ALPHA1', 'G-ALPHA2')


def test_one_row_per_sequence(wide):
    assert list(wide.index) == ['A2', '4MS8', 'IE', 'DRA', 'b2m', 'tcr', 'pep']
    assert list(wide.columns) == ['sequence', 'gene', 'allele', 'allele_identity', 'mhc_class',
                                  'score', 'span',
                                  'domains', 'domain_spans', 'helix', 'regions',
                                  'signal_peptide', 'extracellular', 'transmembrane',
                                  'cytoplasmic', 'error']
    assert wide.loc['A2', ['gene', 'allele', 'mhc_class', 'span', 'domains', 'domain_spans',
                           'helix']].tolist() == \
        ['A', 'A0201', 'I', f'1-{len(HLA_A2)}', 'G-ALPHA1, G-ALPHA2', '1-90, 91-182',
         '50-86, 140-176']
    assert wide.loc['4MS8', 'domains'] == 'G-ALPHA1, G-ALPHA2'
    assert wide.loc['IE', ['mhc_class', 'domains', 'helix']].tolist() == \
        ['II', 'G-BETA', '51-61, 64-88']
    assert wide.loc['DRA', ['gene', 'domains']].tolist() == ['DRA', 'G-ALPHA']


def test_allele_identity_says_how_close_the_allele_is():
    import tcrkit

    # point mutations, as for an allele that is not bundled
    mutant = HLA_A2[:60] + 'W' + HLA_A2[61:100] + 'A' + HLA_A2[101:]
    out = tcrkit.annotate_mhc({'exact': HLA_A2, 'tagged': 'MGSSHHHHHH' + HLA_A2,
                               'mutant': mutant})
    assert out.loc['exact', 'allele_identity'] == 1.0
    assert out.loc['tagged', 'allele_identity'] == 1.0      # the tag is outside the span
    assert 0.95 < out.loc['mutant', 'allele_identity'] < 1.0
    assert out.loc['mutant', 'gene'] == 'A'


def test_topology_is_carried_onto_the_sequence():
    import tcrkit

    full = tcrkit.mhc_sequence('A*0201', region='full')['sequence'][0]
    out = tcrkit.annotate_mhc(full).iloc[0]
    assert out['regions'] == ('signal_peptide, extracellular, transmembrane, '
                              'cytoplasmic')
    assert (out['signal_peptide'], out['extracellular'], out['transmembrane'],
            out['cytoplasmic']) == ('1-24', '25-307', '308-332', '333-365')

    # a fragment is marked for the part it covers, with no collapsed leftovers
    part = tcrkit.annotate_mhc(full[:200]).iloc[0]
    assert (part['signal_peptide'], part['extracellular']) == ('1-24', '25-200')
    assert part['transmembrane'] is None and part['cytoplasmic'] is None


def test_annotate_mhc_can_mark_regions():
    import tcrkit

    seq = tcrkit.mhc_sequence('A0101', region='full')['sequence'][0]
    out = tcrkit.annotate_mhc(seq, mark_regions=True).iloc[0]
    assert len(out['topology_seq']) == len(seq)
    assert out['topology_seq'].startswith('S' * 24 + 'O')
    assert set(out['topology_seq']) == set('SOMI')
    assert 'topology_seq' not in tcrkit.annotate_mhc(seq).columns


# ---------------------------------------------------------------------------
# MHC genes
# ---------------------------------------------------------------------------

def test_mhc_database_order():
    db = mhc.mhc_database()
    assert list(db.columns) == ['allele', 'sequence', 'gene', 'class']
    assert db['gene'].iloc[-1] == 'B2M' and {'A', 'E', 'DRB1'} <= set(db['gene'])


def test_best_mhc():
    gene, allele, score = mhc.best_mhc(HLA_A2)
    assert (gene, allele, round(score, 2)) == ('A', 'A0201', 0.53)
    assert mhc.best_mhc(B2M)[0] == 'B2M'


def test_mhc_span_leaves_out_a_tag():
    assert mhc.mhc_span('GGGGSGGGGS' + HLA_A2) == (10, 10 + len(HLA_A2))


def test_mhc_chain():
    assert [mhc.mhc_chain(g) for g in ['A', 'K', 'DRA', 'DQA1', 'AA', 'EA']] == \
        ['heavy', 'heavy', 'alpha', 'alpha', 'alpha', 'alpha']
    assert [mhc.mhc_chain(g) for g in ['DRB1', 'DQB1', 'AB', 'EB', 'B2M']] == \
        ['beta', 'beta', 'beta', 'beta', 'B2M']
    assert mhc.mhc_chain('H2-K') is None and mhc.mhc_chain(None) is None


def test_mouse_genes():
    db = mhc.mhc_database()
    seq = dict(zip(db['allele'], db['sequence']))
    assert mhc.best_mhc(seq['H2-ABb'])[0] == 'AB' and mhc.best_mhc(seq['H2-AAb'])[0] == 'AA'
    assert mhc.mhc_chain(mhc.best_mhc(seq['H2-Kb'])[0]) == 'heavy'


def test_two_stage_search_matches_a_full_scan():
    """The gene/allele search skips most of the database; it must not skip the answer."""
    from tcrkit.mhc import _best_match, _local_aligner

    db = mhc.mhc_database()
    aligner = _local_aligner()

    def full_scan(seq):
        best_i, best = None, 0.0
        for i, ref in enumerate(db['sequence']):
            score = aligner.score(seq, ref) / max(len(seq), len(ref))
            if score > best:
                best_i, best = i, score
        return db['gene'].iat[best_i], db['allele'].iat[best_i]

    probes = list(db['sequence'][::17]) + [s[50:200] for s in db['sequence'][::41]]
    for seq in probes:
        assert _best_match(seq, db)[:2] == full_scan(seq)


def test_the_leader_is_numbered_backwards_from_imgt_1():
    import tcrkit

    full = tcrkit.mhc_sequence('A0101', region='full')['sequence'][0]
    res = tcrkit.annotate_mhc(full, per_residue=True)

    leader = res[res['topology'] == 'signal_peptide']
    assert list(leader['imgt']) == [str(-n) for n in range(24, 0, -1)]   # -24 ... -1
    assert res.iloc[24]['imgt'] == '1'            # and IMGT 1 follows immediately

    # nothing after the G domains is numbered: those are other IMGT schemes
    after = res[res['topology'].isin(['transmembrane', 'cytoplasmic'])]
    assert (after['imgt'] == '-').all()


def test_the_leader_is_not_numbered_without_a_signal_peptide():
    import tcrkit

    # the mature chain starts at IMGT 1, so there is nothing before it to number
    mature = tcrkit.mhc_sequence('A0101', region='mature')['sequence'][0]
    res = tcrkit.annotate_mhc(mature, per_residue=True)
    assert not [v for v in res['imgt'] if v.startswith('-') and v != '-']


def test_per_residue_topology_labels_every_residue():
    import tcrkit

    full = tcrkit.mhc_sequence('A*0201', region='full')['sequence'][0]
    res = tcrkit.annotate_mhc(full, per_residue=True)
    assert res['topology'].value_counts().to_dict() == {
        'extracellular': 283, 'cytoplasmic': 33, 'transmembrane': 25,
        'signal_peptide': 24}
    assert res.iloc[0]['topology'] == 'signal_peptide'
    assert res.iloc[307]['topology'] == 'transmembrane'


def test_b2m_and_non_mhc(wide):
    assert wide.loc['b2m', ['gene', 'mhc_class']].tolist() == ['B2M', '-']
    assert pd.isna(wide.loc['b2m', 'domains']) and pd.isna(wide.loc['b2m', 'error'])
    assert pd.isna(wide.loc['tcr', 'gene'])
    assert 'no IMGT G-ALPHA1 / G-ALPHA2 / G-ALPHA / G-BETA domain' in wide.loc['tcr', 'error']
    assert 'too short' in wide.loc['pep', 'error']
    assert wide['error'].isna().sum() == 5


def test_per_residue_and_squeeze():
    res = tcrkit.annotate_mhc([HLA_A2, B2M], per_residue=True)
    assert list(res.columns) == ['position', 'aa', 'gene', 'allele', 'domain', 'imgt', 'region',
                                 'helix', 'topology']
    assert len(res) == len(HLA_A2) + len(B2M) and list(res.index.unique()) == [0, 1]
    assert res.loc[0].iloc[49][['aa', 'gene', 'imgt', 'region']].tolist() == \
        ['P', 'A', '50', 'helix']
    assert set(res.loc[1, 'imgt']) == {'-'} and res.loc[1, 'gene'].iloc[0] == 'B2M'
    assert 'Homsap' in res.attrs['templates'][0] and res.attrs['templates'][1] == ''
    one = tcrkit.annotate_mhc(HLA_A2, squeeze=True)
    assert isinstance(one, pd.Series) and one['allele'] == one.allele == 'A0201'
    df = pd.DataFrame({'mhc': [HLA_A2]})
    assert tcrkit.annotate_mhc(df, col_mapping={'sequence': 'mhc'})['domain_spans'][0] == \
        '1-90, 91-182'
    assert tcrkit.annotate_mhc('not a seq!')['error'][0] == 'not an amino-acid sequence'


def test_tpm_bench_helices():
    """The same residues as 'imgt' on class I; a few residues apart on class II."""
    a2 = tcrkit.annotate_mhc(HLA_A2, helix='tpm-bench')
    assert a2['helix'][0] == '50-86, 140-176'
    assert tcrkit.annotate_mhc(IE_BETA, helix='tpm-bench')['helix'][0] == '55-65, 68-92'
    res = tcrkit.annotate_mhc([HLA_A2, B2M], helix='tpm-bench', per_residue=True)
    assert res.loc[0, 'helix'].sum() == 74 and res.loc[1, 'helix'].sum() == 0
    with pytest.raises(ValueError, match='helix must be one of'):
        tcrkit.annotate_mhc(HLA_A2, helix='other')


def test_cli_annotate_mhc(capsys):
    assert main(['annotate-mhc', '--seq', HLA_A2]) == 0
    printed = capsys.readouterr().out
    assert 'A0201' in printed and '50-86, 140-176' in printed
    assert main(['annotate-mhc', '--seq', HLA_A2, '--per-residue']) == 0
    assert 'helix' in capsys.readouterr().out


def test_collapse_into_segments():
    full = tcrkit.mhc_sequence('A*02:01', region='full', squeeze=True)['sequence']
    topo = tcrkit.annotate_mhc(full, collapse='topology')
    assert list(topo['topology']) == ['signal_peptide', 'extracellular', 'transmembrane',
                                      'cytoplasmic']
    assert ''.join(topo['sequence']) == full                    # segments tile the chain
    assert list(topo['start']) == [1, 25, 308, 333]
    assert topo['imgt_end'].iloc[0] == '-1'                     # the leader counts back

    seg = tcrkit.annotate_mhc(HLA_A2, collapse=['domain', 'region'])
    a1 = seg[seg['domain'] == 'G-ALPHA1'].set_index('region')
    assert list(a1.index) == ['A', 'AB', 'B', 'BC', 'C', 'CD', 'D', 'helix']
    assert a1.loc['A', ['imgt_start', 'imgt_end', 'length']].tolist() == ['1', '14', 14]

    with pytest.raises(ValueError, match='collapse takes'):
        tcrkit.annotate_mhc(HLA_A2, collapse='aa')
