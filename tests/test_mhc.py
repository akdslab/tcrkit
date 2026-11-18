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
    assert list(wide.columns) == ['sequence', 'gene', 'allele', 'allele_identity',
                                  'allele_ties', 'allele_source', 'mhc_class',
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
    # human and mouse first, ending with their B2Ms; then the other species
    assert list(db['allele'].iloc[182:184]) == ['B2M_human', 'B2M_mouse']


def test_best_mhc():
    gene, allele, score = mhc.best_mhc(HLA_A2)
    assert (gene, allele, round(score, 2)) == ('A', 'A0201', 0.53)
    assert mhc.best_mhc(B2M)[0] == 'B2M'


#: MHC chains of other species, as deposited: pig SLA-1*04 (3QQ3), chicken BF2*04 (4G42)
#: and chicken BLB2 (7APZ: a peptide-linker in front of the beta chain).
OTHER_SPECIES = {
    'SLA-1': ('pig', 'heavy', 'GPHSLSYFYTAVSRPDRGDSRFIAVGYVDDTQFVRFDNYAPNPRMEPRVPWIQQEGQEYWDRETRNVKE'
              'TAQTYGVGLNTLRGYYNQSEAGSHTLQSMYGCYLGPDGLLLHGYRQDAYDGADYIALNEDLRSWTAADMAAQ'
              'ITKRKWEAADEAERRRSYLQGLCVESLRRYLEMGKDTLQRAEPPKTHVTRHPSSDLGVTLRCWALGFYPKEIS'
              'LTWQREGQDQSQDMELVETRPSGDGTFQKWAALVVPPGEEQSYTCHVQHEGLQEPLTLRWD'),
    'BF2': ('chicken', 'heavy', 'MEFELHTLRYIRTAMTDPGPGQPWFVTVGYVDGELFVHYNSTARRYVPRTEWIAANTDQQYWDGQTQ'
            'IGQLNEQINRENLGIRQRRYNQTGGSHTVQWMFGCDILEDGTIRGYRQSAYDGRDFIALDKDMKTFTAAVPEA'
            'VPTKRKWEEESEPERWKNYLEETCVEWLRRYVEYGKAELGRRERPEVRVWGKEADGILTLSCRAHGFYPRPIV'
            'VSWLKDGAVRGQDAHSGGIVPNGDGTYHTWVTIEAQPGDGDKYQCRVEHASLPQPGLYSWKL'),
    'BLB2': ('chicken', 'beta', 'KMSMSTMNMPMAMKVGGGGSGGGGSGGGGSSAFFFCGAIFECHYLNGTERVRYLQRYIYNRQQLVHF'
             'DSDVGKFVADTPLGEPQAEYWNSNAELLENIMNIADGSCRHNYGILESFTVQRSVEPKVRVSALQSGSLPETDR'
             'LACYVTGFYPPEIEVKWFLNGREETERVVSTDVMQNGDWTYQVLVVLETVPRRGDSYVCRVEHASLRQPISQA'
             'WEPPADAGRSK'),
}


@pytest.mark.parametrize('gene', OTHER_SPECIES)
def test_best_mhc_other_species(gene):
    species, chain, seq = OTHER_SPECIES[gene]
    got, allele, score = mhc.best_mhc(seq)
    assert got == gene and score > 0.07 and mhc.mhc_species(allele) == species
    assert mhc.mhc_chain(gene) == chain


def test_mhc_species():
    assert mhc.mhc_species('A0201') == 'human' and mhc.mhc_species('H2-Kb') == 'mouse'
    assert mhc.mhc_species('B2M_chicken') == 'chicken' and mhc.mhc_species('nope') is None
    # human chains still match human alleles
    assert mhc.mhc_species(mhc.best_mhc(HLA_A2)[1]) == 'human'
    assert mhc.mhc_species(mhc.best_mhc(B2M)[1]) == 'human'


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


# ------------------------------------------- the leader, domains=True (C-LIKE)

def _negatives(res):
    return [v for v in res['imgt'] if v.startswith('-') and v != '-']


def test_class_ii_leader_counts_back_from_the_first_imgt_position():
    # DRA's G-ALPHA starts at IMGT -2 (IMGT's own label), so the leader goes on to -3 ...
    full = tcrkit.mhc_sequence('DRA*01:01', region='full', squeeze=True)['sequence']
    res = tcrkit.annotate_mhc(full, per_residue=True).reset_index(drop=True)
    assert list(res['imgt'][:28]) == [str(-n) for n in range(27, 0, -1)] + ['1']
    assert (res['topology'][:25] == 'signal_peptide').all()

    # DRB1: four mature residues before IMGT 1, numbered on the same run as the leader
    full = tcrkit.mhc_sequence('DRB1*04:01', region='full', squeeze=True)['sequence']
    res = tcrkit.annotate_mhc(full, per_residue=True).reset_index(drop=True)
    assert list(res['imgt'][:34]) == [str(-n) for n in range(33, 0, -1)] + ['1']
    assert res['topology'][:33].value_counts().to_dict() == {'signal_peptide': 29,
                                                             'extracellular': 4}
    mature = tcrkit.mhc_sequence('DRB1*04:01', squeeze=True)['sequence']
    res = tcrkit.annotate_mhc(mature, per_residue=True)
    assert _negatives(res) == ['-4', '-3', '-2', '-1'] and res['aa'].iloc[:4].tolist() == \
        list('GDTR')


def test_the_leader_numbering_leaves_helices_and_domains_alone():
    """Only residues before the first domain residue gain a number; the rest, the
    helix flags included, are what the mature chain gets."""
    for allele in ('DRA*01:01', 'DQA1*05:01', 'DRB1*04:01', 'A*02:01'):
        full = tcrkit.mhc_sequence(allele, region='full', squeeze=True)['sequence']
        mature = tcrkit.mhc_sequence(allele, squeeze=True)['sequence']
        f = tcrkit.annotate_mhc(full, per_residue=True).reset_index(drop=True)
        m = tcrkit.annotate_mhc(mature, per_residue=True).reset_index(drop=True)
        first = int((f['domain'] != '-').to_numpy().argmax())
        assert f['imgt'][:first].str.fullmatch(r'-\d+').all()
        assert (f['domain'][:first] == '-').all() and not f['helix'][:first].any()
        off = full.index(mature)
        cols = ['aa', 'domain', 'imgt', 'region', 'helix']
        assert f[cols].iloc[off:off + len(m)].reset_index(drop=True).equals(m[cols])


def test_b2m_is_numbered_on_imgt_c_like():
    res = tcrkit.annotate_mhc(B2M, per_residue=True, domains=True).reset_index(drop=True)
    cys = res[res['aa'] == 'C']
    assert cys['imgt'].tolist() == ['1023', '1104']    # the C-LIKE disulfide
    assert cys['region'].tolist() == ['B', 'F'] and set(cys['domain']) == {'C-LIKE'}
    assert res['imgt'].iloc[:3].tolist() == ['-2', '-1', '1001']   # IQ before R1
    assert '1045.1' in set(res['imgt']) and '1085.1' in set(res['imgt'])
    assert res.attrs['templates'][0] == 'Homsap B2M 1a1m_B'
    wide = tcrkit.annotate_mhc(B2M, domains=True)
    assert wide[['domains', 'domain_spans']].iloc[0].tolist() == ['C-LIKE', '3-95']
    assert pd.isna(wide['helix'][0]) and pd.isna(wide['error'][0])
    mouse = tcrkit.mhc_sequence('B2M_mouse', squeeze=True)['sequence']
    num = tcrkit.annotate_mhc(mouse, per_residue=True, domains=True)
    assert num.loc[num['aa'] == 'C', 'imgt'].tolist() == ['1023', '1104']


def test_domains_numbers_alpha3_and_class_ii_alpha2_beta2():
    for allele, groove in [('A*02:01', ['G-ALPHA1', 'G-ALPHA2']), ('DRA*01:01', ['G-ALPHA']),
                           ('DRB1*01:01', ['G-BETA'])]:
        seq = tcrkit.mhc_sequence(allele, region='full', squeeze=True)['sequence']
        plain = tcrkit.annotate_mhc(seq, per_residue=True).reset_index(drop=True)
        res = tcrkit.annotate_mhc(seq, per_residue=True, domains=True).reset_index(drop=True)
        assert list(dict.fromkeys(res['domain'].replace('-', None).dropna())) == groove + ['C-LIKE']
        c = res[res['domain'] == 'C-LIKE']
        assert c.loc[c['aa'] == 'C', 'imgt'].tolist() == ['1023', '1104']
        assert c['imgt'].str.fullmatch(r'1\d{3}(\.\d+)?').all()
        # the groove, helices and leader are what they are without domains=True
        other = res['domain'] != 'C-LIKE'
        assert res[other].equals(plain[other])
        assert (plain.loc[~other, 'imgt'] == '-').all()
        wide = tcrkit.annotate_mhc(seq, domains=True)
        assert wide['domains'][0] == ', '.join(groove + ['C-LIKE'])
        assert wide['helix'][0] == tcrkit.annotate_mhc(seq)['helix'][0]
    seg = tcrkit.annotate_mhc(HLA_A2, collapse='domain', domains=True)
    assert seg['domain'].tolist() == ['G-ALPHA1', 'G-ALPHA2']   # the chain stops at alpha2
    assert main(['annotate-mhc', '--seq', HLA_A2, '--domains']) == 0


def test_cdomain_regions_and_templates():
    from tcrkit import cdomain
    assert [cdomain.cdomain_region(x) for x in ['1.2', '15', '15.1', '45.3', '77', '84.2',
                                                '85.2', '96.1', '111.1', '128', '50']] == \
        ['A', 'A', 'AB', 'CD', 'D', 'DE', 'DE', 'EF', 'FG', 'G', '']
    t = cdomain.imgt_cdomain_templates()
    assert {'B2M', 'TRAC', 'TRBC2', 'HLA-A'} <= set(t['gene'])
    assert all(len(s) == len(labs) for s, labs in zip(t['sequence'], t['imgt']))
    slots = cdomain._slots()
    assert slots.index('1.1') < slots.index('1') < slots.index('15') < slots.index('15.1')
    assert slots.index('85.1') < slots.index('85') and '46' not in slots
    with pytest.raises(ValueError, match='no IMGT C-domain template'):
        cdomain.number_cdomain(B2M, genes=['nope'])
