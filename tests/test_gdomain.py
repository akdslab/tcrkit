"""IMGT numbering of MHC groove domains (tcrkit.gdomain)."""

import pytest

from tcrkit import gdomain

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
B2M = ('IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYA'
       'CRVNHVTLSQPKIVKWDRDM')


# ---------------------------------------------------------------------------
# the IMGT data
# ---------------------------------------------------------------------------

def test_imgt_templates_are_shipped_and_ordered():
    t = gdomain.imgt_gdomain_templates()
    assert set(t['domain']) == set(gdomain.GDOMAINS)
    assert set(t['species']) == {'Homsap', 'Musmus', 'Ratnor'} and len(t) == 42
    for row in t.itertuples():
        assert len(row.imgt) == len(row.sequence)
        keys = [gdomain._slot_key(lab) for lab in row.imgt]
        assert keys == sorted(keys) and len(set(keys)) == len(keys)   # IMGT order, no repeats
    a2 = t[(t['gene'] == 'MH1-A') & (t['domain'] == 'G-ALPHA1')].iloc[0]
    # HLA-A2 alpha1: P50 starts the helix, 54A-C are gaps, Q55 follows E54
    assert dict(zip(a2['imgt'], a2['sequence']))['50'] == 'P'
    assert a2['imgt'][a2['imgt'].index('54') + 1] == '55'


@pytest.mark.parametrize('label, region', [('1', 'A'), ('7A', 'A'), ('16', 'AB'), ('49', 'D'),
                                           ('49.3', ''), ('50', 'helix'), ('92B', 'helix'),
                                           ('-2', ''), ('-', '')])
def test_gdomain_region(label, region):
    assert gdomain.gdomain_region(label) == region


@pytest.mark.parametrize('domain, label, on', [('G-ALPHA1', '50', True), ('G-ALPHA1', '86', True),
                                               ('G-ALPHA1', '87', False), ('G-ALPHA2', '51', False),
                                               ('G-ALPHA2', '52', True), ('G-BETA', '61A', True),
                                               ('G-BETA', '61B', False), ('G-BETA', '62', False),
                                               ('G-BETA', '63', True), ('-', '-', False)])
def test_is_helix(domain, label, on):
    assert gdomain.is_helix(domain, label) is on


# ---------------------------------------------------------------------------
# numbering one sequence
# ---------------------------------------------------------------------------

def test_number_gdomains_reproduces_a_template():
    a2 = gdomain.imgt_gdomain_templates().query("gene == 'MH1-A' and domain == 'G-ALPHA1'").iloc[0]
    num = gdomain.number_gdomains(a2['sequence'])
    assert list(num['imgt']) == a2['imgt'] and set(num['domain']) == {'G-ALPHA1'}
    assert num.loc[num['imgt'] == '50', 'region'].item() == 'helix'


def test_number_gdomains_class_i():
    num = gdomain.number_gdomains(HLA_A2)
    # HLA-A2 numbering is IMGT numbering in alpha1; alpha2 starts again at 1
    assert num['imgt'].iloc[49] == '50' and num['imgt'].iloc[90] == '1'
    assert list(num['domain'].iloc[[0, 89, 90]]) == ['G-ALPHA1', 'G-ALPHA1', 'G-ALPHA2']
    helix = num.index[num['helix']] + 1                     # 1-based
    assert list(helix) == list(range(50, 87)) + list(range(140, 177))


def test_number_gdomains_is_not_shared_between_calls():
    num = gdomain.number_gdomains(HLA_A2)
    num.loc[0, 'imgt'] = 'changed'
    assert gdomain.number_gdomains(HLA_A2)['imgt'].iloc[0] == '1'


def test_no_g_domain_in_b2m():
    num = gdomain.number_gdomains(B2M)
    assert set(num['domain']) == {'-'} and max(num.attrs['scores'].values()) < 0.4
    with pytest.raises(ValueError, match='unknown G domains'):
        gdomain.number_gdomains('ACDE', ['G-GAMMA'])


def test_residues_two_domains_claim_are_split():
    """4MS8 lacks residues where alpha1 meets alpha2. Both alignments claim G, Y and H
    there; G and Y are alpha1's 83-84 and H is alpha2's 3."""
    num = gdomain.number_gdomains(H2_LD_GAPPED, ['G-ALPHA1', 'G-ALPHA2'])
    assert set(num['domain']) == {'-', 'G-ALPHA1', 'G-ALPHA2'}
    part = num.iloc[77:80]
    assert ''.join(part['aa']) == 'GYH'
    assert list(zip(part['domain'], part['imgt'])) == [('G-ALPHA1', '83'), ('G-ALPHA1', '84'),
                                                       ('G-ALPHA2', '3')]


def test_numbering_continues_past_a_divergent_end():
    """7RDV's beta helix differs from every IMGT entry at its end; the numbering still
    reaches 86, the last groove-helix residue."""
    num = gdomain.number_gdomains(IE_BETA, ['G-BETA'])
    helix = num[num['region'] == 'helix']
    assert helix['imgt'].iloc[-1] == '89' and num['helix'].sum() == 36


def test_tpm_bench_helix():
    assert gdomain.tpm_bench_helix(HLA_A2).sum() == 74
    with pytest.raises(ValueError, match='too short'):
        gdomain.tpm_bench_helix('GSHSMRYF')


@pytest.mark.parametrize('start', [40, 50, 60, 100, 150])
def test_tpm_bench_helix_of_a_chain_starting_inside_a_helix(start):
    """A chain that starts inside a helix keeps the part of it that it has; the alpha1
    helix (HLA-A2 50-86) used to be dropped whole."""
    full = gdomain.tpm_bench_helix(HLA_A2)
    assert (gdomain.tpm_bench_helix(HLA_A2[start:]) == full[start:]).all()
