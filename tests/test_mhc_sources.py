"""The MHC sequence reference and its sources: bundled UniProt, and DecoderTCR.

The DecoderTCR tests run offline: a small stand-in for its ``hla_reference.json`` is
written to a temporary cache, so nothing is downloaded.
"""

import json

import pandas as pd
import pytest

import tcrkit
from tcrkit import _decodertcr
from tcrkit._build import annotate_mhc_deeptmhmm as deeptmhmm
from tcrkit._build import fetch_mhc_sequences as fetch


# ---------------------------------------------------------------------------
# the bundled reference: mhc_uniprot.fasta + mhc_annotations.csv
# ---------------------------------------------------------------------------

def test_bundled_table_joins_the_fasta_and_the_annotations():
    t = tcrkit.load_mhc_table()
    assert len(t) == 186 and t['allele'].is_unique and t['accession'].is_unique
    assert set(t['source']) == {'uniprot', 'uniparc'}
    assert t['sequence'].notna().all()
    assert (t['length'] == t['sequence'].str.len()).all()
    assert t['in_detect'].sum() == 184
    # what annotate_mhc / mhc_sequence / the tpm-bench helices read
    for col in ('signal_spans', 'outside_spans', 'tmhelix_spans', 'inside_spans',
                'prediction', 'detect_rank', 'detect_class', 'helix_ind',
                'helix_ref_sequence', 'mature_start'):
        assert col in t.columns, col
    a2 = t.set_index('allele').loc['A0201']
    assert a2['sequence'].startswith('MAVMAPRTL') and a2['mature_start'] == 24
    assert t['helix_ref_sequence'].notna().sum() == 18


def test_bundled_path_gives_the_bundled_table():
    assert tcrkit.load_mhc_table(tcrkit.bundled_mhc_table_path()).equals(tcrkit.load_mhc_table())


def test_unknown_source_is_refused():
    with pytest.raises(ValueError, match='source must be one of'):
        tcrkit.load_mhc_table(source='nope')
    with pytest.raises(ValueError, match='not both'):
        tcrkit.mhc_sequence('A*02:01', table=tcrkit.load_mhc_table(), source='uniprot')


# ---------------------------------------------------------------------------
# DecoderTCR, offline
# ---------------------------------------------------------------------------

@pytest.fixture
def decodertcr(tmp_path, monkeypatch):
    """A stand-in DecoderTCR reference: A*02:01/B2M and DRA*01:01/DRB1*01:01, mature."""
    bundled = tcrkit.load_mhc_table().set_index('allele')
    mature = {a: bundled.at[a, 'sequence'][int(bundled.at[a, 'mature_start']):]
              for a in ('A0201', 'B2M_human', 'DRA0101', 'DRB1*0101')}
    fake = {'A0201__B2M': {'HLA_a': mature['A0201'], 'HLA_b': mature['B2M_human']},
            'DRA0101__DRB10101': {'HLA_a': mature['DRA0101'], 'HLA_b': mature['DRB1*0101']}}
    path = tmp_path / 'hla_reference.json'
    path.write_text(json.dumps(fake))
    monkeypatch.setattr(_decodertcr, 'fetch_decodertcr', lambda force=False: path)
    tcrkit.pmhc.load_mhc_table.cache_clear()
    yield fake
    tcrkit.pmhc.load_mhc_table.cache_clear()


def test_decodertcr_pairs_become_single_chains(decodertcr):
    t = tcrkit.load_mhc_table(source='decodertcr')
    assert sorted(t['allele']) == ['A0201', 'B2M_human', 'DRA0101', 'DRB1*0101']
    assert set(t['source']) == {'decodertcr'} and not t['in_detect'].any()
    assert t.set_index('allele').at['DRB1*0101', 'pairs'] == 'DRA0101__DRB10101'


def test_decodertcr_chains_alone_or_paired(decodertcr):
    one = tcrkit.mhc_sequence('DRB1*01:01', source='decodertcr').iloc[0]
    assert one['allele'] == 'DRB1*0101' and pd.isna(one['error'])
    pair = tcrkit.mhc_sequence('HLA-A*02:01', paired=True, source='decodertcr').iloc[0]
    assert (pair['mhc_a_allele'], pair['mhc_b_allele']) == ('A0201', 'B2M_human')
    assert pair['mhc_b_sequence']


@pytest.mark.parametrize('region', ['extracellular', 'mature', 'transmembrane', 'cytoplasmic'])
def test_regions_mean_the_same_in_both_sources(decodertcr, region):
    """DecoderTCR's chains get their topology from the bundled allele, so a region cut
    from either source is the same stretch of the chain."""
    alleles = ['A0201', 'DRA0101', 'DRB1*0101']
    ours = tcrkit.mhc_sequence(alleles, region=region)['sequence'].tolist()
    theirs = tcrkit.mhc_sequence(alleles, region=region, source='decodertcr')['sequence'].tolist()
    assert ours == theirs


def test_decodertcr_imgt_numbering_matches_the_bundled_sequence(decodertcr):
    full = tcrkit.mhc_sequence('A0201', region='full', squeeze=True)['sequence']
    mature = tcrkit.mhc_sequence('A0201', region='full', source='decodertcr',
                                 squeeze=True)['sequence']
    a = tcrkit.annotate_mhc(full, collapse=['domain', 'region'])
    b = tcrkit.annotate_mhc(mature, collapse=['domain', 'region'])
    a, b = a[a['domain'] != '-'], b[b['domain'] != '-']
    assert a[['domain', 'region', 'imgt_start', 'imgt_end', 'sequence']].values.tolist() == \
        b[['domain', 'region', 'imgt_start', 'imgt_end', 'sequence']].values.tolist()


def test_decodertcr_download_is_checked_and_cached(tmp_path, monkeypatch):
    """A file with the wrong checksum is refused; a good one is downloaded once."""
    monkeypatch.setenv('TCRKIT_CACHE', str(tmp_path))
    payload = b'{"A0201__B2M": {"HLA_a": "GSH", "HLA_b": "IQR"}}'
    calls = []

    class Response:
        def __init__(self, data): self.data = data
        def read(self): return self.data
        def __enter__(self): return self
        def __exit__(self, *a): return False

    def fake_urlopen(url, timeout=None):
        calls.append(url)
        return Response(payload)

    monkeypatch.setattr(_decodertcr.urllib.request, 'urlopen', fake_urlopen)
    with pytest.raises(RuntimeError, match='SHA-256'):
        _decodertcr.fetch_decodertcr()
    assert not (_decodertcr.cache_dir() / 'hla_reference.json').exists()

    import hashlib
    monkeypatch.setattr(_decodertcr, 'SHA256', hashlib.sha256(payload).hexdigest())
    path = _decodertcr.fetch_decodertcr()
    assert path.read_bytes() == payload and _decodertcr.COMMIT in str(path)
    _decodertcr.fetch_decodertcr()
    assert len(calls) == 2                    # the refused one, then once - not again


# ---------------------------------------------------------------------------
# the _build scripts' parsers
# ---------------------------------------------------------------------------

def test_build_scripts_parse_their_sources():
    assert fetch.accession('>sp|P61769|B2MG_HUMAN Beta-2-microglobulin') == 'P61769'
    assert fetch.accession('>UPI00005FF33A status=active') == 'UPI00005FF33A'
    recs = fetch.parse_fasta('>sp|P1|X\nMAA\nGG\n>UPI1 status=active\nKK\n')
    assert [(a, r) for a, _, r in recs] == [('P1', '>sp|P1|X\nMAA\nGG\n'),
                                           ('UPI1', '>UPI1 status=active\nKK\n')]


def test_build_deeptmhmm_parser(tmp_path):
    path = tmp_path / 'predicted_topologies.3line'
    path.write_text('>sp|P1|X | SP+TM\nMAAAGGGKKL\nSSSOOOMMII\n'
                    '>UPI2 status=active | GLOB\nKKKK\nOOOO\n')
    out = deeptmhmm.parse_3line(path).set_index('accession')
    assert out.loc['P1', ['signal_spans', 'outside_spans', 'tmhelix_spans', 'inside_spans',
                          'mature_start']].tolist() == ['[(1, 3)]', '[(4, 6)]', '[(7, 8)]',
                                                         '[(9, 10)]', 3]
    assert out.loc['UPI2', 'deeptmhmm_prediction'] == 'GLOB'
    assert out.loc['UPI2', 'mature_start'] == 0


# ---------------------------------------------------------------------------
# list_mhc, print_residues, collapse=True
# ---------------------------------------------------------------------------

def test_list_mhc_bundled():
    alleles = tcrkit.list_mhc()
    assert len(alleles) == 186 and list(alleles.columns) == ['allele', 'gene', 'mhc_class', 'species']
    genes = tcrkit.list_mhc('genes').set_index('gene')
    assert genes['n_alleles'].sum() == 186
    assert 'A0201' in genes.at['A', 'alleles']
    with pytest.raises(ValueError, match='single chains'):
        tcrkit.list_mhc('pairs')
    with pytest.raises(ValueError, match="what must be"):
        tcrkit.list_mhc('loci')


def test_list_mhc_decodertcr_pairs(decodertcr):
    pairs = tcrkit.list_mhc('pairs', source='decodertcr')
    assert pairs.values.tolist() == [['A0201__B2M', 'A0201', 'B2M_human', 'I'],
                                     ['DRA0101__DRB10101', 'DRA0101', 'DRB1*0101', 'II']]
    assert sorted(tcrkit.list_mhc(source='decodertcr')['allele']) == \
        ['A0201', 'B2M_human', 'DRA0101', 'DRB1*0101']


def test_collapse_true_uses_the_default_columns():
    seq = tcrkit.mhc_sequence('A*02:01', region='full', squeeze=True)['sequence']
    assert tcrkit.annotate_mhc(seq, collapse=True).equals(
        tcrkit.annotate_mhc(seq, collapse=['domain', 'region']))
    assert len(tcrkit.annotate_mhc(seq, collapse=False)) == 1        # the summary row


def test_print_residues(capsys):
    seq = tcrkit.mhc_sequence('A*02:01', region='full', squeeze=True)['sequence']
    tcrkit.print_residues(tcrkit.annotate_mhc({'A2': seq}, per_residue=True), width=50)
    out = capsys.readouterr().out
    lines = out.splitlines()
    assert lines[0] == f'A2  {len(seq)} residues'
    assert lines[2] == 'sequence  ' + seq[:50]
    assert lines[3] == 'topology  ' + 'S' * 24 + 'O' * 26             # signal peptide 1-24
    assert 'topology: O extracellular, M transmembrane' in out or 'M transmembrane' in out
    with pytest.raises(ValueError, match='per-residue table'):
        tcrkit.print_residues(tcrkit.annotate_mhc(seq))


def test_print_residues_color(capsys):
    seq = tcrkit.mhc_sequence('A*02:01', squeeze=True)['sequence']
    res = tcrkit.annotate_mhc({'A2': seq}, per_residue=True)
    tcrkit.print_residues(res)
    plain = capsys.readouterr().out
    assert '\033[' not in plain                                   # off by default
    tcrkit.print_residues(res, color=True)
    colored = capsys.readouterr().out
    assert '\033[34m' in colored and '\033[31m' in colored        # strands, helix
    import re
    assert re.sub(r'\033\[\d+m', '', colored) == plain            # same text otherwise
    with pytest.raises(ValueError, match='no .topology. line'):
        tcrkit.print_residues(res, columns=['region'], color='topology')


def test_imgt_regions_lists_the_ranges_behind_the_labels():
    t = tcrkit.imgt_regions()                       # by default: CDRs and helices only
    assert t['region'].tolist() == ['CDR1', 'CDR2', 'CDR3'] + ['helix'] * 4
    assert set(t.loc[t['scheme'] == 'mhc', 'column']) == {'helix'}
    full = tcrkit.imgt_regions(all_regions=True)
    assert set(full['scheme']) == {'tcr', 'mhc'} and full['source'].str.len().gt(0).all()
    tcr = tcrkit.imgt_regions('tcr').set_index('region')
    assert tuple(tcr.loc['CDR3', ['imgt_start', 'imgt_end']]) == (105, 117)
    b = tcrkit.imgt_regions('mhc', all_regions=True).query("domain == 'G-BETA' and region == 'helix'")
    assert b[['column', 'imgt_start', 'imgt_end', 'skipped']].values.tolist() == \
        [['region', 50, 92, ''], ['helix', 52, 86, '61B, 62']]
    with pytest.raises(ValueError, match='scheme'):
        tcrkit.imgt_regions('bcr')
