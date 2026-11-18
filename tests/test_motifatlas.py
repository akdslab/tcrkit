"""The MHC Motif Atlas source, allele detection, anchors and download_data, offline.

A small stand-in for the four Motif Atlas files is written into a temporary cache, with
motifs planted in its ligands, so nothing is downloaded.
"""

import random
import warnings

import pandas as pd
import pytest

import tcrkit
from tcrkit import motifatlas as ma

AA = 'ACDEFGHIKLMNPQRSTVWY'


def _clear():
    for f in (ma.available, ma.load_motifatlas_table, ma._by_gene,
              ma.load_motifatlas_peptides, ma._background):
        f.cache_clear()
    ma._warned = False


def _peptides(n, length, fixed, rng):
    """`n` random peptides with residues `fixed` ({0-based position: aa}) planted."""
    out = []
    for _ in range(n):
        p = [rng.choice(AA) for _ in range(length)]
        for k, a in fixed.items():
            p[k] = a
        out.append(''.join(p))
    return out


@pytest.fixture
def atlas(tmp_path, monkeypatch):
    rng = random.Random(0)
    b = tcrkit.load_mhc_table().set_index('allele')
    a2 = b.at['A0201', 'sequence'][25:206]              # mature residues 2-182, 181 aa
    b7 = b.at['B0702', 'sequence'][25:206]
    kb = b.at['H2-Kb', 'sequence'][22:203]
    dra, drb = b.at['DRA0101', 'sequence'], b.at['DRB1*0101', 'sequence']
    gapped = drb[:100] + '-' + drb[100:120] + ('A' if drb[120] != 'A' else 'G') + drb[121:]

    folder = tmp_path / 'motifatlas'
    folder.mkdir()
    (folder / 'sequences_I.txt').write_text(
        'Allele\tSequence\n' + ''.join(f'{n}\t{s}\n' for n, s in
                                       [('A0201', a2), ('A0299', a2), ('B0702', b7),
                                        ('H2-Kb', kb)]))
    (folder / 'sequences_II.txt').write_text(
        'Allele\tSequence\n' + ''.join(f'{n}\t{s}\n' for n, s in
                                       [('DRA*01:01:01:01', dra), ('DRB1*01:01', drb),
                                        ('DRB1*01:02', gapped)]))
    rows = [('A0201', p) for p in _peptides(80, 9, {1: 'L', 8: 'V'}, rng)]
    rows += [('A0201', p) for p in _peptides(10, 10, {1: 'L', 9: 'V'}, rng)]
    rows += [('B0702', p) for p in _peptides(80, 9, {1: 'P', 8: 'L'}, rng)]
    (folder / 'peptides_I.txt').write_text(
        'Allele\tPeptide\n' + ''.join(f'{a}\t{p}\n' for a, p in rows))
    rows = []
    for core in _peptides(80, 9, {0: 'Y', 3: 'L', 5: 'A', 8: 'L'}, rng):
        rows.append(('DRB1_01_01', 'GS' + core + 'KRE', core))
    for name in ('DPA1_01_03__DPB1_04_01', 'DPA1_02_01__DPB1_04_01',
                 'DPA1_02_01__DPB1_04_01__reverse'):
        for core in _peptides(60, 9, {0: 'F', 8: 'K'}, rng):
            rows.append((name, 'AA' + core + 'GG', core))
    (folder / 'peptides_II.txt').write_text(
        'Allele\tPeptide\tCore\n' + ''.join(f'{a}\t{p}\t{c}\n' for a, p, c in rows))

    monkeypatch.setenv('TCRKIT_CACHE', str(tmp_path))
    _clear()
    tcrkit.pmhc.load_mhc_table.cache_clear()
    yield {'a2': a2, 'drb': drb}
    _clear()
    tcrkit.pmhc.load_mhc_table.cache_clear()


# ---------------------------------------------------------------------------
# the table, and mhc_sequence(source='motifatlas')
# ---------------------------------------------------------------------------

def test_table_names_normalized_once_with_the_original_kept(atlas):
    t = tcrkit.load_mhc_table(source='motifatlas').set_index('allele')
    assert t.at['DRA0101', 'motifatlas_name'] == 'DRA*01:01:01:01'
    assert t.at['DRB1*0101', 'motifatlas_name'] == 'DRB1*01:01'
    assert t.at['A0201', 'covers'] == 'groove' and t.at['DRB1*0101', 'covers'] == 'full'
    assert '-' not in t.at['DRB1*0102', 'sequence']                 # alignment gap gone
    assert bool(t.at['A0201', 'has_ligands']) and not bool(t.at['A0299', 'has_ligands'])
    assert (ma.cache_dir() / 'table.csv').exists()                  # built once, cached


def test_class_i_only_has_the_groove(atlas):
    out = tcrkit.mhc_sequence(['A*02:01', 'DRB1*01:01'], source='motifatlas')
    assert pd.isna(out.loc[0, 'sequence']) and "region='groove'" in out.loc[0, 'error']
    assert isinstance(out.loc[1, 'sequence'], str) and pd.isna(out.loc[1, 'error'])
    g = tcrkit.mhc_sequence('A*02:01', source='motifatlas', region='groove', squeeze=True)
    assert g['sequence'] == atlas['a2'] and g['length'] == 181
    bundled = tcrkit.mhc_sequence('A*02:01', region='groove', squeeze=True)['sequence']
    assert bundled[1:] == atlas['a2']                    # the Motif Atlas starts at residue 2


@pytest.mark.parametrize('region', ['extracellular', 'mature', 'signal_peptide',
                                    'transmembrane', 'groove', 'full'])
def test_class_ii_regions_are_the_bundled_ones(atlas, region):
    ours = tcrkit.mhc_sequence(['DRA*01:01', 'DRB1*01:01'], region=region)['sequence']
    theirs = tcrkit.mhc_sequence(['DRA*01:01', 'DRB1*01:01'], region=region,
                                 source='motifatlas')['sequence']
    assert ours.tolist() == theirs.tolist()


def test_b2m_comes_from_the_bundled_reference(atlas):
    # the Motif Atlas has no B2M; class I pairs name it, and it is looked up in ours
    pair = tcrkit.mhc_sequence('A*02:01', paired=True, source='motifatlas',
                               region='groove').iloc[0]
    assert pair['mhc_b_allele'] == 'B2M_human'
    assert pd.isna(pair['mhc_b_sequence'])                     # B2M has no groove
    b2m = tcrkit.mhc_sequence('B2M_human', source='motifatlas', region='mature').iloc[0]
    assert b2m['sequence'] == tcrkit.mhc_sequence('B2M_human', region='mature')['sequence'][0]


# ---------------------------------------------------------------------------
# annotate_mhc's allele
# ---------------------------------------------------------------------------

def test_annotate_mhc_reports_the_motif_atlas_allele(atlas):
    full = tcrkit.mhc_sequence('A*02:01', region='full', squeeze=True)['sequence']
    r = tcrkit.annotate_mhc(full, squeeze=True)
    assert (r['allele'], r['allele_source'], r['allele_identity']) == ('A0201', 'motifatlas', 1.0)
    assert r['allele_ties'] == 1                          # A0299 has the same groove
    # numbering, helices and membrane regions do not depend on where the allele is from
    u = tcrkit.annotate_mhc(full, allele_source='uniprot', squeeze=True)
    assert u['allele_source'] == 'uniprot'
    for col in ('gene', 'domains', 'helix', 'extracellular', 'transmembrane'):
        assert r[col] == u[col], col


def test_annotate_mhc_closest_when_nothing_matches_exactly(atlas):
    drb = atlas['drb']
    mutant = drb[:60] + ('W' if drb[60] != 'W' else 'Y') + drb[61:]
    r = tcrkit.annotate_mhc(mutant, squeeze=True)
    assert r['allele'] == 'DRB1*0101' and 0.99 < r['allele_identity'] < 1.0


def test_annotate_mhc_falls_back_offline(tmp_path, monkeypatch):
    monkeypatch.setenv('TCRKIT_CACHE', str(tmp_path))               # empty, no downloads
    _clear()
    seq = tcrkit.mhc_sequence('A*02:01', squeeze=True)['sequence']
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        r = tcrkit.annotate_mhc(seq, squeeze=True)
    assert r['allele_source'] == 'uniprot' and r['allele'] == 'A0201'
    assert any('Motif Atlas not available' in str(x.message) for x in w)
    _clear()


def test_allele_source_is_checked():
    with pytest.raises(ValueError, match='allele_source'):
        tcrkit.annotate_mhc('GSHSMRYF', allele_source='ipd')


# ---------------------------------------------------------------------------
# anchors and motifs
# ---------------------------------------------------------------------------

def test_class_i_anchors_from_the_ligands(atlas):
    out = tcrkit.peptide_anchors(['LLFGYPVYV', 'RPRGEVRFL', 'ALAAAAAAAV'],
                                 ['HLA-A*02:01', 'B*07:02', 'A*02:01'])
    assert out['anchors'].tolist()[:2] == ['P2, P9', 'P2, P9'] and pd.isna(out.loc[2, 'anchors'])
    assert out.loc[0, 'anchor_residues'] == 'L, V' and out.loc[0, 'n_ligands'] == 80
    assert 'needs 50' in out.loc[2, 'error']             # only 10 ten-mers


def test_class_ii_core_and_anchors(atlas):
    r = tcrkit.peptide_anchors('KKYAALTAWLGE', 'DRB1*01:01', squeeze=True)
    assert (r['core'], r['core_start']) == ('YAALTAWLG', 3) or r['core'].startswith('Y')
    assert r['anchors'] == 'P1, P4, P6, P9'
    assert r['anchor_positions'] == ', '.join(str(r['core_start'] + k) for k in (0, 3, 5, 8))


def test_class_ii_pairs_must_be_named_and_reverse_sets_are_left_out(atlas):
    out = tcrkit.peptide_anchors(['AAFKKKKKKKGG'] * 2,
                                 ['DPB1*04:01', 'DPA1*02:01__DPB1*04:01'])
    assert 'name the pair' in out.loc[0, 'error']
    assert out.loc[1, 'motif_allele'] == 'DPA1*0201__DPB1*0401'
    assert out.loc[1, 'n_ligands'] == 60                 # not 120: __reverse excluded


def test_unknown_alleles_are_errors_not_guesses(atlas):
    out = tcrkit.peptide_anchors(['LLFGYPVYV'] * 2, ['B*51:250', 'not-an-allele'])
    assert out['anchors'].isna().all()
    assert 'no Motif Atlas ligands' in out.loc[0, 'error']
    assert 'not an allele name' in out.loc[1, 'error']


def test_peptide_motif(atlas):
    m = tcrkit.peptide_motif('A*02:01')
    assert m.shape == (9, 21) and m.attrs['n_ligands'] == 80
    assert m['ic'].idxmax() in (2, 9)
    assert m[list(AA)].sum(axis=1).round(6).eq(1).all()
    assert m.loc[2].drop('ic').idxmax() == 'L'


def test_plot_motif_needs_logomaker(atlas):
    try:
        import logomaker  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match=r"tcrkit\[logos\]"):
            tcrkit.plot_motif('A*02:01')
    else:
        import matplotlib
        matplotlib.use('Agg')
        assert tcrkit.plot_motif('A*02:01') is not None


# ---------------------------------------------------------------------------
# download_data
# ---------------------------------------------------------------------------

def test_download_data_reports_each_source(atlas, monkeypatch):
    from tcrkit import download

    monkeypatch.setattr(download, '_stitchr',
                        lambda sp, force: {'status': 'cached', 'path': f'/x/{sp}'})

    def boom(force):
        raise OSError('no network')

    monkeypatch.setattr(download, '_decodertcr', boom)
    out = tcrkit.download_data()
    assert out['source'].tolist() == ['stitchr', 'stitchr', 'motifatlas', 'decodertcr']
    assert out['status'].tolist() == ['cached', 'cached', 'cached', 'failed']
    assert 'no network' in out['error'].iat[-1]
    with pytest.raises(ValueError, match='sources'):
        tcrkit.download_data(['imgt'])


def test_cli_download_parses():
    from tcrkit.cli import build_parser

    args = build_parser().parse_args(['download', '--only', 'motifatlas', '--force'])
    assert args.only == ['motifatlas'] and args.force
