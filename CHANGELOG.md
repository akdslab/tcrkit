# Changelog

## 0.1.9 - 2026-11-14

- MHC chains of 20 more species (macaque, pig, cattle, chicken, duck, ...: those of
  TCR3d's structures): 76 UniProt entries, matched at gene and species level. New
  `tcrkit.mhc.mhc_species(allele)`; `load_mhc_table()` has an `organism` column.
  Human and mouse matches are unchanged. `list_mhc(species=...)` keeps only some species.
- New `germline_species` column in `anarci`: the V germline's species. `species` is
  the HMM's, which can differ (2E7L's mouse beta chain: HMM human, germline mouse).
- `anarci(..., include_imgt_index=True)` always has the `imgt_dict` column, also when
  no sequence is numbered (a pMHC).

## 0.1.8 - 2026-11-13

- The docs and the three user skills come with the wheel: `tcrkit.docs()` / `tcrkit docs
  [page]`, `tcrkit.install_skills()` / `tcrkit skills install`. Nothing else changes.

## 0.1.7 - 2026-11-12

- `print_residues` / `show=` legend: strands grouped, items in sequence order, wrapped at
  80 characters instead of one long line per annotation.
- README: coloured `show=True` images of a TCR and an MHC; shorter citation table (full
  references in the data README).

## 0.1.6 - 2026-11-12

- **Class II leaders numbered**: everything before the first IMGT position counts back
  from it (DRA `-27..-3`, DRB1 `-33..-1`), as class I already did. Domains, regions and
  helices are unchanged.
- New `domains=True` on `annotate_mhc` (CLI `--domains`): α3 / α2 / β2 and B2M numbered
  on IMGT C-LIKE + 1000 (`domain='C-LIKE'`); new `tcrkit.cdomain`, `data/imgt_cdomains.csv`.
- New `domains=True` on `anarci` (CLI `--domains`): `domain` (V/C) and `gene_segment`
  (L/V/J/N-D/C, L = the V gene's leader) per residue, the TCR constant domain on IMGT C-DOMAIN + 1000, and
  `domains` / `domain_spans` per sequence.
- New `mark_regions=True` on `anarci` (CLI `--mark-regions`): `region_seq` (and
  `segment_seq`), one character per residue. `print_residues` shows `gene_segment`.
- New `show=True` on `anarci` / `annotate_mhc`: prints the `print_residues` view
  (returns None).
  `print_residues(width=None)`: one line per sequence.

## 0.1.5 - 2026-11-04

- **Stitched nucleotides are opt-in and renamed**: `sequence` → `sequence_nt`
  (`{chain}_sequence_nt`), returned only with `include_nt=True` on `stitch`,
  `stitch_pair` and `build_tcr_pmhc` (CLI `--include-nt`).
- `squeeze=` on `to_gene`, `stitch_pair` and `build_tcr_pmhc`.
- New `imgt_regions`: the IMGT ranges behind `region` / `helix`, with sources.
- `print_residues(color=...)`: residues and a line coloured by its values.

## 0.1.4 - 2026-11-03

- **MHC sequences now come from UniProt**: 186 alleles, each identical to its entry; 65
  without an exact entry are dropped. Gene, domain and helix calls are unchanged.
- **`squeeze=True` returns a Series, not a dict**, for one input with several columns.
- New `source='decodertcr'` on `mhc_sequence` / `load_mhc_table`, fetched once and cached.
- New `list_mhc`, `print_residues`, `anarci(per_residue=True)`, and `collapse=` on
  `anarci` / `annotate_mhc` (also `--collapse`).
- Data rebuild scripts moved to `tcrkit._build`.

## 0.1.3 - 2026-11-03

- **Removed `detect_chains`**; it lives in tcrstruct now.
- Fixed: `helix='tpm-bench'` dropped a helix when the chain started inside it.

## 0.1.2 - 2026-11-03

- **`mhc_sequence(region=...)` replaces `full=True`.**
- `annotate_mhc` reports membrane regions and `allele_identity`.
- MHC sequences include their signal peptide; the MHC data is one bundled table.
- New CLI subcommands: `diversity`, `morisita-horn`, `generation-probability`.
- Fixed: mouse class II alpha chains were numbered as G-BETA.

## 0.1.1 - 2026-10-30

- New `annotate_mhc` and `detect_chains`, also on the command line.

## 0.1.0

First release.
