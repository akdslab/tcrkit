# Changelog

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
