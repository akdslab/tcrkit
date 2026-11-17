---
name: tcrkit-mhc
description: Using tcrkit on MHC/HLA - normalize allele names, get MHC protein sequences or regions, annotate an MHC chain's gene, IMGT numbering, groove helices and membrane regions, choose a sequence source (bundled UniProt or DecoderTCR). Use when writing analysis code with tcrkit for MHC.
---

# tcrkit for MHC

```python
tcrkit.normalize_mhc(['A*0201', 'H2-Kb'])                     # allele names
tcrkit.mhc_sequence('HLA-A*02:01')                            # extracellular domain
tcrkit.mhc_sequence('A*02:01', region='full', paired=True)    # + B2M; regions below
tcrkit.mhc_sequence('DRA*01:02', source='decodertcr')         # other source
tcrkit.annotate_mhc(seq)                                      # gene, allele, domains, helix
tcrkit.annotate(seq)['chain_type']     # any chain: TCR_alpha/beta, MHC_I, MHC_II_alpha/beta, B2M, peptide
tcrkit.annotate_mhc(seq, collapse=True)                       # strands A-D, turns, helix
tcrkit.annotate_mhc(seq, collapse='topology')                 # signal peptide .. tail
tcrkit.annotate_mhc(seq, per_residue=True, domains=True)      # + alpha3/beta2/B2M, C-LIKE + 1000
tcrkit.annotate_mhc(seq, show=True)                           # prints it (colour) (print_residues); returns None
tcrkit.list_mhc('genes', species='human'); tcrkit.list_mhc('pairs', source='decodertcr')
tcrkit.imgt_regions()          # IMGT ranges of CDRs + helices, with sources (all_regions=True)
```

- `region=`: `extracellular` (default), `full`, `mature`, `signal_peptide`,
  `transmembrane`, `cytoplasmic`. Regions mean the same stretch of the chain for both
  sources (DecoderTCR chains get the bundled topology).
- Sources: bundled = 186 human / mouse alleles from UniProt, plus 76 chains of TCR3d's
  20 other species (gene and species level: `mhc_species(best_mhc(seq)[1])`), shipped; `decodertcr` = 87 human chains
  / 88 pairs, downloaded once and cached (`~/.cache/tcrkit`). Misses are in `error`.
- `annotate_mhc`: trust `gene`; `allele` is the closest bundled one
  (`allele_identity` < 1.0 = not bundled). `helix='imgt'` (default) or `'tpm-bench'`.
  The `helix` column (docking-angle residues, IMGT 50/52-86) is narrower than
  `region == 'helix'` (IMGT's whole 50-92). Use the `helix` column as the reference
  for the groove helices; "helix" elsewhere means the column.
  The N-terminus is numbered backwards from the first IMGT position (class I leader
  `-24`..`-1`; DRA continues from IMGT -2; DRB1's 4 extra mature residues are `-4`..`-1`);
  `topology` says which are signal peptide. α3 / α2 / β2 and B2M stay `'-'`
  unless `domains=True`: then IMGT C-LIKE + 1000 (`domain == 'C-LIKE'`, Cys
  `1023`/`1104`). The membrane anchor and tail always stay `'-'`.
- Full reference: `tcrkit docs api` (prints `docs/api.md`).
