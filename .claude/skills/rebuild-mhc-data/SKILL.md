---
name: rebuild-mhc-data
description: Refresh tcrkit's bundled MHC sequences from UniProt and their DeepTMHMM topology. Use when asked to update, rebuild or add MHC alleles/sequences in src/tcrkit/data.
---

# Rebuild the MHC data

Files: `src/tcrkit/data/mhc_uniprot.fasta` (UniProt, verbatim) and `mhc_annotations.csv`
(one row per allele; the `source`/`accession` columns drive everything).

1. Add or change alleles by editing rows of `mhc_annotations.csv` - `accession` must be a
   UniProtKB (`source=uniprot`) or UniParc (`source=uniparc`) entry *identical* to the
   allele. No IPD-IMGT/HLA data in the package.
2. `python -m tcrkit._build.fetch_mhc_sequences --check`, then without `--check`.
3. `python -m tcrkit._build.annotate_mhc_deeptmhmm --missing-only` (needs `pybiolib`;
   DeepTMHMM is pinned to 1.0.24). Full rerun: drop `--missing-only`.
4. New allele that should be matched by `annotate_mhc`: set `in_detect`, `detect_rank`,
   `detect_class`.
5. `uv run pytest`, then re-run `notebooks/tcrkit_showcase.ipynb`. Update counts in
   `docs/api.md`, `pmhc.py` and the data README if the allele count changed.
