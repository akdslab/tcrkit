---
name: tcrkit-tcr
description: Using tcrkit on TCR sequences - stitch full chains from V/J/CDR3, repair CDR3 junctions, ANARCI/IMGT annotation, gene-name normalization, building whole TCR-pMHC tables. Use when writing analysis code with tcrkit for TCRs.
---

# tcrkit for TCRs

All functions take a value, list, dict or DataFrame and return a DataFrame
(`col_mapping={'v': 'v_call', ...}` for your column names; `squeeze=True` for the value).
Check the `error` column - bad rows are reported, not raised.

```python
import tcrkit
tcrkit.normalize_tcr(['TCRBV20S1'])                       # names -> IMGT ('TRBV20-1')
tcrkit.add_cdr3_junctions('ASSIRSSYEQY', 'TRBJ2-7')       # restore C104/F118 first...
tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta')   # ...then stitch
tcrkit.anarci(seq, allow={'A', 'B'})                      # V/J genes, CDR1-3
tcrkit.anarci(seq, allow={'A', 'B'}, collapse=True)       # FR1..FR4 / CDR segments
tcrkit.anarci(seq, allow={'A', 'B'}, collapse='gene_segment', domains=True)  # leader/V/N-D/J/C
tcrkit.anarci(seq, allow={'A', 'B'}, domains=True, show=True)   # prints seq with V/J/C, CDR lines; returns None
tcrkit.anarci(seq, allow={'A', 'B'}, domains=True, mark_regions=True)  # region_seq, segment_seq
tcrkit.build_tcr_pmhc(df)                                 # all steps, self-checked
```

- `stitch` needs `stitchrdl -s human` once; it repairs nothing - run
  `add_cdr3_junctions` first. `chain=` is not inferred (default beta). Output is
  `sequence_aa`; `include_nt=True` adds `sequence_nt` (also on `stitch_pair`,
  `build_tcr_pmhc` as `{chain}_sequence_nt`).
- `anarci`: if looking for alpha/beta chains pass `allow={'A','B'}`, else an alpha
  chain can come back as delta. `CDR3_aa` is IMGT 105-117 (no C/F).
  `domains=True`: per residue `domain` V/C, `gene_segment` L(eader)/V/J/N-D/C; per sequence
  `domains` / `domain_spans`. The C domain is IMGT C-DOMAIN + 1000 (`'1023'`); V/J
  split in the CDR3 = longest germline match (Stitchr's full J if installed), can be a
  residue off. Genes next to per-residue rows: `res.join(anarci(seq)[['v_gene']])`.
- `normalize_tcr`: a valid but wrong `species` silently changes the gene; leave it unset.
- `build_tcr_pmhc` expects `alpha_v/alpha_j/alpha_cdr3`, `beta_*`, `mhc_a` (`mhc_b`).
- Full reference: `tcrkit docs api` (prints `docs/api.md`).
