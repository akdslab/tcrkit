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
tcrkit.build_tcr_pmhc(df)                                 # all steps, self-checked
```

- `stitch` needs `stitchrdl -s human` once; it repairs nothing - run
  `add_cdr3_junctions` first. `chain=` is not inferred (default beta).
- `anarci`: if looking for alpha/beta chains pass `allow={'A','B'}`, else an alpha
  chain can come back as delta. `CDR3_aa` is IMGT 105-117 (no C/F).
- `normalize_tcr`: a valid but wrong `species` silently changes the gene; leave it unset.
- `build_tcr_pmhc` expects `alpha_v/alpha_j/alpha_cdr3`, `beta_*`, `mhc_a` (`mhc_b`).
- Full reference: `docs/api.md`.
