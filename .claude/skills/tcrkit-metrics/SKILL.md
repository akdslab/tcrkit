---
name: tcrkit-metrics
description: Using tcrkit's repertoire metrics - diversity/clonality, Morisita-Horn overlap, OLGA generation probability, tcrdist. Use when computing repertoire-level statistics with tcrkit.
---

# tcrkit repertoire metrics

Needs the extra: `pip install 'tcrkit[metrics]'` (or `uv sync --extra metrics`).

```python
tcrkit.diversity(df, by='sample', col_mapping={'clone': 'cdr3'})     # one row per repertoire
tcrkit.morisita_horn(df, by='sample', col_mapping={'clone': 'cdr3'}) # square matrix
tcrkit.morisita_horn(rep_a, rep_b, squeeze=True)                     # one number
tcrkit.generation_probability('CASSIRSSYEQYF', 'TRBV19', 'TRBJ2-7')  # log10_pgen
tcrkit.tcrdist(df, chains=('alpha', 'beta'))                         # beta_cdr3/beta_v/beta_j...
```

- A `count` column is used if present; otherwise rows (cells) are counted.
- `shannon` is in bits; `clonality` = 1 - Pielou evenness; `simpson` = 1 - Σp²,
  `dominance` = Σp².
- OLGA wants the CDR3 *with* C/F (`add_cdr3_junctions` output).
- tcrdist3 silently drops rows with unknown genes - check the shape.
- CLI: `tcrkit diversity | morisita-horn | generation-probability --help`.
