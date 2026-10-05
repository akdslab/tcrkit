# tcrkit

Build full-length TCR sequences from V/J gene calls and a CDR3, number and annotate them
in the IMGT scheme, and standardize messy gene and allele names. Usable as a Python
library or a command-line tool.

| Capability           | What it does                                                  | Backend                                               |
| -------------------- | ------------------------------------------------------------- | ----------------------------------------------------- |
| `stitch`             | build a full TCR nt/aa sequence from V/J calls + CDR3         | [Stitchr](https://jamieheather.github.io/stitchr/)    |
| `anarci`             | annotate V/J genes and CDR1/2/3 regions of a protein sequence | [ANARCI](https://github.com/oxpig/ANARCI)             |
| `normalize_tcr`      | standardize TCR gene/allele names to IMGT                     | [tidytcells](https://tidytcells.readthedocs.io/)      |
| `normalize_mhc`      | standardize HLA/MHC allele names                              | [mhcgnomes](https://github.com/pirl-unc/mhcgnomes)    |
| `add_cdr3_junctions` | put the IMGT junction residues back on a CDR3                 | bundled reference table                               |

## How the API works

tcrkit wraps Stitchr, ANARCI, tidytcells and mhcgnomes behind one consistent interface, so
you call them the same way instead of learning four input formats.

- **One function per capability.** Each takes a single value or any collection of them — a
  list, dict, Series or DataFrame — and returns a DataFrame with your input echoed back
  beside the result. `squeeze=True` gives the bare value instead.
- **No assumed schema.** `col_mapping=` says where your columns live.
- **One `error` column** everywhere; a bad row is reported, never raised.
- **Species are `human` or `mouse`**, in and out, and spellings like `Mus musculus` are
  understood.

📓 Two runnable notebooks with outputs saved: `notebooks/tcrkit_showcase.ipynb` (a tour)
and `notebooks/tcrkit_metrics.ipynb` (the repertoire metrics).

---

## Install

```bash
pip install git+https://github.com/akdslab/tcrkit.git
```

or clone it and install in editable mode, which is what you want if you will be changing
the code:

```bash
git clone https://github.com/akdslab/tcrkit.git
cd tcrkit

uv sync                      # if you use uv
source .venv/bin/activate    # activate the env uv just created
pip install -e .             # or plain pip, into your own conda/venv
```

Prefixing commands with `uv run` works without activating, e.g. `uv run tcrkit --help`.

### Data it needs

**`stitch`** needs Stitchr's IMGT germline data, downloaded once per environment:

```bash
stitchrdl -s human
stitchrdl -s mouse          # and any other species you need
```

If `stitchrdl` fails with a download or server error, that is IMGT being unavailable
rather than anything local — retry later.

**`build_tcr_pmhc`**'s MHC step needs a table of allele sequences, which is not shipped
with the package. Pass your own with `mhc_table=` (see below). Without one it leaves the
MHC sequence columns empty and builds the TCR half anyway — it does not fail.

---

## The capabilities

### 1. `stitch` — build a full TCR sequence from V/J + CDR3

```bash
tcrkit stitch --v TRBV19 --j TRBJ2-7 --cdr3 CASSIRSSYEQYF --chain beta
```

```python
import tcrkit

tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta')
#         v        j           cdr3              sequence_aa  error
# 0  TRBV19  TRBJ2-7  CASSIRSSYEQYF  MSNQVLCCVVLCFLGANTVD...      -
```

You also get `sequence` (nucleotide) and the gene calls Stitchr actually used.

Mouse works the same way — pass `species='mouse'`, having run `stitchrdl -s mouse` once:

```bash
tcrkit stitch --v TRBV12-1 --j TRBJ2-7 --cdr3 CASSRANYEQYF --chain beta --species mouse
```

```python
# the OT-I receptor, both chains at once
tcrkit.stitch([('TRAV14-1', 'TRAJ33',  'CAASDNYQLIW'),
               ('TRBV12-1', 'TRBJ2-7', 'CASSRANYEQYF')],
              chain=['alpha', 'beta'], species='mouse')
#           v        j          cdr3 error  aa_len
# 0  TRAV14-1   TRAJ33   CAASDNYQLIW     -     275
# 1  TRBV12-1  TRBJ2-7  CASSRANYEQYF     -     313
```

Any keyword that varies per row — `chain`, `species`, `constant` — takes one value for the
whole input or one per row, and a column of the same name is used when you pass neither.
So a mixed alpha/beta frame stitches in one call.

`stitch` stitches exactly the V, J and CDR3 it is given — it repairs nothing. A CDR3
reported without its conserved C104/F118 will not stitch and says
`Unable to locate C terminus of CDR3`; complete it with
[capability 5](#5-add_cdr3_junctions--restore-cdr3-junction-residues) **first**, as its
own step, so you can see what the repair did before anything is stitched.

For large inputs, `max_workers=8` fans rows out over processes and `engine='thimble'` runs
Stitchr's CLI once for the whole input. `tcrkit.stitch_pair` does both chains of a receptor
in a single call.

#### Which constant region?

`constant=None` (the default) uses `DEFAULT_CONSTANT` — this package's fixed convention,
carried over from the data prep these wrappers came from. Stitchr's own rule is **J-gene
aware**, and the two disagree for beta:

|                             | `TRBJ1*`   | `TRBJ2*`   | alpha     |
| --------------------------- | ---------- | ---------- | --------- |
| `constant=None` (tcrkit)    | `TRBC1*01` | `TRBC1*01` | `TRAC*01` |
| `constant='auto'` (Stitchr) | `TRBC1*01` | `TRBC2*01` | `TRAC*01` |

For a TRBJ2 rearrangement that is a 5-residue, 2-length difference in the stitched chain,
so pick deliberately: `constant='auto'` for Stitchr's rule, `constant='TRBC2*01'` for an
explicit allele, or `constant={'beta': 'TRBC2*01'}` per chain.

### 2. `anarci` — annotate V/J and CDR regions of a protein sequence

```bash
tcrkit anarci --seq MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGTTAKGEVPDGYNVSRLKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL --allow A B
```

```python
tcrkit.anarci('MNAGVTQTPKFRVLKTGQSMTLLCAQDMNHDYMYWYRQDPGMGLRLIHYSVGEGT'
              'TAKGEVPDGYNVSRLKKQNFLLGLESAAPSQTSVYFCASSFTDTQYFGPGTRLTVL',
              allow={'A', 'B'})          # ← always pass allow for TCRs; see below
#   chain_type      v_gene      j_gene    CDR3_aa  CDR3_aa_start  error
# 0          B  TRBV6-2*01  TRBJ2-3*01  ASSFTDTQY             93      -
```

`CDR3_aa_start` is 1-based into the sequence you passed. `CDR3_aa` is IMGT 105–117, so it
stops *inside* the anchors and carries no C104/F118.

**Always pass** `allow=` **for TCRs.** ANARCI scores against every chain type it knows, and
the alpha and delta loci overlap, so an alpha chain comes back as `chain_type` `D` with
`j_gene` `TRDJ4*01` instead of `TRAJ29*01` — wrong, with no error. Use `{'A', 'B'}` for
alpha/beta, `{'G', 'D'}` for gamma/delta.

For the numbering itself rather than a summary, `include_imgt_index=True` adds an
`imgt_dict` column of position → residue, and `separate_index=True` spreads it into one
column per IMGT position — which lines several chains up on shared columns, since an IMGT
position means the same thing in each.

### 3. `normalize_tcr` — standardize TCR gene/allele names to IMGT

```bash
tcrkit normalize-tcr TCRBV20S1 TRAV23/6 --species human
```

```python
tcrkit.normalize_tcr(['TCRBV20S1', 'TRAV23/6'], species='human')
#       symbol standardized species error
# 0  TCRBV20S1     TRBV20-1   human  None
# 1   TRAV23/6   TRAV23/DV6   human  None
```

`tcrkit.to_gene` drops the allele (`TRBV20-1*01` → `TRBV20-1`) and `tcrkit.genes_match`
compares two calls allowing for nomenclature differences. Leave `species` unset and it is
detected, reporting which won in the `species` column. Unresolvable symbols come back
empty rather than raising.

> **A valid but wrong species changes the gene silently.** tidytcells coerces rather than
> refusing, so `normalize_tcr('TRBV20-1', species='mouse')` returns `'TRBV20'` — a
> different gene — with no error. An *unrecognised* species (`'zebra'`) does raise. Leave
> `species` unset, or compare `standardized` against `symbol`.

### 4. `normalize_mhc` — standardize HLA/MHC allele names

```bash
tcrkit normalize-mhc 'A*0201' H2-Kb
```

```python
tcrkit.normalize_mhc(['A*0201', 'H2-Kb'])
#   allele_input mhc_prefix species gene allele mhc_class error
# 0       A*0201        HLA   human    A  A0201         I  <NA>
# 1        H2-Kb         H2   mouse    K  H2-Kb         I  <NA>
```

`allele_input` is your input echoed back; everything after it is parsed out of it.
`tcrkit.mhc_chain_of('HLA-DRB1*04:01')` is a plain predicate returning `'beta'`.

Class-II pairs come back in canonical order regardless of how you wrote them, flagged in
`chain_order_swapped`. `paired=True` returns the alpha/beta layout for *every* class — a
class-I heavy chain is paired with its invariant B2M light chain — so the column shape does
not depend on which allele you passed. `errors='coerce'` records failures in `error`
instead of aborting the run.

### 5. `add_cdr3_junctions` — restore CDR3 junction residues

```bash
tcrkit add-cdr3-junctions --cdr3 ASSYSGNTEAF --j TRBJ1-1
```

```python
tcrkit.add_cdr3_junctions(['ASSYSGNTEAF', 'AVMDSNYQLI'], ['TRBJ1-1', 'TRAJ33*01'])
#           cdr3          j     cdr3_fixed  changed junction  purity       n level error
# 0  ASSYSGNTEAF    TRBJ1-1  CASSYSGNTEAFF     True      C/F     1.0  211265     J  None
# 1   AVMDSNYQLI  TRAJ33*01   CAVMDSNYQLIW     True      C/W     1.0   65517     J  None
```

The result carries the lookup it used — `junction`, `purity` (the share of reference rows
agreeing), `n` and which `level` hit — so the table's reasoning is visible, not hidden.

Only the CDR3 and the J gene are required; V and species just make the match more specific.
The J gene determines the 118 residue, which is why it is not optional — `TRBJ1-1` ends in
F, `TRAJ33` in W, `TRAJ35` in C, so this is table-driven rather than a blind
`'C' + cdr3 + 'F'`. Genes unknown to the table leave `cdr3_fixed` empty, or keep the
original with `keep_unknown=True`. `tcrkit.build_junction_table` learns a table from your
own reference data.

The table is derived from the [OTS database](https://opig.stats.ox.ac.uk/webapps/ots)
(CC-BY 4.0; per-gene aggregates only, no sequences redistributed) — please cite
[Raybould et al., *Cell Reports* 43(9):114704, 2024](https://doi.org/10.1016/j.celrep.2024.114704)
if you use it.

---

## All of it at once: `build_tcr_pmhc`

Give it a table of TCR pairs and their MHC restriction and it runs every capability in the
order these have to happen, then checks itself:

```python
out = tcrkit.build_tcr_pmhc(df)
out[out.ok]                      # the rows that came through clean
out.loc[~out.ok, 'error']        # and why the rest did not
```

| #   | step                 | does                                           |
| --- | -------------------- | ---------------------------------------------- |
| 1   | `normalize_tcr`      | V and J calls → IMGT                           |
| 2   | `add_cdr3_junctions` | CDR3s → with their 104 C / 118 F               |
| 3   | `stitch`             | V + J + CDR3 → full-length alpha and beta      |
| 4   | MHC lookup           | MHC allele(s) → both chains of the heterodimer |
| 5   | `anarci`             | the stitched sequences → V/J/CDR annotations   |

The check is step 5: ANARCI re-reads each stitched sequence knowing nothing of the input,
and its CDR3 is compared against the one that went to Stitchr — `{chain}_cdr3_ok`.

Step 4 needs a table of MHC sequences, which the package does not ship. Give it one with
`mhc_table=` — a DataFrame or CSV path with `allele` and `sequence` columns, whose allele
spellings are normalized for you, so `HLA-A*02:01` in your file is found by `A*0201` in
your data:

```python
out = tcrkit.build_tcr_pmhc(df, mhc_table='my_mhc_sequences.csv')
```

Without one, the MHC sequence columns come back empty with a warning and the rest of the
build still happens — the alleles are still normalized, just not resolved to sequences.

Input columns are looked for as `alpha_v` / `alpha_j` / `alpha_cdr3`, the `beta_`
equivalents, and `mhc_a` (plus `mhc_b` for an explicit class-II pair); `col_mapping=` says
where yours live. A chain whose columns are absent is skipped, so single-chain tables work.
Anything else you brought along — an epitope, a study id, an affinity — is carried through
untouched. `run_anarci=False` skips the slow step.

---

## Repertoire metrics

Measures over a *set* of TCRs rather than one at a time, wrapping the established package
in each case. Two need extra dependencies, so they are an optional install:

```bash
pip install 'tcrkit[metrics]'      # or: uv sync --extra metrics
```

| function                 | answers                                                             | wraps           |
| ------------------------ | ------------------------------------------------------------------- | --------------- |
| `diversity`              | how evenly is the repertoire spread? (entropy, clonality, richness) | scikit-bio      |
| `morisita_horn`          | how much do two repertoires overlap?                                | *computed here* |
| `generation_probability` | how likely was V(D)J recombination to make this CDR3?               | OLGA            |
| `tcrdist`                | how similar are two TCRs?                                           | tcrdist3        |

```python
tcrkit.diversity(df, by='sample', col_mapping={'clone': 'cdr3'})   # a row per repertoire
tcrkit.morisita_horn(df, by='sample', col_mapping={'clone': 'cdr3'})     # square matrix
tcrkit.generation_probability('CASSIRSSYEQYF', 'TRBV19', 'TRBJ2-7')     # log10_pgen
tcrkit.tcrdist(df, chains=('alpha', 'beta'))                            # square matrix
```

`by=` splits one long frame into repertoires; without it the whole input is one. Raw labels
work too — `tcrkit.diversity(['A', 'A', 'A', 'B', 'C'])` counts them for you. Both matrix
functions take `long=True` for one row per pair, and `morisita_horn(rep_a, rep_b,
squeeze=True)` is the single number for two repertoires.

Notes worth knowing:

- `shannon` is in **bits**, so `pielou_evenness` and `clonality = 1 - pielou_evenness` land
on 0–1 and are comparable between samples of different richness. Raw entropy is not — it
rises just from having more clones.
- Mind scikit-bio's naming: `simpson` **is 1 − Σp²** and `dominance` **is Σp²**. Some
immunology papers call the latter "Simpson".
- `morisita_horn` is the one metric computed here, since no maintained Python package
provides it; the usual reference is R's `vegan::vegdist(method = "horn")`, which it matches
to within 2e-15, as does every `diversity` column.
- OLGA wants the **junction** (CDR3 with its conserved C and F) — exactly what
`add_cdr3_junctions` produces, so the two compose. V and J act independently, so passing
one alone still narrows the result.
- tcrdist3 matches V/J at the **allele** level; `*01` is appended where your calls have
none. It silently drops rows whose genes it cannot resolve, so check the returned shape.

---

## CLI notes

Every subcommand takes values as flags, as shown above, *or* a table: `--in-csv` plus the
relevant `--*-col` flags, writing with `--out-csv` (or `--json`, or a readable table on
stdout). The two-step order of the Python API holds here too:

```bash
tcrkit add-cdr3-junctions --in-csv tcrs.csv --cdr3-col cdr3_aa --j-col j_call \
                          --out-col cdr3_aa --out-csv fixed.csv
tcrkit stitch --in-csv fixed.csv --v-col v_call --j-col j_call --cdr3-col cdr3_aa \
              --chain beta --out-csv stitched.csv
```

`tcrkit <command> --help` lists each command's options and examples.

## License

[MIT](LICENSE) — use, modify and redistribute it freely, including commercially, as long
as the copyright notice travels with it. The bundled CDR3 junction table is derived from
OTS and carries its own CC-BY 4.0 attribution, noted above.

Written and maintained by Ali Davari.

## Development

```bash
uv sync --extra dev --extra metrics
uv run pytest            # 114 tests, ~12 s
```

Or without `uv`: `pip install -e '.[dev]' && pytest`.

The showcase notebook doubles as an integration test — re-run it after a change and it
fails loudly if any capability regressed:

```bash
uv run jupyter nbconvert --to notebook --execute --inplace notebooks/tcrkit_showcase.ipynb
```
