# progenitor-score

Call hematopoietic progenitor types in single-cell data by scoring each cell against the **BoneMarrowMap** atlas (Zeng et al., *Blood Cancer Discovery* 2025;6:307-24; healthy adult bone marrow, 55 states). Works on any dataset with raw counts and human gene symbols.

* A linear model (one weight per gene per state) stored as **plain tables**. No pickle is read. Python needs only numpy/pandas/scipy/anndata; **base R** scores with the same tables; a CellTypist model can be rebuilt from them.
* Per cell: best-matching state, runner-up, `top1` and `margin` (decision values), a group, and a call.
* Status: **0.1.1, validated on held-out bone-marrow donors only** (see Caveats). **Please cite the atlas paper** (see Licence and attribution).

## What goes in, what comes out

| Input | Requirement |
|---|---|
| Counts | **raw** UMI counts (not normalised, not log); cells x genes in Python, genes x cells in R |
| Genes | **HGNC symbols** (rows of the matrix, or a `var` column via `symbol_col`). Ensembl ids stop with an error |
| QC flag | optional boolean per cell; `FALSE` cells are returned as `not_scored` |

Output columns: `bm_state`, `runner_up`, `top1`, `margin`, `bm_group`, `progenitor_call`.

## Install

```bash
pip install "git+https://github.com/finlay-immuno-lab/progenitor-score.git"        # Python >= 3.10; or `pip install -e .` in a checkout
```
```r
remotes::install_github("finlay-immuno-lab/progenitor-score")                        # R package `progenitorscore`; or source("R/progenitor_score.R") in a checkout
```

## How it is called

| You have | Call |
|---|---|
| an `.h5ad` on disk | `progenitor-score data.h5ad --symbol-col gene_symbol --qc-col qc_pass -o scores.csv.gz` (one CSV row per cell) |
| an AnnData in memory | `ps.annotate(adata, symbol_col="gene_symbol", qc_col="qc_pass")` adds `adata.obs["bmps_*"]` and `adata.uns["progenitor_score"]` |
| a matrix and gene names | `ps.score_counts(X_cells_by_genes, symbols, cell_names, qc=...)` returns a DataFrame |
| a Seurat object (v4/v5) | `obj <- bmps_score_seurat(obj, qc_col = "qc_pass")` adds `bmps_*` columns to the metadata |
| an R matrix / dgCMatrix | `res <- bmps_score(counts, qc = qc)` (genes x cells, symbols as row names) |

```python
import anndata as ad, progenitor_score as ps
adata = ad.read_h5ad("data.h5ad")                      # X = RAW counts
ps.annotate(adata, qc_col="qc_pass")                   # defaults: min_margin=100, min_top1=0
adata.obs.bmps_progenitor_call.value_counts()
```
```r
library(progenitorscore)
res <- bmps_score(counts)                              # or: obj <- bmps_score_seurat(obj)
table(res$progenitor_call)
```

The tool does not alter your object beyond the `bmps_*` columns, never normalises it in place, and runs in seconds (about 10 s for 130,000 cells on a laptop). Scoring is a matrix product, so it does not need a GPU or the 3.7 GB atlas.

**Worked examples on simulated data with known labels:** [`examples/progenitor_score_demo.ipynb`](examples/progenitor_score_demo.ipynb) (Python) and [`examples/progenitor_score_demo.Rmd`](examples/progenitor_score_demo.Rmd) (R, rendered to [`.md`](examples/progenitor_score_demo.md)). Both read the same files in `examples/data/` and finish by checking that R and Python give identical calls. Regenerate the data with `python examples/simulate_counts.py`. The simulation tests plumbing and behaviour, not accuracy.

## What the calls cover

`bm_state` (best of the 55 bone-marrow states), `runner_up`, `top1` and `margin` are general: they describe any cell, in any lineage. The `progenitor_call` column is built around myeloid and dendritic-cell progenitors, plus multipotent and lymphoid-shared progenitors for review; erythroid, megakaryocyte and other lineages are returned as `not_progenitor`, so for those read `bm_state` directly. The grouping lives in `src/progenitor_score/groups.py` and `docs/bm_state_groups_v1.csv` and can be changed.

## Calls

| Call | Meaning |
|---|---|
| `progenitor_in` | best state is GMP-like (Early GMP, GMP-Cycle) or a DC precursor (Pre-cDC, Pre-pDC), with a clear margin |
| `progenitor_review` | HSC/MPP, LMPP, MLP, neutrophil-primed GMP, with a clear margin |
| `monocyte_precursor` | Early/Late ProMono: part of the monocyte continuum, not a progenitor call |
| `mature_mnp` | mature monocyte, macrophage-lineage or dendritic-cell state (CD14/CD16 monocyte, cDC1, cDC2, AS-DC, pDC) |
| `uncertain` | a progenitor best state without a clear margin or good match |
| `not_progenitor` | any other state (T, B, erythroid, ...) |
| `not_scored` | failed the QC flag |

The group table (`docs/bm_state_groups_v1.csv`) is a proposal, not a ruling. **Use `margin` and `top1`, not probabilities**: this model's probabilities are saturated.

### Thresholds: what `margin` and `top1` do

A call needs `margin >= 100` (a clear lead over the runner-up) **and** `top1 >= 0` (a good match to the state). They do different jobs:

| Gate (held-out bone-marrow donors, in-progenitor groups) | recall | precision | `progenitor_in` calls on a non-marrow myeloid test set (210,839 cells, ~no progenitors) |
|---|---|---|---|
| `margin >= 100` only | 0.62 | 0.89 | 962 |
| `margin >= 100`, `top1 >= -500` | 0.58 | 0.91 | 153 |
| `margin >= 100`, `top1 >= 0` (default) | 0.34 | 0.93 | 12 |

`top1 >= 0` costs recall in marrow but is what stops macrophages and other non-marrow cells from being forced onto progenitor states. Keep the default for tissue data; for marrow-like input, `min_top1=-500` recovers much of the recall (`validation_in_progenitor_margin_top1.csv`, `validation_top1_gate_non_marrow.csv`).

## How well it works

Held-out donors (9 of 45), at most 300 cells per state. Group-level (best state) recall / precision: myeloid_progenitor 0.68 / 0.85, dc_precursor 0.77 / 0.85, monocyte_precursor 0.89 / 0.79, multipotent 0.97 / 0.83, mature MNP 0.93 / 0.92, **granulocyte_primed 0.91 / 0.76 and lymphoid_shared 0.70 / 0.63 (weakest)**. A gene score from the atlas's top-50 markers per state reaches group AUC 0.96-0.97 for the in-progenitor groups, against 0.98-0.99 for this classifier. Details: `docs/validation/`, `docs/BRIEF_bone_marrow_atlas.md`, `MODEL_CARD.md`.

## Caveats

* Closed set, adult bone marrow only: cells the atlas lacks are forced onto the nearest state (microglia read as Mature B; cycling macrophages read as cycling precursors, usually `uncertain`). `not_progenitor` is not lineage evidence outside marrow-like tissue. No yolk sac, fetal liver or other developmental source is covered.
* Validated on bone-marrow donors the model never saw, not on an independent progenitor-labelled dataset. Thresholds are starting values.
* Pro-monocyte stages cannot be checked with CD34/SPINK2: those genes are gone by then.
* The shipped demos use simulated counts built from atlas averages; real cells are noisier, so expect lower recall than the demos show.

## Verification

`pytest` (13 tests) checks, against decision values from the original CellTypist model on synthetic cells: the numpy scorer (max abs difference 3e-4 on values up to 15,000; identical best states), the same with 300 genes missing, the clip at 10, the rebuilt CellTypist model, tamper detection on the tables, and end-to-end output on a synthetic count matrix (dense, sparse, AnnData, CLI). `Rscript tests/R/test_golden.R` and `Rscript tests/R/test_e2e.R` do the same for the R scorer and for `bmps_score` / `bmps_score_seurat`, including agreement with the Python output. On 129,449 real cells of one public dataset the tables-only scorer gave the same best state and call as the CellTypist route for every cell.

## Layout

`src/progenitor_score/` Python package and model tables · `R/` and `inst/extdata/` R package (identical model copy, checked by a test) · `tests/` fixtures (synthetic) and tests · `examples/` simulated data, notebook, Rmd · `training/` scripts that built and validated the model (example paths to edit; the atlas h5ad is not in this repo) · `docs/`.

## Licence and attribution

* **Code** (Python, R, tests, scripts): MIT, see `LICENSE`.
* **Model tables, state-group table, validation tables and the example state profiles**: CC BY-NC 4.0 with attribution, see `LICENSE-MODEL.md`. They are derived from the BoneMarrowMap atlas; the Zeng et al. article is CC BY-NC-ND 4.0 and no source states that derived weights may be used commercially, so this repository takes the more restrictive choice.
* **Please cite both** (also in `ATTRIBUTION.md`, `CITATION.cff`, `progenitor-score --cite`, `bmps_citation()` and the `uns["progenitor_score"]` entry):
  1. Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation. *Blood Cancer Discov* 2025;6:307-24. doi:10.1158/2643-3230.BCD-24-0342. Atlas and R package: https://github.com/andygxzeng/BoneMarrowMap (GPL >= 3; none of its code is used here).
  2. Conor Finlay Lab, Trinity College Dublin. progenitor-score.
