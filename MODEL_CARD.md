# Model card: bone_marrow_zeng2025_v1

| | |
|---|---|
| Task | assign each cell to one of 55 BoneMarrowMap states; group the state into progenitor / monocyte-precursor / mature monocyte/DC / other |
| Source data | the annotated BoneMarrowMap h5ad (263,159 cells, 45 donors, 6 datasets, 10x v2, healthy bone marrow) |
| Training set | 58,908 cells, cap 300 per state x study; 2,997 genes = the atlas's variable genes plus the top-50 markers of every state (Supp. Table S2) |
| Model | one-vs-rest logistic regression trained with stochastic gradient descent (CellTypist `use_SGD`, alpha 1e-4, 200 iterations); stored as `coef.csv.gz` (genes x states), `states.csv` (intercepts, groups), `scaler.csv` (mean, scale) |
| Input | raw counts; library size taken over the model genes present; log1p(CP10k); (x - mean)/scale; clipped above 10; decision = scaled x coef + intercept |
| Verified | `model.json` carries sha256 of each table; `load_model` refuses a changed file |
| Reproducibility | decisions are deterministic given the tables; retraining is a new model id |

**Intended use:** flagging progenitor-like cells in single-cell datasets for review. **Not for:** clinical use; calling cell types in tissues absent from bone marrow; calling progenitors without checking the margin and the neighbouring state.

**Limits:** see README Caveats. Validation: 9 held-out bone-marrow donors (`docs/validation/`). Group-level best-state recall/precision for the weak groups (granulocyte_primed, lymphoid_shared) is 0.91/0.76 and 0.70/0.63.

**Known incompatibility:** CellTypist 1.7.1's default trainer fails with scikit-learn >= 1.8 (`multi_class`), which is why the SGD route was used and why the tables, not a pickle, are the source of truth.

## Licence and attribution
Model tables: CC BY-NC 4.0 (`LICENSE-MODEL.md`); code: MIT (`LICENSE`). Cite Zeng et al. *Blood Cancer Discov* 2025;6:307-24, doi:10.1158/2643-3230.BCD-24-0342 (see `ATTRIBUTION.md`). The article itself is CC BY-NC-ND 4.0.
