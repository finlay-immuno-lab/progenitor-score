# Attribution

Using `progenitor-score` in work you publish or share? Please cite **both** of these.

1. **The atlas the model is trained on (required)**
   Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation. *Blood Cancer Discovery* 2025;6:307-24. doi:10.1158/2643-3230.BCD-24-0342. Atlas and R package: https://github.com/andygxzeng/BoneMarrowMap

2. **This tool**
   Conor Finlay Lab, Trinity College Dublin. progenitor-score (version in `CITATION.cff`).

The source datasets inside the BoneMarrowMap atlas are HCA bone marrow, Oetjen et al. 2018, Granja et al. 2019, Setty et al. 2019, Mende et al. 2022 and Ainciburu et al. 2022; cite them too if you discuss the training data.

| What | Licence |
|---|---|
| Code (Python, R, tests, scripts) | MIT, see `LICENSE` |
| Model tables, state-group table, validation tables | CC BY-NC 4.0 with the attribution above, see `LICENSE-MODEL.md` |

Where the attribution travels with the tool: `progenitor-score --cite`, `progenitor_score.citation()`, `bmps_citation()` in R, the `uns["progenitor_score"]` entry written by `annotate()`, and the `attribution` field of `model.json`.
