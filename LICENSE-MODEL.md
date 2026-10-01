# Licence of the model tables: CC BY-NC 4.0

**Covered files:** `src/progenitor_score/models/**`, `inst/extdata/models/**` (identical copy for R), `docs/bm_state_groups_v1.csv`, `docs/validation/**`.

**Licence:** Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0). Legal code: https://creativecommons.org/licenses/by-nc/4.0/legalcode

You may copy, redistribute and adapt these files for non-commercial purposes, provided you give the attribution below, link the licence, and indicate changes.

## Required attribution

> Model weights derived from the BoneMarrowMap atlas of Zeng et al., trained by the Conor Finlay Lab (Trinity College Dublin) as `progenitor-score`.
> Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation. *Blood Cancer Discov* 2025;6:307-24. doi:10.1158/2643-3230.BCD-24-0342

## Why non-commercial

The weights are a statistical summary of the BoneMarrowMap annotated atlas (cells, labels and Supplementary Table S2 markers), which was assembled from six public 10x datasets (HCA bone marrow, Oetjen 2018, Granja 2019, Setty 2019, Mende 2022, Ainciburu 2022). At the time of writing:

* the Zeng et al. article is published under CC BY-NC-ND 4.0;
* the BoneMarrowMap repository (https://github.com/andygxzeng/BoneMarrowMap) is GPL >= 3 for its code and states no separate licence for the atlas objects;
* each source dataset has its own terms.

No source states that derived model weights may be used commercially, so this repository takes the more restrictive choice. A restriction can later be loosened if the authors confirm in writing that they permit it; a permissive licence cannot be tightened once the tables are public. None of the BoneMarrowMap code is used in this repository.

*This is a working choice made by the lab, not legal advice. If anyone wants commercial use, ask the BoneMarrowMap authors and the TCD technology-transfer office first.*
