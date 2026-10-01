# Changelog

## 0.1.1 (2026-10-01)
* Fixed: `bmps_score_seurat()` could surface SeuratObject's "the `slot` argument ... is now defunct" message instead of a useful error. It now calls `GetAssayData()` correctly for SeuratObject >= 5.0.0 (`layer=`, never `slot=`) and, when a layer is missing or split across samples (unjoined `Assay5`), raises a clear error naming the available layers and pointing at `SeuratObject::JoinLayers()`. Added a regression test (`tests/R/test_e2e.R`) that merges two Seurat objects with split counts layers and checks the error message, then checks that scoring works after `JoinLayers()`.

## 0.1.0
* First tagged version. Python package `progenitor_score` and R package `progenitorscore`; model distributed as plain CSV tables (no pickle); CLI, `annotate()`/`score_counts()` (Python), `bmps_score()`/`bmps_score_seurat()` (R); simulated-data demos (notebook and Rmd).
