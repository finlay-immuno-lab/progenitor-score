# Run from the repo root:  Rscript R/test_e2e.R   (end-to-end: bmps_score / bmps_score_seurat versus the Python output on the same synthetic counts)
source("R/progenitor_score.R")
df <- read.csv(gzfile("tests/e2e_counts.csv.gz"), check.names = FALSE); cts <- as.matrix(df[, -1]); rownames(cts) <- df[[1]]   # contains one duplicated symbol (must be dropped, as in Python)
qc <- read.csv("tests/e2e_qc.csv", row.names = 1)$qc; exp <- read.csv("tests/e2e_expected.csv", row.names = 1, na.strings = c("", "NA"))
check <- function(res, label) {
  stopifnot(identical(rownames(res), rownames(exp)), identical(res$progenitor_call, exp$progenitor_call))
  ok <- !is.na(exp$bm_state); stopifnot(all(res$bm_state[ok] == exp$bm_state[ok]), all(res$runner_up[ok] == exp$runner_up[ok]),
    max(abs(res$top1[ok] - exp$top1[ok])) < 1e-2, max(abs(res$margin[ok] - exp$margin[ok])) < 1e-2, all(is.na(res$bm_state[!ok])))
  cat(sprintf("%-28s OK (%d cells, %d not_scored)\n", label, nrow(res), sum(!ok)))
}
check(bmps_score(cts, qc = qc, quiet = TRUE), "dense matrix")
if (requireNamespace("Matrix", quietly = TRUE)) check(bmps_score(Matrix::Matrix(cts, sparse = TRUE), qc = qc, quiet = TRUE), "sparse dgCMatrix")
check(bmps_score(cts, qc = qc, chunk = 7, quiet = TRUE), "chunk = 7")
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  dd <- rownames(cts) %in% rownames(cts)[duplicated(rownames(cts))]    # Seurat forbids duplicate rownames; dropping both copies = what the scorer does
  obj <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cts[!dd, ], sparse = TRUE), min.cells = 0, min.features = 0)
  obj$qc <- qc; obj <- bmps_score_seurat(obj, qc_col = "qc"); md <- obj[[]]
  stopifnot(identical(md$bmps_progenitor_call, exp$progenitor_call), all(md$bmps_bm_state[!is.na(exp$bm_state)] == exp$bm_state[!is.na(exp$bm_state)]))
  cat("Seurat object                OK (columns:", paste(grep("^bmps_", names(md), value = TRUE), collapse = ", "), ")\n")

  # regression: an un-joined Assay5 (split counts.1/counts.2 layers) must give a clear, actionable error -- never the
  # "slot argument ... defunct" message that leaked through when GetAssayData() was called with a deprecated argument
  m1 <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cts[!dd, 1:20], sparse = TRUE)); m1$qc <- qc[1:20]
  m2 <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(cts[!dd, 21:40], sparse = TRUE)); m2$qc <- qc[21:40]
  merged <- merge(m1, m2)
  err <- tryCatch({ bmps_score_seurat(merged, qc_col = "qc"); NULL }, error = function(e) conditionMessage(e))
  stopifnot(!is.null(err), grepl("JoinLayers", err), !grepl("defunct|deprecat", err, ignore.case = TRUE))
  cat("Split-layer Seurat object    OK (clear error, mentions JoinLayers, no deprecated-argument message)\n")
  merged <- SeuratObject::JoinLayers(merged)
  rj <- bmps_score_seurat(merged, qc_col = "qc")[[]]
  stopifnot(identical(unname(rj$bmps_progenitor_call), exp$progenitor_call[match(rownames(rj), rownames(exp))]))
  cat("Split-layer after JoinLayers OK (40 cells)\n")
}
bmps_citation()
