# progenitor_score.R: score cells against BoneMarrowMap with the exported model tables. Base R; Matrix for sparse input; Seurat optional.
# Model: CC BY-NC 4.0, trained on the BoneMarrowMap atlas (Zeng et al. Blood Cancer Discov 2025;6:307-24). Code: MIT. See bmps_citation().
#
# One call on a counts matrix (genes x cells, RAW counts, gene SYMBOLS as rownames):
#   source("R/progenitor_score.R")  # or library(progenitorscore) after installing the package
#   res <- bmps_score(counts)                       # data.frame: bm_state, runner_up, top1, margin, bm_group, progenitor_call
# On a Seurat object (v4 or v5; join split layers first):
#   obj <- bmps_score_seurat(obj)                   # adds bmps_* columns to obj@meta.data
# Low-level pieces (cells x genes): bmps_normalise, bmps_decision, bmps_calls, bmps_decide.

bmps_citation <- function() {
  cat("Please cite:\n  1. Zeng AGX, Iacobucci I, Shah S, Mitchell A, Wong G, Bansal S, Chen D, Gao Q, Kim H, Kennedy JA, Arruda A, Minden MD, Haferlach T, Mullighan CG, Dick JE. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation. Blood Cancer Discov 2025;6:307-24. doi:10.1158/2643-3230.BCD-24-0342\n  2. Conor Finlay Lab, Trinity College Dublin. progenitor-score (code MIT; model tables CC BY-NC 4.0).\n", sep = "")
  invisible(NULL)
}

bmps_default_model_dir <- function() {
  d <- system.file("extdata", "models", "bone_marrow_zeng2025_v1", package = "progenitorscore")
  if (nzchar(d)) return(d)
  for (root in c(".", "..")) {                                                                 # running from a repo checkout (root or examples/)
    d <- file.path(root, "src", "progenitor_score", "models", "bone_marrow_zeng2025_v1")
    if (dir.exists(d)) return(d)
  }
  stop("model tables not found: install the package, or run from the repository root, or pass the model directory")
}

bmps_read_model <- function(dir = bmps_default_model_dir()) {
  need <- c("coef.csv.gz", "scaler.csv", "states.csv", "model.json")
  stopifnot(all(file.exists(file.path(dir, need))))
  coef <- read.csv(gzfile(file.path(dir, "coef.csv.gz")), row.names = 1, check.names = FALSE)
  sc <- read.csv(file.path(dir, "scaler.csv")); st <- read.csv(file.path(dir, "states.csv"), check.names = FALSE)
  stopifnot(identical(rownames(coef), sc$gene), identical(colnames(coef), st$state))
  list(genes = sc$gene, states = st$state, groups = st$group, coef = t(as.matrix(coef)), intercept = st$intercept, mean = sc$mean, scale = sc$scale)
}

bmps_normalise <- function(counts, target_sum = 1e4) {          # cells x genes
  counts <- as.matrix(counts); tot <- rowSums(counts); f <- ifelse(tot > 0, target_sum / tot, 0)
  log1p(counts * f)
}

bmps_decision <- function(X, genes, model) {                    # X: cells x genes, log1p(CP10k); returns cells x states
  X <- as.matrix(X); idx <- match(genes, model$genes); ok <- !is.na(idx)
  if (!any(ok)) stop("no query genes overlap the model; gene symbols are required")
  X <- X[, ok, drop = FALSE]; idx <- idx[ok]
  S <- sweep(sweep(X, 2, model$mean[idx], "-"), 2, model$scale[idx], "/"); S[S > 10] <- 10   # CellTypist clips above 10 only
  D <- S %*% t(model$coef[, idx, drop = FALSE]); sweep(D, 2, model$intercept, "+")
}

bmps_calls <- function(D, model) {                              # best state, runner-up margin, group
  o <- t(apply(D, 1, function(r) order(r, decreasing = TRUE)[1:2])); r <- seq_len(nrow(D))
  data.frame(bm_state = model$states[o[, 1]], runner_up = model$states[o[, 2]], top1 = D[cbind(r, o[, 1])], margin = D[cbind(r, o[, 1])] - D[cbind(r, o[, 2])],
             bm_group = model$groups[o[, 1]], row.names = rownames(D), stringsAsFactors = FALSE)
}

bmps_decide <- function(calls, min_margin = 100, min_top1 = 0) {
  clear <- calls$margin >= min_margin & calls$top1 >= min_top1; g <- calls$bm_group
  inp <- g %in% c("myeloid_progenitor", "dc_precursor"); rev <- g %in% c("granulocyte_primed", "multipotent", "lymphoid_shared")
  ifelse(inp & clear, "progenitor_in", ifelse(rev & clear, "progenitor_review", ifelse((inp | rev) & !clear, "uncertain",
         ifelse(g == "monocyte_precursor", "monocyte_precursor", ifelse(g == "mature_mnp", "mature_mnp", "not_progenitor")))))
}

# ---- high-level entry points -------------------------------------------------------------------------------------------------------------
bmps_score <- function(counts, symbols = rownames(counts), qc = NULL, model = NULL, min_margin = 100, min_top1 = 0, chunk = 20000, quiet = FALSE) {
  # counts: genes x cells RAW counts (matrix or Matrix sparse). symbols: gene symbols for the rows. qc: logical per cell (FALSE -> "not_scored").
  if (is.null(model)) model <- bmps_read_model()
  if (is.null(symbols)) stop("gene symbols are required: pass `symbols` or give the matrix rownames")
  symbols <- as.character(symbols); stopifnot(length(symbols) == nrow(counts))
  dup <- symbols %in% symbols[duplicated(symbols)]; feat <- model$genes[model$genes %in% symbols[!dup]]
  if (length(feat) < 0.5 * length(model$genes)) stop(sprintf("only %d of %d model genes present: are gene symbols used?", length(feat), length(model$genes)))
  pos <- match(feat, symbols); n <- ncol(counts); cells <- colnames(counts); if (is.null(cells)) cells <- as.character(seq_len(n))
  keep <- if (is.null(qc)) rep(TRUE, n) else as.logical(qc); stopifnot(length(keep) == n)
  out <- data.frame(bm_state = rep(NA_character_, n), runner_up = rep(NA_character_, n), top1 = rep(NA_real_, n), margin = rep(NA_real_, n),
                    bm_group = rep(NA_character_, n), progenitor_call = rep("not_scored", n), stringsAsFactors = FALSE); rownames(out) <- cells
  idx <- which(keep)
  for (a in seq(1, length(idx), by = chunk)[seq_len(ceiling(length(idx) / chunk))]) {
    ii <- idx[a:min(a + chunk - 1, length(idx))]
    L <- bmps_normalise(t(as.matrix(counts[pos, ii, drop = FALSE])))         # library size over the model genes present, as in training
    cl <- bmps_calls(bmps_decision(L, feat, model), model); cl$progenitor_call <- bmps_decide(cl, min_margin, min_top1)
    out[ii, c("bm_state", "top1", "margin", "bm_group", "progenitor_call")] <- cl[, c("bm_state", "top1", "margin", "bm_group", "progenitor_call")]
    out$runner_up[ii] <- cl$runner_up
  }
  if (!quiet) message("progenitor-score: model trained on BoneMarrowMap (Zeng et al. 2025, doi:10.1158/2643-3230.BCD-24-0342); model tables CC BY-NC 4.0. bmps_citation() for references.")
  out
}

bmps_score_seurat <- function(obj, assay = NULL, layer = "counts", symbols = NULL, qc_col = NULL, prefix = "bmps_", add = TRUE, ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required for bmps_score_seurat")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(obj)
  # SeuratObject >= 5.0.0: GetAssayData() takes `layer`; its `slot` argument is deprecated there and defunct (errors) as of 5.4.0 , so it
  # must never be passed on that version. Older SeuratObject (< 5.0.0) has no `layer` argument and needs `slot` instead.
  old_api <- tryCatch(utils::packageVersion("SeuratObject") < "5.0.0", error = function(e) FALSE)
  cts <- tryCatch(if (old_api) SeuratObject::GetAssayData(obj, assay = assay, slot = layer) else SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
                   error = function(e) NULL)
  if (is.null(cts) || ncol(cts) != ncol(obj)) {
    avail <- tryCatch(paste(SeuratObject::Layers(obj, assay = assay), collapse = ", "), error = function(e) NA)
    stop(sprintf("layer '%s' of assay '%s' is missing or incomplete (available layers: %s); if counts are split across samples, run SeuratObject::JoinLayers() first",
                 layer, assay, if (is.na(avail)) "unknown" else avail))
  }
  qc <- if (is.null(qc_col)) NULL else obj[[qc_col, drop = TRUE]]
  res <- bmps_score(cts, symbols = if (is.null(symbols)) rownames(cts) else symbols, qc = qc, ...)
  if (!add) return(res)
  names(res) <- paste0(prefix, names(res)); SeuratObject::AddMetaData(obj, res)
}
