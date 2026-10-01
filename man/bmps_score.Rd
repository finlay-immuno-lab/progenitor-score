\name{bmps_score}
\alias{bmps_score}
\title{Score cells against the BoneMarrowMap model and call progenitor types}
\description{Score cells against the BoneMarrowMap model and call progenitor types.}
\usage{bmps_score(counts, symbols = rownames(counts), qc = NULL, model = NULL, min_margin = 100, min_top1 = 0, chunk = 20000, quiet = FALSE)}
\arguments{
  \item{counts}{RAW counts, genes x cells (base matrix or Matrix sparse). Not normalised.}
  \item{symbols}{Gene symbols for the rows (default: row names). Duplicated symbols are dropped.}
  \item{qc}{Optional logical per cell; FALSE cells get \code{progenitor_call = "not_scored"}.}
  \item{model}{A model from \code{bmps_read_model()}; default is the bundled bone-marrow model.}
  \item{min_margin}{Minimum \code{margin} (best minus runner-up decision value) for a progenitor call.}
  \item{min_top1}{Minimum \code{top1} (best decision value) for a progenitor call.}
  \item{chunk}{Cells scored per block (limits memory).}
  \item{quiet}{Suppress the attribution message.}
}
\value{A data frame with one row per cell: \code{bm_state}, \code{runner_up}, \code{top1}, \code{margin}, \code{bm_group}, \code{progenitor_call} (progenitor_in, progenitor_review, monocyte_precursor, mature_mnp, uncertain, not_progenitor, not_scored).}
\details{The library size is taken over the model genes present, as in training. Use \code{margin} and \code{top1}, not probabilities. The model is a closed-set classifier trained on healthy adult bone marrow; see the model card.}
\examples{
m <- bmps_read_model()
set.seed(1)
x <- matrix(rpois(length(m$genes) * 20, 0.05), nrow = length(m$genes), dimnames = list(m$genes, paste0("c", 1:20)))
res <- bmps_score(x, quiet = TRUE)
table(res$progenitor_call)
}
