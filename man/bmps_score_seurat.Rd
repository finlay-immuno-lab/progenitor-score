\name{bmps_score_seurat}
\alias{bmps_score_seurat}
\title{Score a Seurat object and add the results to its metadata}
\description{Score a Seurat object and add the results to its metadata.}
\usage{bmps_score_seurat(obj, assay = NULL, layer = "counts", symbols = NULL, qc_col = NULL, prefix = "bmps_", add = TRUE, ...)}
\arguments{
  \item{obj}{A Seurat object (v4 or v5). Join split layers first.}
  \item{assay}{Assay name (default: the default assay).}
  \item{layer}{Counts layer (or slot in Seurat v4).}
  \item{symbols}{Gene symbols if the row names are not symbols.}
  \item{qc_col}{Name of a logical metadata column; FALSE cells are not scored.}
  \item{prefix}{Prefix for the added metadata columns.}
  \item{add}{If FALSE, return the data frame instead of the object.}
  \item{...}{Passed to \code{bmps_score}.}
}
\value{The Seurat object with \code{bmps_*} metadata columns, or a data frame when \code{add = FALSE}.}
\details{Requires the SeuratObject package.}
