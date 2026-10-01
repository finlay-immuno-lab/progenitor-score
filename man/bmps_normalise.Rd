\name{bmps_normalise}
\alias{bmps_normalise}
\title{Raw counts to log1p(CP10k)}
\description{Raw counts to log1p(CP10k).}
\usage{bmps_normalise(counts, target_sum = 1e4)}
\arguments{
  \item{counts}{Raw counts, cells x genes (pass only the model genes present).}
  \item{target_sum}{Library size to scale to.}
}
\value{A dense matrix of log1p counts per target_sum.}
