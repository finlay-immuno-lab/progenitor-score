\name{bmps_decision}
\alias{bmps_decision}
\title{Decision values of the linear model}
\description{Decision values of the linear model.}
\usage{bmps_decision(X, genes, model)}
\arguments{
  \item{X}{log1p(CP10k) values, cells x genes.}
  \item{genes}{Gene symbols for the columns of X.}
  \item{model}{A model from \code{bmps_read_model()}.}
}
\value{A cells x states matrix of decision values. Model genes absent from X contribute 0; the scaled value is clipped above 10, as in CellTypist.}
