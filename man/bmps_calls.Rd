\name{bmps_calls}
\alias{bmps_calls}
\title{Best state, runner-up and margin from decision values}
\description{Best state, runner-up and margin from decision values.}
\usage{bmps_calls(D, model)}
\arguments{
  \item{D}{Decision matrix from \code{bmps_decision}.}
  \item{model}{A model from \code{bmps_read_model()}.}
}
\value{A data frame: bm_state, runner_up, top1, margin, bm_group.}
