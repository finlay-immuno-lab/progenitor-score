\name{bmps_read_model}
\alias{bmps_read_model}
\title{Read the model tables}
\description{Read the model tables.}
\usage{bmps_read_model(dir = bmps_default_model_dir())}
\arguments{
  \item{dir}{Directory holding coef.csv.gz, scaler.csv, states.csv and model.json.}
}
\value{A list with genes, states, groups, coef (states x genes), intercept, mean and scale.}
\examples{
m <- bmps_read_model()
length(m$genes); length(m$states)
}
