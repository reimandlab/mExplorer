\name{mExplorer}
\alias{mExplorer}
\title{Selection of process-specific regulators from high-throughput data using multinomial regression models.}
\usage{
  mExplorer(dframe, response, interactions = F,
    significance = 0.05, n_cores = 1, multitest = "BY")
}
\arguments{
  \item{dframe}{Data frame of predictors. Row and column
  names are required for identifying samples (genes) and
  predictors (gene regulators), respectively.}

  \item{response}{Vector of factors. Names of vector need
  to correspond to rownames in dframe.}

  \item{interactions}{If enabled, pairs of predictors as
  interactions will be evaluated (much slower).}

  \item{significance}{Significance cutoff for p-values from
  log likelihood ratio tests.}

  \item{n_cores}{Number of processor cores to engage in
  computation. Use all available cores by default
  (n_cores=0).}

  \item{multitest}{Method to perform multiple testing
  correction for p-values from predictor evaluation. See
  p.adjust() for details.}
}
\value{
  Vector of scores, with names corresponding to predictors.
}
\description{
  Selection of process-specific regulators from
  high-throughput data using multinomial regression models.
}
\examples{
data(yeastCCgenes)
	data(yeastTFdata)
 mExplorer(yeastTFdata, yeastCCgenes)
}
\author{
  Juri Reimand <juri.reimand@utoronto.ca>
}
\references{
  TBA
}

