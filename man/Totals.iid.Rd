\name{Totals.iid}
\alias{Totals.iid}
\title{Compute distribution of total populations after fixed time}
\usage{
  Totals.iid(n.0, matrices, p, t.max, N)
}
\arguments{
  \item{n.0}{initial number in stages}

  \item{matrices}{is a d^2 x s matrix, where d number of
  stages, s the number of matrixes}

  \item{p}{vector of probabilities for drawing each matrix}

  \item{t.max}{time to stop and compute population size}

  \item{N}{number of replicates to compute the distribution
  over}
}
\description{
  Compute distribution of total populations after fixed
  time
}
\details{
  Should include a weighting/masking vector for totals in
  one class only
}

