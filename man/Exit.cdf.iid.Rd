\name{Exit.cdf.iid}
\alias{Exit.cdf.iid}
\title{Compute cdf for probability of extinction versus time}
\usage{
  Exit.cdf.iid(n.0, matrices, p, t.max, N, qe = 0)
}
\arguments{
  \item{n.0}{initial number in stages}

  \item{matrices}{is a d^2 x s matrix, where d number of
  stages, s the number of matrixes}

  \item{p}{vector of probabilities for drawing each matrix}

  \item{t.max}{time to stop and compute population size}

  \item{N}{number of replicates to compute the distribution
  over}

  \item{qe}{quasi-extinction threshold}
}
\description{
  Compute cdf for probability of extinction versus time
}
\details{
  wrapper for ecdf function returns a function
}

