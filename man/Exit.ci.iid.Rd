\name{Exit.ci.iid}
\alias{Exit.ci.iid}
\title{Compute confidence interval for probability of extinction within specified timeframe}
\usage{
  Exit.ci.iid(time, ci.width, reps, n.0, matrices, p,
    t.max, N, qe = 0)
}
\arguments{
  \item{time}{time at which to calculate CI}

  \item{ci.width}{width (in probability units) of
  confidence interval to produce}

  \item{reps}{number of distinct cdfs to compute, don't
  make this too large}

  \item{n.0}{initial number in stages}

  \item{matrices}{is a d^2 x s matrix, where d number of
  stages, s the number of matrixes}

  \item{p}{vector of probabilities for drawing each matrix}

  \item{t.max}{time to stop and compute population size}

  \item{N}{number of replicates to compute the exit time
  distribution over within each run}

  \item{qe}{quasi-extinction threshold}
}
\description{
  Compute confidence interval for probability of extinction
  within specified timeframe
}
\details{
  This MAY NOT fully represent the variablity in the system
  this function does not use distribution of vital rates,
  but rather repeated draws of complete matrices. See
  Fieberg and Ellner (2003, Ecology) and Morris and Doak
  (2003).
}

