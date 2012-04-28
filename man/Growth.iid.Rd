\name{Growth.iid}
\alias{Growth.iid}
\title{Compute the stochastic growth rate in an iid stochastic environment}
\usage{
  Growth.iid(matrices, p, t.max)
}
\arguments{
  \item{matrices}{list of matrices or a d^2 x s matrix,
  where d number of stages, s the number of matrixes}

  \item{p}{vector of probabilities for drawing each matrix}

  \item{t.max}{timeframe over which to compute growth rate}
}
\description{
  Compute the stochastic growth rate in an iid stochastic
  environment
}
\examples{
A1 = matrix(c(0,1,.5,0), nrow=2)
A2 = matrix(c(0,2,.5,0), nrow=2)
A3 = matrix(c(0,7,.5,0), nrow=2)
## d=2 stages, to generating appropriate d^2 by 3 matrix for input...
d = 2
mats = c(A1, A2, A3) # s = 3
dim(mats) = c(d^2, length(mats)/d^2)
Growth.iid(mats, rep(1/3, 3), 10)
}

