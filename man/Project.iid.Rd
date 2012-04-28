\name{Project.iid}
\alias{Project.iid}
\title{Project population growth of all age classes in a stochastic environment after fixed time}
\usage{
  Project.iid(n.0, matrices, p, t.max, nreps = 100)
}
\arguments{
  \item{n.0}{initial number in stages}

  \item{matrices}{is a list of matrices or d^2 x s matrix,
  where d number of stages, s the number of matrixes}

  \item{p}{vector of probabilities for drawing each matrix}

  \item{t.max}{time to stop and compute population size}

  \item{nreps}{number of times to compute}
}
\description{
  Project population growth of all age classes in a
  stochastic environment after fixed time
}
\examples{
A1 = matrix(c(0,1,.5,0), nrow=2)
A2 = matrix(c(0,2,.5,0), nrow=2)
A3 = matrix(c(0,7,.5,0), nrow=2)
## d=2 stages, to generating appropriate d^2 by 3 matrix for input...
d = 2
mats = c(A1, A2, A3) # s = 3
dim(mats) = c(d^2, length(mats)/d^2)
Project.iid(matrix(1, d, 1), mats, rep(1/3, 3), 10)
}

