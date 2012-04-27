\name{GenF}
\alias{GenF}
\title{Generate a fecundity matrix}
\usage{
  GenF(d, fec)
}
\arguments{
  \item{d}{dimension of matrix (number of stages)}

  \item{fec}{vector of fecundites starting from largest
  class}
}
\description{
  Generate a fecundity matrix
}
\examples{
f.v <- c(1, 7)
GenF(4, f.v)
}

