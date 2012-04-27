\name{RepVal}
\alias{RepVal}
\title{Compute reproductive value vector}
\usage{
  RepVal(A)
}
\arguments{
  \item{A}{projection matrix for population model x(t+1) =
  A x(t)}
}
\description{
  Compute reproductive value vector
}
\details{
  The reproductive value of a stage/age class is the
  expected future contribution to population growth of an
  individual in that class. Mathematically, this
  corresponds to the right eigenvector of the dominant
  eigenvalue of A.
}
\examples{
A = matrix(c(0,1,.5,0), nrow=2)
RepVal(A)
}

