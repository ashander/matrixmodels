\name{StabAge}
\alias{StabAge}
\title{Compute stable age distribution}
\usage{
  StabAge(A)
}
\arguments{
  \item{A}{projection matrix for population model x(t+1) =
  A x(t)}
}
\description{
  Compute stable age distribution
}
\details{
  The stable age distribution is the left eigenvector
  corresponding to the dominant eigenvalue of A.
}
\examples{
A = matrix(c(0,1,.5,0), nrow=2)
StabAge(A)
}

