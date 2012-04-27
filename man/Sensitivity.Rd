\name{Sensitivity}
\alias{Sensitivity}
\title{Compute a matrix of sensitivities}
\usage{
  Sensitivity(A)
}
\arguments{
  \item{A}{projection matrix for population model x(t+1) =
  A x(t)}
}
\description{
  Compute a matrix of sensitivities
}
\details{
  The sensitivity is defined as the change in dominant
  eigenvalue, i.e., the growth rate, with changes in an
  entry in matrix A. See Caswell 2001.
}
\examples{
A = matrix(c(0,1,.5,0), nrow=2)
Sensitivity(A)
}

