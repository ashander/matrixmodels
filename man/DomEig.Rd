\name{DomEig}
\alias{DomEig}
\title{Compute dominant eigenvalue}
\usage{
  DomEig(A)
}
\arguments{
  \item{A}{projection matrix for population model x(t+1) =
  A x(t)}
}
\description{
  Compute dominant eigenvalue
}
\details{
  The dominant eigenvalue is the eigenvalue of A with the
  largest real part.
}
\examples{
A = matrix(c(0,1,.5,0), nrow=2)
DomEig(A)
}

