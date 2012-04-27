\name{Flip}
\alias{Flip}
\title{Flip a matrix so plotting with image makes visual sense}
\usage{
  Flip(A)
}
\arguments{
  \item{A}{projection matrix for population model x(t+1) =
  A x(t)}
}
\description{
  Flip a matrix so plotting with image makes visual sense
}
\details{
  R plotting function image can plot matrices, but default
  layout won't correspond to the way we're used to seeing
  matrices. This function transposes then reorders the
  columns. See help for function image.
}
\examples{
A = matrix(c(0,1,.5,0), nrow=2)
Sensitivity(A)
image(Flip(A))
}

