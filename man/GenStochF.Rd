\name{GenStochF}
\alias{GenStochF}
\title{Generate many  matrices from a stochastic forcing of vital rate: fecundity}
\usage{
  GenStochF(force.vec, vitals, element.mask, force.fxn)
}
\arguments{
  \item{force.vec}{sequence of numbers to force with}

  \item{vitals}{list of vital rates upon which to perform
  forcing}

  \item{element.mask}{elements of fecundity}

  \item{force.fxn}{forcing function: see details}
}
\value{
  list (of length equal to force.vec) of matrices square
  matrices with dimension specified in vitals
}
\description{
  Generate many matrices from a stochastic forcing of vital
  rate: fecundity
}
\details{
  The function force.fxn(forcing, elements) will take the
  elements of force.vec and elements as its arguments, and
  return a new vector representing the forced elements
}

