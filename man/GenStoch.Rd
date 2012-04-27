\name{GenStoch}
\alias{GenStoch}
\title{Generate many  matrices from a stochastic forcing}
\usage{
  GenStoch(force.vec, base.matrix, element.mask, force.fxn)
}
\arguments{
  \item{force.vec}{sequence of numbers to force with}

  \item{base.matrix}{matrix upon which to perform forcing}

  \item{element.mask}{elements on which to operate can be
  matrix of logicals or ONE index i, j}

  \item{force.fxn}{forcing function: see details}
}
\description{
  Generate many matrices from a stochastic forcing
}
\details{
  The function force.fxn(forcing, elements) will take the
  elements of force.vec and elements as its arguments, and
  return a new vector representing the forced elements
}

