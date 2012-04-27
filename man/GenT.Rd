\name{GenT}
\alias{GenT}
\title{Generate a transition matrix}
\usage{
  GenT(growth.prob, surv.prob)
}
\arguments{
  \item{growth.prob}{vector of transition probabilities}

  \item{surv.prob}{vector of survival probabilites}
}
\description{
  Generate a transition matrix
}
\examples{
g.p <- c(1,0.208,0.268, 0.280)
s.p <- c(.05,.3,.716, 0.839)
GenT(g.p, s.p)
}

