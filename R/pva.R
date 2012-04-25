#' Compute growth in all age classes in a stochastic environment after fixed time
#' @param n.0 initial number in stages
#' @param matrices is a d^2 x s matrix, where d number of stages, s the number of matrixes
#' @param p vector of probabilities for drawing each matrix
#' @param t.max time to stop and compute population size
#' @param qe quasi-extinction threshold 
#' @export
Grow.iid <- function(n.0, matrices, p, t.max, qe=0){
  mat.dim <- dim(matrices)
  s <- mat.dim[2]
  d <- sqrt(mat.dim[1])
  if (length(p) != s)
    stop("need value p for each matrix")
  if (sum(p) != 1)
    stop("p must sum to 1")
  # generate t.max random matrix indices
  A.indices <- sample(1:s, size=t.max, prob=p, replace=TRUE)
  n <- n.0
  for(t in seq_along(A.indices)){
    i <- A.indices[t] # extract appropriate index
    A <- matrices[,i]
    dim(A) <- c(d,d) #reshape
    n <- A%*%n # project population 1 period ahead
    if(sum(n) < qe)
      return(t) #return the time of extinction
  }
  if (qe == 0) 
    return(n) ##
  if (qe > 0)
    return(t.max+1)
}

#' Compute total popuation in stochastic environment after fixed time
#' @param n.0 initial number in stages
#' @param matrices is a d^2 x s matrix, where d number of stages, s the number of matrixes
#' @param p vector of probabilities for drawing each matrix
#' @param t.max time to stop and compute population size
#' @export
Total.iid <- function(n.0, matrices, p, t.max){
  return(sum(Grow.iid(n.0, matrices, p, t.max)))
}

#' Compute distribution of total populations after fixed time
#' @param n.0 initial number in stages
#' @param matrices is a d^2 x s matrix, where d number of stages, s the number of matrixes
#' @param p vector of probabilities for drawing each matrix
#' @param t.max time to stop and compute population size
#' @param N number of replicates to compute the distribution over 
#' @export
Rep.total.iid <- function(n.0, matrices, p, t.max, N){
  out<- numeric(N)
  out <- sapply(out, function(z){ z=0;  sum(Grow.iid(n.0, matrices, p, t.max))})
  return(out)
}

#' Compute distribution of times to quasi-extinction
#' @param n.0 initial number in stages
#' @param matrices is a d^2 x s matrix, where d number of stages, s the number of matrixes
#' @param p vector of probabilities for drawing each matrix
#' @param t.max time to stop and compute population size
#' @param N number of replicates to compute the distribution over
#' @param qe quasi-extinction threshold 
#' @export
Rep.exit.time.iid <- function(n.0, matrices, p, t.max, N, qe=0){
  if (qe == 0)
    stop("quasi extinction threshold must be greater than 0")
  if (qe >= sum(n.0))
    stop("quasi extinction threshold must be below starting population sum(n.0)")
  out <- numeric(N)
  out <- sapply(out, function(z){ z=0; Grow.iid(n.0, matrices, p, t.max, qe)})
  return(out)
}


