#' Generate a transition matrix 
#' @param growth.prob vector of transition probabilities
#' @param surv.prob vector of survival probabilites 
#' @examples
#' g.p <- c(1,0.208,0.268, 0.280)
#' s.p <- c(.05,.3,.716, 0.839)
#' GenT(g.p, s.p)
#' @export
GenT <- function(growth.prob, surv.prob){
  d <- length(growth.prob)
  if( d != length(surv.prob))
    stop('growth and survival probability vectors must be same length')
  T <- matrix(nrow=d,ncol=d,0)
  for (i in 1:d){
    T[i,i] <- surv.prob[i]*(1-growth.prob[i])  # survival and not grow 
    if(i < d){T[i+1,i] <- surv.prob[i]*growth.prob[i]}
  }
  T[d,d] <- surv.prob[d]
  return(T)
}

#' Generate a fecundity matrix
#' @param d dimension of matrix (number of stages)
#' @param fec vector of fecundites starting from largest class
#' @examples
#' f.v <- c(1, 7)
#' GenF(4, f.v)
#' @export
GenF <- function(d, fec){
  F <- matrix(nrow=d,ncol=d,0)
  stage.at.mat <- d-length(fec) + 1
  F[1, stage.at.mat:d] <- fec
  return(F)
}

#' Compute dominant eigenvalue
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The dominant eigenvalue is the eigenvalue of A with the largest real part.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' DomEig(A)
#' @export
DomEig<-function(A){ 
  lambda <- Re(eigen(A)$values[1]) # dominant eigenvalue
  return(lambda)}


#' Compute stable age distribution
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The stable age distribution is the left eigenvector corresponding to
#' the dominant eigenvalue of A.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' StabAge(A)
#' @export
StabAge<-function(A){
  v <- Re(eigen(A)$vector[,1])
  v <- v/sum(v)  # stable age distribution
  return(v)}

#' Compute reproductive value vector
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The reproductive value of a stage/age class is the expected 
#' future contribution to population growth of an individual in that class.
#' Mathematically, this corresponds to the right eigenvector of the dominant
#' eigenvalue of A.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' RepVal(A)
#' @export
RepVal<-function(A){
  v <- StabAge(A)
  w <- Re(eigen(t(A))$vector[,1])
  w <- w/sum(v*w)  # reproductive values
  return(w)
}

#' Compute a matrix of sensitivities
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The sensitivity is defined as the change in dominant eigenvalue,
#' i.e., the growth rate, with changes in an entry in matrix A.
#' See Caswell 2001.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' Sensitivity(A)
#' @export
Sensitivity <- function(A){
  v.A <- StabAge(A)
  w.A <- RepVal(A)
  return(w.A%*%t(v.A))
}

#' Compute a matrix of sensitivities
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The sensitivity is defined as the change in dominant eigenvalue,
#' i.e., the growth rate, with changes in an entry in matrix A.
#' See Caswell 2001.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' Elasticity(A)
#' @export
Elasticity <- function(A){
  s.A <-  Sensitivity(A)
  lam <- DomEig(A)
  return(s.A*A/lam * sign(A))
}

#' Flip a matrix so plotting with image makes visual sense 
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details R plotting function image can plot matrices, but default layout
#' won't correspond to the way we're used to seeing matrices.
#' This function transposes then reorders the columns. 
#' See help for function image.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' Sensitivity(A)
#' image(Flip(A))
#' @export
Flip<-function(A){t(A)[,ncol(A):1]}


