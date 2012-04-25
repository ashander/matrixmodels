
#' Compute dominant eigenvalue
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The dominant eigenvalue is the eigenvalue of A with the largest real part.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' dom.eig(A)
#' @export
dom.eig<-function(A){ 
  lambda<-Re(eigen(A)$values[1]) # dominant eigenvalue
  return(lambda)}


#' Compute stable age distribution
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The stable age distribution is the left eigenvector corresponding to
#' the dominant eigenvalue of A.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' stab.age(A)
#' @export
stab.age<-function(A){
  v<-Re(eigen(A)$vector[,1])
  v<-v/sum(v)  # stable age distribution
  return(v)}

#' Compute reproductive value vector
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The reproductive value of a stage/age class is the expected 
#' future contribution to population growth of an individual in that class.
#' Mathematically, this corresponds to the right eigenvector of the dominant
#' eigenvalue of A.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' repro.val(A)
#' @export
repro.val<-function(A){
  v<-stab.age(A)
  w<-Re(eigen(t(A))$vector[,1])
  w<-w/sum(v*w)  # reproductive values
  return(w)
}

#' Compute a matrix of sensitivities
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The sensitivity is defined as the change in dominant eigenvalue,
#' i.e., the growth rate, with changes in an entry in matrix A.
#' See Caswell 2001.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' sensitivity(A)
#' @export
sensitivity <- function(A){
  v.A <- stab.age(A)
  w.A <- repro.val(A)
  return(w.A%*%t(v.A))
}

#' Compute a matrix of sensitivities
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details The sensitivity is defined as the change in dominant eigenvalue,
#' i.e., the growth rate, with changes in an entry in matrix A.
#' See Caswell 2001.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' elasticity(A)
#' @export
elasticity <- function(A){
  S.A <-  sensitivity(A)
  lam <- dom.eig(A)
  return(S.A*A/lam * sign(A))
}

#' Flip a matrix so plotting with image makes visual sense 
#' @param A projection matrix for population model x(t+1) = A x(t)
#' @details R plotting function image can plot matrices, but default layout
#' won't correspond to the way we're used to seeing matrices.
#' This function transposes then reorders the columns. 
#' See help for function image.
#' @examples
#' A = matrix(c(0,1,.5,0), nrow=2)
#' sensitivity(A)
#' image(mflip(A))
#' @export
mflip<-function(A){t(A)[,ncol(A):1]}


