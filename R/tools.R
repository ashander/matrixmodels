#' Simple confidence interval
#' @param samples the data
#' @param width the width in probability units of the desired CI
#' @export 
CI <- function(samples, width=0.95){
  a <- (1-width)/2
  return(quantile(samples, probs=c(a, 1-a)))
}
