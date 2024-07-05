#' Function to optimize the Gaussian twisting function parameters
#'
#'This function takes a collection of particle locations and the associated value of the log target 
#'function evaluated at those points and fits a Gaussian twisting function via the controlled SMC approach.
#'
#' @param x Locations
#' @param lfn Log function value
#'
#' @return Twisting function parameters
#' @export
#'
optimize_psi <- function(x, lfn){
  params <- vector()
  d <- dim(x)[2]
  
  coef <- -lm(lfn~., data.frame(cbind(x^2, x)))$coefficients
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  params <- c(b/(-2*a), 1/(2*a))
  
  return(params)
}
