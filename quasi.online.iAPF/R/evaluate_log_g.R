#' Function to evaluate the log-density of the observation density
#'
#'This function evaluates the log-density of the observation density. 
#'It is expected to be updated to evaluate log-density of general potential functions.
#'
#' @param model List containing model parameters
#' @param x State at which to evaluate
#' @param datum Data point at which to evaluate
#'
#' @return Log-density of the observation density
#' @export
#'
evaluate_log_g <- function(model, x, datum){  
  C <- model$C
  D <- model$D
  d <- length(x)
  
  dif <- datum - C%*%x
  
  log_density <- (-d/2)*log(2*pi) - (1/2)*log(prod(diag(D))) + (-1/2)*t(dif)%*%D%*%dif
  return(log_density) 
}

