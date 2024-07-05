#' Function to sample from a twisted transition distribution
#'
#'This function samples from a twisted transition distribution specified as a Gaussian mixture model.

#' @param x Initial states that one-step-ahead of the states we want to sample
#' @param model List containing model parameters
#' @param psi Parameters of the twisting function
#' @param N
#'
#' @return The sample values 
#' @export
#'
sample_twisted_transition <- function(x, model, psi, N){
  A <- model$A
  B <- model$B
  
  params <- list(mean = A%*%x, cov = B)
  
  output <- compute_twisted_params(params, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]

  samples <- rmvn(N, mu, cov)
  
  return(samples)
}

#' @import FKF
