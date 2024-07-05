#' Function to compute the parameters of a twisting Gaussian mixture distribution 
#'
#'This function takes the parameters of a Gaussian mixture distribution and 
#'a Gaussian twisting function and returns the parameters of the corresponding twisted distribution 
#'obtained by multiplying the two together and renormalizing.

#' @param params List of Gaussian mixture parameters
#' @param psi Parameters of the twisting function
#'
#' @return List of twisted Gaussian mixture parameters
#' @export
#'
compute_twisted_params <- function(params, psi){
  d <- length(psi)/2
  dist_mu <- params$mean
  dist_cov <- params$cov
  psi_mu <- psi[1:d]
  psi_cov <- psi[(d+1):(d+d)]
  
  zc <- (psi_cov^(-1) + diag(dist_cov)^(-1))^(-1)

  mu <- zc*(diag(dist_cov)^(-1)*dist_mu + psi_cov^(-1)*psi_mu)
  
  cov <- diag(zc, nrow = d, ncol = d)
  
  return(params = list(mu, cov))
}
