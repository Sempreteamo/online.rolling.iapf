#' Function to provide the log-density of a twisted potential function
#'
#'This function evaluates a twisted likelihood when provided with the likelihood function, the
#'twisting function, the observation and the locations at which to perform the evaluation.
#'
#' @param model List containing model parameters
#' @param psi_pa Parameters of psi function
#' @param x States at which to evaluate
#' @param likelihoods Likelihoods of the data point at which to evaluate
#'
#' @return Twisted potential function evaluated at x
#' @export
#'
eval_twisted_potential <- function(model, psi_pa, x, likelihoods){
  ini_mu <- model$ini_mu
  ini_cov <- as.matrix(model$ini_cov)
  d <- length(x)
  psi_d <- psi_pa[[1]]
  psi_t <- psi_pa[[2]]
  psi <- psi_pa[[3]]

  dif <- ini_mu - psi_d[1:d]

  add <- diag(psi_d[(d+1):(d+d)], nrow=d, ncol=d) + ini_cov

  inv_cov <- solve(add)

  log_det_cov <- log(det(add))
  #if the time t is 1, density should not be 0
  #if time t  is T, set the psi_pa of psi_tilde as NA
  #if t is T+1, set psi_pa of psi as NA

  psi_t_0 <- -(d/2)*log(2*pi) - (1/2)* log_det_cov -
    (1 / 2) * t(dif) %*% inv_cov %*% dif

  if(is.na(psi_t_0)){
    psi_t_0 <- 0
  }

  potential <- likelihoods + evaluate_psi_tilde(x, psi_t, model) + psi_t_0 -
    evaluate_psi(x, psi)


  return(potential)
}
