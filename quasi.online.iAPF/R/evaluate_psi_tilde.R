#' Function to evaluate a twisting function psi_tilde
#'
#'This function evaluates the one-step-ahead expectation of a twisting function with
#'the specified parameters at one or more points.
#'
#' @param x Points to evalute the function at
#' @param psi_pa Parameters of twisting function
#' @param model List containing model parameters
#'
#' @return Twisting function psi-tilde evaluated at specified points
#' @export
#'
evaluate_psi_tilde <- function(x, psi_pa, model){
  d <- length(x)
  A <- model$tran_mu
  B <- model$tran_cov

  dif <- as.vector(A%*%x - psi_pa[1:d])

  full_covariance <- diag(psi_pa[(d + 1):(d + d)], nrow=d, ncol=d) + B

  psi_tilde <- tryCatch({
    (-d / 2) * log(2 * pi) - (1 / 2) * log(det(full_covariance)) -
      (1 / 2) * t(dif) %*% solve(full_covariance)%*% dif
  }, error = function(e) {

    return(0)
  })

  if(is.na(psi_tilde) == TRUE){
    psi_tilde <- 0
  }
  return(psi_tilde)
}
