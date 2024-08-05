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
 
  psi_tilde <- (-d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)] + diag(B), nrow=d, ncol=d))) +
    (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)] + diag(B))^(-1), nrow=d,ncol=d)%*%dif
  
  #if t+1 = Time, we set psi_pa to be NA to make psi_tilde = 0
  if(is.na(psi_tilde)){
    psi_tilde <- 0
  }
  
  return(psi_tilde)
}
