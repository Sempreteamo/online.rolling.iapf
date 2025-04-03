#' Function to learn twisting log_psi function parameters
#'
#'This function takes a collection of particle locations and calculates twisting log_psi function arguments.
#'
#' @param x A collection of particle locations
#' @param obs Observations
#' @param model List containing model parameters
#' @param log_likelihoods log_likelihoods of the particles generated during the iAPF
#'
#' @return Twisting log_psi function parameters
#' @export
#'
#'
#'
learn_psi <- function(s, psi_t, H_prev, model){
  N <- nrow(H_prev$X)
  d <- ncol(H_prev$X)
  X_prev <- H_prev$X
  log_likelihoods <- H_prev$log_li
  
  log_psi <- numeric(N)
  psi_pa <- numeric(2*d)
  
  if(all(is.na(psi_t))){
    
    log_psi <- log_likelihoods
    
  }else{
    for(i in 1:N){
      log_psi[i] <- log_likelihoods[i] + evaluate_psi_tilde(X_prev[i,], psi_t, model)
      
    }
  }
  
  psi_pa <- optimize_psi(X_prev, log_psi)
  
  return(psi_pa = psi_pa)
  
}


learn_psi <- function(x, model, log_likelihoods){
  output <- dim(x)
  Time <- output[1]
  N <- output[2]
  d <- output[3]
  log_psi <- matrix(NA, nrow = Time, ncol = N)
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)

  #calculate log_psi
  for(t in Time:1){

    if(t == Time){
      log_psi[t,] <- log_likelihoods[t,]

    }else{

      for(i in 1:N){
        log_psi[t,i] <- log_likelihoods[t,i] + evaluate_psi_tilde(x[t,i,], psi_pa[t+1, ], model)

      }
    }

    psi_pa[t,] <- optimize_psi(x[t,,], log_psi[t,])



  }
 
  return(psi_pa = psi_pa)

}
