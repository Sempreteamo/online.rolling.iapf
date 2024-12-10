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
learn_psi <- function(x, obs, model, log_likelihoods){
  output <- dim(x)
  Time <- output[1]
  N <- output[2]
  d <- output[3]
  obs <- as.matrix(obs)
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
