#' Function to learn twisting psi function parameters
#'
#'This function takes a collection of particle locations and calculates twisting psi function arguments.
#'
#' @param x A collection of particle locations
#' @param obs Observations
#' @param model List containing model parameters
#' @param likelihoods Likelihoods of the particles generated during the iAPF
#'
#' @return Twisting psi function parameters
#' @export
#'
learn_psi <- function(x, obs, model, likelihoods){
  output <- dim(x)
  Time <- output[1]
  N <- output[2]
  d <- output[3]
  obs <- as.matrix(obs)
  psi <- matrix(NA, nrow = Time, ncol = N)
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)

  #calculate psi
  for(t in Time:1){

    if(t == Time){
      for(i in 1:N){
        dif <- as.vector(x[t, i,] - obs[t, ,drop = FALSE])

        psi[t,i] <- (1 / ((2 * pi)^(d / 2))) *
          exp(-0.5 * t(dif) %*% dif)
      }


    }else{

      for(i in 1:N){

        psi[t,i] <- exp(likelihoods[t,i] + evaluate_psi_tilde(x[t,i,], psi_pa[t+1, ], model))

      }
    }

    psi_pa[t,] <- optimize_psi(x[t,,], log(psi[t,]))


    #print(psi_pa[t, 1:d])
    #print(obs[t, drop = FALSE])

  }
  #print(psi_pa)
  return(params = psi_pa)

}
