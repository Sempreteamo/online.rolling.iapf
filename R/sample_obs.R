#' Function to sample observations
#'
#'This function samples synthetic observations for each time step given the state space model.
#'
#' @param model List containing model parameters
#' @param dist Distribution of the observation density
#' @param Time Total time length
#'
#' @return A list containing the observation sequence
#' @export
#'
sample_obs <- function(tran_mu, tran_cov, den_mean = 0, den_cov = 0, Time, d, dist = 'lg'){
  if(dist == 'lg'){
    X <- matrix(0, nrow = Time, ncol = d)
    data <- matrix(0, nrow = Time, ncol = d)

    X[1,] <- rnorm(d)

    for(t in 2:Time){
      X[t,] <- rmvn(1, tran_mu%*%X[t-1,], tran_cov)
    }

    for(t in 1:Time){
      data[t,] <- rmvn(1, den_mean%*%X[t,], den_cov)
    }
  }else{
    print('provide the observation data directly')
  }

  return(data)
}

#' @import mvnfast
