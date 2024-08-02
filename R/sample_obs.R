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
sample_obs <- function(model, Time, dist = 'lg'){
  if(dist = 'lg'){
    A <- model$A
    B <- model$B
    den_mean <- model$den_mean
    den_cov <- model$den_cov
    d <- model$d

    X <- matrix(0, nrow = Time, ncol = d)
    data <- matrix(0, nrow = Time, ncol = d)

    X[1,] <- rnorm(d)

    for(t in 2:Time){
      X[t,] <- rmvn(1, A%*%X[t-1,], B)
    }

    for(t in 1:Time){
      data[t,] <- rmvn(1, den_mean%*%X[t,], den_cov)
    }
  }

  return(data)
}

#' @import mvnfast
