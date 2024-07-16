#' Function to sample observations
#'
#'This function samples synthetic observations for each time step given the state space model.
#'
#' @param model List containing model parameters
#' @param N Number of observations
#'
#' @return A list containing the observation sequence
#' @export
#'
sample_obs <- function(model){
  A <- model$A
  B <- model$B
  C <- model$C
  D <- model$D
  d <- model$d

  X <- matrix(0, nrow = Time, ncol = d)
  data <- matrix(0, nrow = Time, ncol = d)

  X[1,] <- rnorm(d)

  for(t in 2:Time){
    X[t,] <- chol(B)%*%rnorm(d) + A%*%X[t-1,]
  }
  for(t in 1:Time){
    data[t,] <- rmvn(1, C%*%X[t,], D)
  }
  return(data)
}

#' @import mvnfast
