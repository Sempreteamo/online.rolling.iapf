#' Function to simulate observations
#'
#' @param state State which to evaluate observations at
#' @param params Parameters of the observation distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation <- function(state, params){
  den_mean <- params[[1]]
  den_cov <- params[[2]]
  obs <- rmvn(1, den_mean%*%state, den_cov)
  return(obs)
}
#' @import mvnfast
