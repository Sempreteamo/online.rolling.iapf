#' Function to simulate observations for svm model
#'
#' @param state State which to evaluate observations at
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_msvm <- function(state, params){
  d <- length(state)
  obs <- rnorm(d, 0, exp(state/2))
  return(obs)
}
#' @import mvnfast
