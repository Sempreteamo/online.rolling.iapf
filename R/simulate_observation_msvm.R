#' Function to simulate observations for multivariate svm model
#'
#' @param state State which to evaluate observations at
#' @param params Parameters of the observation distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_msvm <- function(state, params){
  d <- length(state)
  obs <- stats::rnorm(d, 0, exp(state/2))
  return(obs)
}
#' @import stats
