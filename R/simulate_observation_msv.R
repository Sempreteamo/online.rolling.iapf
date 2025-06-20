#' Function to simulate observations for multivariate svm model
#'
#' @param state State which to evaluate observations at
#' @param params Parameters of the observation distribution
#'
#' @return Observations generated at the specific state
#' @export
#'
simulate_observation_msv <- function(x, obs_params){
  d <- length(x)
  const <- obs_params[[1]]
  cov <- obs_params[[2]]

  c_val <- rep(const, d)

  Sigma_xi <- (pi^2 / 2) * cov

  mu_t <- c_val + x


  # Compute multivariate log-likelihood
  obs <- mvtnorm::rmvnorm(1, mean = mu_t, sigma = Sigma_xi)

  return(obs)
}
#' @import stats
