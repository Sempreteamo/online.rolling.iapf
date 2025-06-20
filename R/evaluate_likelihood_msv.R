#' Function to evaluate the log-likelihood of the observations for multivariate svm models
#'
#'This function evaluates the log-likelihood of the observations for stochastic volatility obs_paramss.
#'It is expected to be updated to evaluate log-density of general potential functions.
#'
#' @param x State at which to evaluate
#' @param datum Data point at which to evaluate
#' @param obs_params List containing parameters information of observation density
#'
#' @return Log-density of the observation density
#' @export
#'
evaluate_likelihood_msv <- function(x, y_t, obs_params) {
  # x: vector of latent log-volatilities at time t (length d)
  # y_t: vector of observations at time t (length d)
  # obs_params[[1]]: correlation matrix R

  #R <- obs_params
  d <- length(x)
  const <- obs_params[[1]]
  cov <- obs_params[[2]]

  c_val <- rep(const, d)

  Sigma_xi <- (pi^2 / 2) * cov

  mu_t <- c_val + x


  # Compute multivariate log-likelihood
  log_lik <- mvtnorm::dmvnorm(y_t, mean = mu_t, sigma = Sigma_xi, log = TRUE)

  return(log_lik)
}
#' @import stats
