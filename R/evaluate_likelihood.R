#' Function to evaluate the log-likelihood of the observations
#'
#'This function evaluates the log-likelihood of the observations for Gaussian linear models.
#'It is expected to be updated to evaluate log-density of general potential functions.
#'
#' @param model List containing model parameters
#' @param x State at which to evaluate
#' @param datum Data point at which to evaluate
#' @param dist Distribution of the density function
#'
#' @return Log-density of the observation density
#' @export
#'
evaluate_likelihood <- function(model, x, datum, dist = 'lg') {
  if(dist == 'lg'){
    obs_mean <- model$obs_mean
    obs_cov <- model$obs_cov
    d <- length(x)

    dif <- as.vector(datum) - obs_mean %*% x

    likelihood <- (-d / 2) * log(2 * pi) - (1 / 2) * log(prod(diag(obs_cov))) -
      (1 / 2) * t(dif) %*% solve(obs_cov) %*% dif
  }

  return(likelihood)
}
