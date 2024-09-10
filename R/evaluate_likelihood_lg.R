#' Function to evaluate the log-likelihood of the observations for linear gaussians
#'
#'This function evaluates the log-likelihood of the observations for Gaussian linear obs_paramss.
#'It is expected to be updated to evaluate log-density of general potential functions.
#'
#' @param x State at which to evaluate
#' @param datum Data point at which to evaluate
#' @param obs_params List containing parameters information of observation density
#'
#' @return Log-density of the observation density
#' @export
#'
evaluate_likelihood <- function(x, datum, obs_params) {


    obs_mean <- obs_params[[1]]
    obs_cov <- obs_params[[2]]
    d <- length(x)

    dif <- as.vector(datum) - obs_mean %*% x

    likelihood <- (-d / 2) * log(2 * pi) - (1 / 2) * log(prod(diag(obs_cov(x)))) -
      (1 / 2) * t(dif) %*% solve(obs_cov(x)) %*% dif


  return(likelihood)
}

