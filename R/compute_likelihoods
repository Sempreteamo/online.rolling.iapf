#' Function to compute likelihoods for particles
#'
#' @param particles Particles to compute likelihoods
#' @param datum Observations
#' @param model Model information
#'
#' @return likelihoods are the likelihoods of the particles
#' @export
#'
compute_likelihoods <- function(particles, datum, model) {
  likelihoods <- numeric(nrow(particles))

  for (i in 1:nrow(particles)) {
    likelihoods[i] <- evaluate_likelihood(model, particles[i, ], datum)
  }

  return(likelihoods)
}
