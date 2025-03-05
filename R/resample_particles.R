#' Function to resample
#'
#' @param log_weights log weights
#'
#' @return A list containing:
#' 
#'
#' @export
#'
resample_particles <- function(log_weights) {
  weights <- exp(log_weights - log_sum_exp(log_weights))  # Convert to normalized probabilities
  return(sample(1:length(weights), size = length(weights), replace = TRUE, prob = weights))
}
