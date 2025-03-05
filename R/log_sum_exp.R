
#' 
#'
#' Log-Sum-Exp function to prevent underflow
#'
#' @param log_weights 
#'
#' @return A list containing:
#' 
#'
#' @export
#'
log_sum_exp <- function(log_weights) {
  max_logW <- max(log_weights)
  max_logW + log(sum(exp(log_weights - max_logW)))
}

# Main function implementing ψ-APF
