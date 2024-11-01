#' Function to generate the ESS
#'
#' @param logW_un Unnomalised log weights
#'
#' @return Effective sample size
#' @export
#'
compute_ESS_log <- function(logW_un){
    logW <- normalise_weights_in_log_space(logW_un)[[1]]
    ess <- 1/sum(logW^2)
    return(ess)
}
