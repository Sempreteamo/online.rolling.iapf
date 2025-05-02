#' Rolling window function
#'
#'
#'
#' @param Z normalizing constant estimate
#' @param l number of iteration
#' @param k user-specified constant
#'
#' @return Number of particles next iteration needed
#' @export
Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

