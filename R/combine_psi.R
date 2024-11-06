#' Function generating preferred psi function arguments
#' 
#' This function takes the parameters of the psi functions produced by each 
#' of the two sets of iAPF and combines them to give a single preferred value at each
#' time. Easily extensible to more than 2 iAPFs if needed.
#'
#' @param psi Parameters of the psi functions
#' @param index Indicator of which psi is preferred at each time 
#'
#' @return Preferred psi function arguments
#' 
#' @export
#'
combine_psi <- function(psi, index){
  dims <- dim(psi[[1]])
  Time <- length(index)
  d <- dims[2]
  
  combined_psi <- matrix(NA, Time, d)
    
   for(t in 1: Time){
     combined_psi[t,] <- psi[[index[t]]][t,]
   }
    
    return(psi = combined_psi)
}

