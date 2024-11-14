#' Function to run the whole iAPF processes and generate final results
#'
#' @param model Model information
#' @param data Observations
#' @param lag Length of time intervals
#' @param Napf Number of particles used in the smoothing process
#' @param N Number of particles used in the filtering process
#'
#' @return A list contains
#' X is the particles
#' w is the weights of particles
#' logZ is the estimated normalising constant
#' Xs is the smoothing particles
#' log_ratio is the log_ratio between estimates and real normalising constant
#' dist is the distance between real estimates and smoothing distribution
#'
#' @export
#'
run_quasi_online_pf <- function(model, data, lag, Napf, N){
  breaks <- data$breaks
  psi_index <- data$psi_index
  obs <- data$obs
  Time <- nrow(obs)
  Xs <- array(NA, dim = c(nrow(obs), N, d))
  d = ncol(obs)
  
  X_apf <- NULL
  w_apf <- NULL
  
  
  #step 1: run iAPF whenever the end of a block is reached
  for(index in 1:length(breaks)){
    psi_pa1 <- NULL
    
    for( b in  2:length(breaks[[index]])){
      data$run_block <- c(index, b)  #which layer, which block iAPF runs
      data$past <- list(X_apf[dim(X_apf)[1],,], w_apf[nrow(w_apf),])
      
      output <- run_iAPF(model, data, Napf)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      psi_pa <- output[[3]]
      
      psi_pa1 <- rbind(psi_pa1, psi_pa)
    }
    psi_final[[index]] <- psi_pa1
  }
  
  #combine psi
  psi_final1 <- combine_psi(psi_final, psi_index)

  #step2: advance the ψ-APF
  output1 <- run_psi_APF(model, list(obs, breaks[[1]][1], 0, 0), Napf, psi_final1, init = FALSE)
  X <- output1[[1]]
  w <- output1[[2]]
  logZ <- output1[[3]]
  ancestors <- output1[[4]]
  resample_time <- output1[[5]]
  #avg <- output1[[7]]

  for(i in 1:N){
    Xs[Time, i,] <- X[Time, i,]
    a <- ancestors[Time, i]
    for (t in (Time-1):1) {
      Xs[t, i,] <- X[t, a,]
      a <- ancestors[t, a]
    }
  }

  logZ = 0
  for(t in c(resample_time)){
    logZ <- logZ + normalise_weights_in_log_space(w[t,])[[2]]
  }

  return(list(X = X, w = w, logZ = logZ, psi_final = psi_final1))
}



