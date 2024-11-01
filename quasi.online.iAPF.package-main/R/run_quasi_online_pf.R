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
  index <- data$psi_index
  obs <- data$obs
  Time <- nrow(obs)
  d = ncol(obs)
  Xs <- array(NA, dim = c(nrow(obs), N, d))
  output <- run_iAPF(model, data, Napf)
  X_record <- output[[1]]
  w_record <- output[[2]]
  psi_pa <- output[[3]]
  #logZ <- output[[4]]
  #ancestors <- output[[5]]

  psi_final1 <- combine_psi(psi_pa, index)

  output1 <- run_psi_APF(model, list(obs, breaks[[1]][1], 0, 0), Napf, psi_final1, init = FALSE)
  X <- output1[[1]]
  w <- output1[[2]]
  logZ <- output1[[3]]
  ancestors <- output1[[4]]
  resample_time <- output1[[5]]
  avg <- output1[[7]]

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
  #logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]

  return(list(X = X, w = w, logZ = logZ, Xs = Xs, avg = avg))
}



