
#' Function to run the whole iAPF processes and generate final results
#'
#' @param model Model information
#' @param data Observations
#' @param lag Length of time intervals
#' @param Napf Number of particles used in the smoothing process
#' @param N Number of particles used in the filtering process
#' @param filter The filter method used. For linear Gaussian models, it is set to be Kalman filter
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
run_quasi_online_pf <- function(model, data, lag, Napf, N, filter){
  breaks <- data$breaks
  index <- data$psi_index
  obs <- data$obs
  fkf.obj <- filter$fkf.obj
  fks.obj <- filter$fks.obj
  Xs <- array(NA, dim = c(nrow(obs), N, model$d))
  output <- run_iAPF(model, data, Napf)
  #X <- output[[1]]
  #w <- output[[2]]
  psi_pa <- output[[3]]
  logZ <- output[[4]]
  #ancestors <- output[[5]]

  psi_final <- combine_psi(psi_pa, index)

  output1 <- run_psi_APF(model, list(obs, breaks[[1]][1], 0, 0), Napf, psi_final, init = FALSE)
  X <- output1[[1]]
  w <- output1[[2]]
  logZ <- output1[[3]]
  ancestors <- output1[[4]]
  resample_time <- output1[[5]]

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

  log_ratio <- compute_log_ratio(logZ, fkf.obj)

  dist <- compute_dKS(X, w, fks.obj)

  return(list(X = X, w = w, logZ = logZ, Xs = Xs, log_ratio = log_ratio, dist = dist))
}



