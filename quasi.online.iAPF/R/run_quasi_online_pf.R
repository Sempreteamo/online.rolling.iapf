
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
  
  for(i in 1:N){
    Xs[i, Time,] <- X[i, Time,]
    a <- ancestors[Time, i] 
    for (t in (Time-1):1) {
      Xs[i, t, ] <- X[i, a, ]
      a <- ancestors[t, a]
    }
  }
  
  # Call the function to find rows where elements differ from the previous row
  different_rows <- find_different_rows(ancestors)
  logZ = 0
  for(t in c(different_rows)){
    logZ <- logZ + normalise_weights_in_log_space(w[t-1,])[[2]] 
  }
  logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]] 
  
  log_ratio <- compute_log_ratio(logZ, fkf.obj)
  
  KS_dist <- compute_dKS(X, w, fks.obj)  
  
  return(list(X = X, w = w, logZ = logZ, Xs = Xs, log_ratio = log_ratio, KS_dist = KS_dist))
}



