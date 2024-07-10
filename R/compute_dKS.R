compute_dKS <- function(x, w, fks.obj){
  d <-  dim(x)[2]
  dist <- vector()
  for(t in 1:Time){
    w_ <- normalise_weights_in_log_space(w[t,])[[1]]
    weighted_mean <- colSums(w_*x[t,,])
    
    fks_mean <- fks.obj$ahatt[,t]
    fks_cov <- fks.obj$Vt[,,t]
    
    dist[t] <- mahalanobis(weighted_mean, fks_mean, fks_cov)/d
  }
  
  plot(dist)
  return(dist)
}
