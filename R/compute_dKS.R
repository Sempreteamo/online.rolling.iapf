compute_dKS <- function(x, w, fks.obj){
  d <-  ncol(as.matrix(x[1,,]))
  dist <- vector()
  for(t in 1:Time){
    w_ <- normalise_weights_in_log_space(w[t,])[[1]]
    if(d == 1){
      weighted_mean <- sum(w_*x[t,,])
    }else{
      weighted_mean <- colSums(w_*x[t,,])
    }

    fks_mean <- fks.obj$ahatt[,t]
    fks_cov <- fks.obj$Vt[,,t]

    dist[t] <- mahalanobis(weighted_mean, fks_mean, fks_cov)
  }

  plot(dist)
  return(dist)
}
