compute_log_ratio <- function(Z, fkf.obj){
  
  fkf.obj_Z <- fkf.obj$logLik
  cat('NC = ', fkf.obj_Z, 'est = ', Z, 'log_ratio = ', exp(Z-fkf.obj_Z))
  return(exp(Z-fkf.obj_Z))
}
