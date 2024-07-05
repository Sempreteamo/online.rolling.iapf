compute_filtering <- function(model, data){
  a0 <- model$a0
  P0 <- model$P0
  dt <- model$dt
  ct <- model$ct
  Tt <- model$Tt
  Zt <- model$Zt
  Ht <- model$Ht
  Gt <- model$Gt
  obs <- data$obs
  fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))
  fks.obj <- fks(fkf.obj)
  
  return(list(fkf.obj, fks.obj))
}
