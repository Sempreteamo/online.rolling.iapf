#' Function to compute filtering results for Gaussian linear model
#'
#' @param params Information of the model
#' @param data Observations
#'
#' @return A list contains:
#' fkf.obj is the filtering information
#' fks.obj is the smoothing information
#'
#' @export
#'
compute_fkf_filtering <- function(params, obs){
  a0 <- params$a0
  P0 <- params$P0
  dt <- params$dt
  ct <- params$ct
  Tt <- params$Tt
  Zt <- params$Zt
  Ht <- params$Ht
  Gt <- params$Gt
  fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))
  fks.obj <- fks(fkf.obj)

  return(list(fkf.obj, fks.obj))
}

#' @import FKF
