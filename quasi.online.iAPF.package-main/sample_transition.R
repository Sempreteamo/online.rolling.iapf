sample_transition <- function(x, model, N){
  A <- model$tran_mu
  B <- model$tran_cov
  ini <- model$ini_mu
  samples <- mvnfast::rmvn(N, ini + A%*%(x - ini), B)

  return(samples)
}
