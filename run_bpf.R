run_bpf <- function(model, data, lag, Napf){
  obs <- data$obs
  Time <- nrow(obs)
  d = ncol(obs)

  Xs <- array(NA, dim = c(nrow(obs), N, d))
  X <- array(NA, dim = c(Time, Napf, d))
  w <- matrix(NA, Time, Napf)
  logZ <- 0
  ancestors <- matrix(NA, Time, N)
  ini <- model$ini_mu
  ini_c <- model$ini_cov
  tran_mu <- model$tran_mu
  tran_cov <- model$tran_c

  X[1,,] <- mvnfast::rmvn(Napf, ini, ini_c)
  for(i in 1:Napf){
    w[1,i] <- model$eval_likelihood(X[1,i,], obs[1,], model$obs_params)
  }
  logZ = 0

  logZ <- logZ + normalise_weights_in_log_space(w[1,])[[2]]
  mix <- resample(w[1,])
  X[1,,] <- X[1,mix,]
  ancestors[1,] <- seq(1:N)
  for (t in 2:Time) {

    for(i in 1:Napf){

      X[t,i,] <-  mvnfast::rmvn(1, tran_mu%*%X[t-1,i,],  tran_cov)
      w[t,i] <-  model$eval_likelihood(X[t,i,], obs[t,], model$obs_params)
    }


    logZ <- logZ + normalise_weights_in_log_space(w[t,])[[2]]
    mix <- resample(w[t,], mode = 'multi')
    X[t,,] <- X[t,mix,]
    #ancestors[t,] <- mix
    #V[t,] <- V[t-1,]
    #V[t, -unique(I[t-1,mix])] <- t


  }
  #logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]

  return(list(X = X, w = w, logZ = logZ))
}
