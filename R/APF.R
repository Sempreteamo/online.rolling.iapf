#'
#'
#' Log-Sum-Exp function to prevent underflow
#'
#' @param psi_pa psi parameters
#' @param data data that include observations
#' @param N number of particle needed
#' @param model model used
#'
#' @return A list containing:
#'
#'
#' @export
#'
APF <- function(psi_pa, data, N, model){
  #l >= 2
  ini_cov <- model$ini_cov
  ini_mu <- model$ini_mu
  kappa <- model$parameters$kappa
  obs <- data$obs
  d = ncol(obs)
  Time <- nrow(obs)
  log_likelihoods <- matrix(NA, Time, N)
  ancestors <- matrix(NA, Time, N)
  obs_params <- model$obs_params

  X <- array(NA, dim = c(Time, N, d))
  w <- matrix(NA, Time, N)
  logZ <- 0


  X[1, ,] <- sample_twisted_initial(list(mean = ini_mu,
                                         cov = as.matrix(ini_cov)), psi_pa[1,], N)  #particles

  if(Time == 1){
    for(i in 1:N){
      log_likelihoods[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
      w[1,i] <- eval_twisted_potential(model,
                                       list(NA, NA, psi_pa[1,]), X[1,i,],  log_likelihoods[1,i])
      #w[1,i] <- log(g_aux(obs[1,], X[1,i,],1, psi_pa))

    }

  }else{
    for(i in 1:N){
      log_likelihoods[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
      w[1,i] <- eval_twisted_potential(model,
                                       list(psi_pa[1,], psi_pa[2,], psi_pa[1,]), X[1,i,],  log_likelihoods[1,i])
      #w[1,i] <- log(g_aux(obs[1,], X[1,i,],1, psi_pa))

    }
  }



  #re=0
  ancestors[1,] <- seq(1:N)
  #t=2:T
  #2. conditional sample
  if(Time > 2){
    for(t in 2:(Time - 1)){

      #print(t)


      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        ancestors[t,] <- resample(w[t-1,])
        #mx <- max(w[t-1,  ])
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        add <- rep(0, N)

      }else{
        ancestors[t,] <- c(1:N)
        add <- w[t-1,]
      }

      for(i in 1:N){

        X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],],
                                             model, psi_pa[t,], 1)
        log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)

        w[t,i] <- eval_twisted_potential(model,
                                         list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i]) + add[i]
        #w[t,i] <- log(g_aux(obs[t,], X[t,i,], t, psi_pa)) + add[i]
      }


    }
  }

  if(Time >= 2){
    t = Time


    if(compute_ESS_log(w[t-1,]) <= kappa*N){
      ancestors[t,] <- resample(w[t-1,])
      #mx <- max(w[t-1,  ])
      logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
      add <- rep(0, N)

    }else{
      ancestors[t,] <- c(1:N)
      add <- w[t-1,]
    }

    for(i in 1:N){

      X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],],
                                           model, psi_pa[t,], 1)
      log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)

      w[t,i] <- eval_twisted_potential(model,
                                       list(NA, NA, psi_pa[t,]), X[t,i,], log_likelihoods[t,i]) + add[i]
      #w[t,i] <- log(g_aux(obs[t,], X[t,i,], t, psi_pa)) + add[i]
    }

    #mx <- max(w[t,  ])

  }
  logZ = logZ + normalise_weights_in_log_space(w[t,])[[2]]


  #log(mean(exp(w[t,  ]-mx))) + mx


  return(list(X = X, w = w, logZ = logZ, log_likelihoods = log_likelihoods))
}

