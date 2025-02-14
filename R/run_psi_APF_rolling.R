#' Function of psi_APF
#'
#' This function performs sampling and resampling procedures of
#' psi-Auxiliary Particle Filter with $/kappa$-adaptive resampling
#'
#' @param model List containing model parameters
#' @param data List containing required data information on which to run the filter
#' @param N Number of particles to use
#' @param psi_pa Parameters of the twisting function
#' @param init Index to decide whether the algorithm is at its initialized block
#' @param jump_ini Index to decide whether the algorithm is at the middle stage of the psi-APF
#'
#' @return A list containing:
#' X is the particle set generated during the sampling and resampling approach
#' w is the weights of particles
#' logZ is normalizing constant estimate
#' ancestors is indices associated with every resampling event
#' log_likelihoods is the log likelihoods for each particle X given likelihood function
#'
#' @export
#'
run_psi_APF_rolling <- function(model, data, N, psi_pa, init){
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  A <- model$tran_mu
  B <- model$tran_cov
  obs_params <- model$obs_params
  re = 0
  
  obs <- as.matrix(data[[1]])
  run_block <- data[[2]]
  b_s <- run_block[1]
  b_e <- run_block[2]
  
  X_init <- data[[3]]
  w_init <- data[[4]]
  
  Time <- nrow(obs)
  d = nrow(as.matrix(A))
  kappa <- model$parameters$kappa
  ancestors <- matrix(NA, Time, N)
  resample_time <- vector()
  
  X <- array(NA, dim = c(Time, N, d))
  w <- matrix(NA, Time, N)
  logZ <- 0
  log_likelihoods <- matrix(NA, Time, N)
  
  avg <- matrix(nrow = 1, ncol = Time)
  #init controls the the condition where we conduct the initialization

  
  #initialization with standard distribution
  if(init == TRUE){
    for(t in 1:Time){
      ancestors[t,] <- 1:N
      
      for(i in 1:N){
        X[t,i,] <-  mvnfast::rmvn(1, ini_mu + A%*%(X_init[i,drop = FALSE] - ini_mu), B)
        #mvnfast::rmvn(1, A%*%X[t-1, i,], B)
        log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
      }
      w[t,] <- log_likelihoods[t,]
    }
  }
  #loglikelihoods only has 1 row 
  
  #twisted distribution for t = 1#
  if(init != TRUE){
    index = 1
    for(t in b_s:b_e){
      if(t == b_s){
        w_t = w_init
      }else{
        w_t = w[index - 1,]
      }
      
      if(compute_ESS_log(w_t) <= kappa*N){
        #re = re + 1
        ancestors[index,] <- resample(w_t)
        logZ = logZ + normalise_weights_in_log_space(w_t)[[2]]
        add = rep(0, N)
      }else{
        ancestors[index,] <- 1:N
        add = w_t
      }
      
      for(i in 1:N){
        #filtering particles
        X[index,i,] <- sample_twisted_transition(X_init[ancestors[index,i], drop = FALSE], model, psi_pa[t,], 1)
        log_likelihoods[index,i] <- model$eval_likelihood(X[index,i,], obs[index,, drop = FALSE], obs_params)
        
        w[index,i] <- add[i] + eval_twisted_potential(model, list(NA, psi_pa[t,], psi_pa[t,]), X[index,i,], log_likelihoods[index,i])
        
      }
      index = index + 1
    }

  }
  
  
  
  for (t in 1:Time){
    avg[, t] <- length(unique(X[t,,,drop = FALSE]))/d
  }
  
  
  return(list(X, w, logZ, ancestors, log_likelihoods,  normalise_weights_in_log_space(w[Time,])[[2]]))
}
#' @import mvnfast
#' @import stats
