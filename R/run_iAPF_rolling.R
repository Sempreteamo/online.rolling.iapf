#' Function of psi_APF
#'
#' This function performs sampling and resampling procedures of
#' psi-Auxiliary Particle Filter with $/kappa$-adaptive resampling
#'
#' @param t time
#' @param data List containing required data information on which to run the filter
#' @param H_prev previous information
#' @param psi_pa Parameters of the twisting function
#' @param init Index to decide whether the algorithm is at its initialized block
#' @param model A list containing information of model
#'
#' @return A list containing:
#' 
#'
#' @export
#'
run_psi_APF_rolling <- function(data, t, psi_pa, H_prev, model, init) {
  
  # Extract previous values
  obs <- data$obs
  X_prev <- H_prev$X
  logW_prev <- H_prev$logW
  logZ_t <- H_prev$logZ
  #logZ_prev <- H_prev$logZ
  #logZ_t = 0
  kappa <- model$parameters$kappa
  N <- length(logW_prev)  # Number of particles
  
  ini <-model$ini
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  
  A <- model$tran_mu
  B <- model$tran_cov
  obs_params <- model$obs_params
  d <- ncol(A)
  
  X_new <- matrix(NA, N, d)
  
  # Step 1: Compute v^n in log domain
  log_psi_tilde <- rep(NA, N)
  for(i in 1:N){
    log_psi_tilde[i] <-  evaluate_psi_tilde(X_prev[i,, drop = FALSE], psi_pa[t,], model)
  }
  
  log_v <- logW_prev + log_psi_tilde  # log(W_t-1) + log(f_t)
  
  
  # Step 2: Compute logZ_t and normalized logV
  logZ_t <- logZ_t + log_sum_exp(log_v)
  logV <- log_v - log_sum_exp(log_v)
  
  # Step 3: Compute ESS and resample if necessary
  ESS <- exp(-log_sum_exp(2 * logV))
  
  if (ESS < kappa * N) {
    ancestors <- resample(logV)
    cat('re at ', t)
    logV <- rep(-log(N), N)  # Reset logV after resampling
    #add = rep(0, N) #??
  } else {
    ancestors <- 1:N  # No resampling needed
    #add = logW_prev #??
  }
  
  # Step 4: Sample new states using f_t^ψ
  log_likelihoods <- rep(NA, N)
  log_g <- rep(NA, N)
  log_psi_t <- rep(NA, N)
  
  
  for(i in 1:N){
    if (init == TRUE && t == 1) {
      # Case 1: Initialisation step at time t = 1
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu, model$ini_cov)
      
    } else if (init == TRUE && t != 1) {
      # Case 2: Initialisation at time t > 1
      X_new[i, ] <- mvnfast::rmvn(1, ini_mu + A %*% (X_prev[ancestors[i], ] - ini_mu), B)
      
    } else if (init == FALSE && t == 1) {
      # Case 3: Non-initialisation but t = 1 (should rarely happen)
      X_new[i, ] <- sample_twisted_initial(list(mean = ini_mu, cov = as.matrix(ini_cov)[1,1]), psi_pa[t,], 1)
      
    } else {
      # Case 4: Non-initialisation and t > 1 (normal online step)
      X_new[i, ] <- sample_twisted_transition(X_prev[ancestors[i], ], model, psi_pa[t, ], 1)
    }
    
    log_likelihoods[i] <- model$eval_likelihood(X_new[i,], obs[t,, drop = FALSE], obs_params)
    log_psi_t[i] <- evaluate_psi(X_new[i,], psi_pa[t,])
  }
  
  
  
  log_g <- log_likelihoods #+ add 
  
  # Step 5: Compute log w_t
  log_w <- logV + log_g - log_psi_t  # log(V_t) + log(g) - log(ψ_t)
  
  # Step 6: Compute log Z_t and normalize weights
  logZ_t <- logZ_t + log_sum_exp(log_w)
  logW <- log_w - log_sum_exp(log_w)  # Self-normalized log weights
  
  # Return updated particle system
  return(list(X = X_new, logW = logW, logZ = logZ_t, log_likelihoods = log_likelihoods))
}


#' @import mvnfast
