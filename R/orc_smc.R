#' Rolling window function
#'
#' 
#'
#' @param lag 
#' @param data 
#' @param model 
#' @param N
#'
#' @return A matrix of observations with dimensions Time x d.
#' @export
Orc_SMC <- function(lag, data, model, N) {
  obs <- data$obs
  Time <- nrow(obs)
  d <- ncol(obs)
  
  A <- model$tran_mu
  B <- model$tran_cov
  obs_params <- model$obs_params
  K <- model$parameters$k
  
  X0 <- matrix(0, nrow = N, ncol = d)
  w0 <- matrix(log(1/N), 1, N)
  
  logZ_vec <- numeric(Time)
  logZ_tilde <- numeric(Time)
  logZ_s  <- 0
  
  X <- X_apf <- array(NA, c(Time, N, d))
  log_W_apf <- log_W <- matrix(NA, Time, N)
  log_likelihoods <- log_likelihoods_apf <- matrix(NA, Time, N)
  filtering_estimates <- matrix(NA, Time, d)
  
  psi_pa <- matrix(NA, Time, 2*d)
  
  H <- vector('list', Time + 1)
  H_tilde <- vector('list', Time + 1)
  
  H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  H_tilde[[1]] <- list(X = X0, logW = w0, logZ = 0)
  
  for (t in 1:Time) {
    print(t)
    t0 <- max(t - lag + 1, 1)
    
    # Step 3: Init pass with psi ≡ 1
    
    output <- run_psi_APF_rolling(data, t, psi_pa, H_tilde[[t]] , model, init = TRUE)
    H_tilde[[t+1]] <- output$H
    log_likelihoods_apf[t,] <- output$log_likelihoods
    
    # Step 4: Policy Refinement
    for (k in 1:K) {
      
      for (s in t0:t) {X_apf[s,,] <- H_tilde[[s + 1]]$X}
      
      psi_pa[t0:t,] <- learn_psi(X_apf[t0:t,,, drop = FALSE],
                                 model, log_likelihoods_apf[t0:t,, drop = FALSE])
      
      for (s in t0:t) {
        output <- run_psi_APF_rolling(data, s, psi_pa, H_tilde[[s]], model, init = FALSE)
        
        H_tilde[[s+1]] <- output$H
        log_likelihoods_apf[s,] <- output$log_likelihoods
      }
    }
    
    # Step 5: Final pass to get filtering distributions + logZ

    for (s in t0:t) {
      
      output <- run_psi_APF_rolling(data, s, psi_pa, H[[s]], model, init = FALSE)
      
      H[[s+1]] <- output$H
     
      #W_t <- normalise_weights_in_log_space(log_W[s,])[[1]]
      #filtering_estimates[s,] <- colSums(W_t * X[s,,])
      
    }
    
  
    logZ_vec[t] <- H[[t + 1]]$logZ
    
  }
  
  return(list(logZ = logZ_vec, f_means = filtering_estimates))
}
