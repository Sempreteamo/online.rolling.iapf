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
  logZ_s  <- 0
  
  X <- X_apf <- array(NA, c(Time, N, d))
  log_W_apf <- log_W <- matrix(NA, Time, N)
  log_likelihoods <- log_likelihoods_apf <- matrix(NA, Time, N)
  filtering_estimates <- matrix(NA, Time, d)
  
  psi_pa <- matrix(NA, Time, 2*d)
  
  for (t in 1:Time) {
    print(t)
    t0 <- max(t - lag + 1, 1)
    
    # Step 3: Init pass with psi ≡ 1
    if (t == 1) {
      H_prev <- list(X = X0, logW = w0)
    } else {
      H_prev <- list(X = X_apf[t-1,,], logW = log_W_apf[t-1,])
    }
    
    output <- run_psi_APF_rolling(data, t, psi_pa, H_prev, model, init = TRUE)
    X_apf[t,,] <- output$X
    log_W_apf[t,] <- output$logW
    log_likelihoods_apf[t,] <- output$log_likelihoods
    
    # Step 4: Policy Refinement
    for (k in 1:K) {
      psi_pa[t0:t,] <- learn_psi(X_apf[t0:t,,, drop = FALSE], obs[t0:t,, drop = FALSE],
                                 model, log_likelihoods_apf[t0:t,, drop = FALSE])
      
      logZ_t = 0
      for (s in t0:t) {
        if (s == 1) {
          H_prev <- list(X = X0, logW = w0, logZ = 0)
        } else {
          H_prev <- list(X = X_apf[s - 1,,], logW = log_W_apf[s - 1,], logZ = logZ_t)
        }
        
        output <- run_psi_APF_rolling(data, s, psi_pa, H_prev, model, init = FALSE)
        X_apf[s,,] <- output$X
        log_W_apf[s,] <- output$logW
        log_likelihoods_apf[s,] <- output$log_likelihoods
        logZ_t <- output$logZ
      }
    }
    
    # Step 5: Final pass to get filtering distributions + logZ
    H_prev <- if (t0 == 1) list(X = X0, logW = w0, logZ = 0) else list(X 
                    = X[t0 - 1,, ], logW = log_W[t0 - 1,], logZ = 0)
    
    logZ_s = 0
    
    for (s in t0:t) {
      
      output <- run_psi_APF_rolling(data, s, psi_pa, H_prev, model, init = FALSE)
      X[s,,] <- output$X
      log_W[s,] <- output$logW
      logZ_s <- output$logZ
      W_t <- normalise_weights_in_log_space(log_W[s,])[[1]]
      filtering_estimates[s,] <- colSums(W_t * X[s,,])
      
      H_prev <- list(X = X[s,,], logW = log_W[s,], logZ = logZ_s)
    }
    
    logZ_vec[t] <- if (t0 == 1) logZ_s  else logZ_vec[t0 - 1] + logZ_s
  }
  
  return(list(logZ = logZ_vec, f_means = filtering_estimates))
}
