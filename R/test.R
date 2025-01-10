test <- function(model, obs, w, X, N, psi_pa, psi_l, psi_u, logZ, ancestors){
logZ = 0 #logZ - normalise_weights_in_log_space(w[psi_l,])[[2]]
obs_params <- model$obs_params

log_likelihoods <- matrix(NA, Time, N)


    for(t in (psi_l + 1):(psi_u - 1)){
      #print(t)
      if(compute_ESS_log(w[t-1,]) <= 50){
        
        #re = re + 1
        #if(t != (psi_l + 1)){
          ancestors[t,] <- resample(w[t-1,])
          logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        #}
        
        
       
        #resample_time <- c(resample_time, t-1)
        cat('logZ =', logZ, 're= ', t, "add =", normalise_weights_in_log_space(w[t-1,])[[2]])
        
        for(i in 1:N){
          
          X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],], model, psi_pa[t,], 1)
          log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
          if(t == (psi_l + 1)){
           w[t,i] <- eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i])+  w[t - 1,i] 
          }else{
            w[t,i] <- eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
          }
        }
        
        
      }else{
        
        ancestors[t,] <- 1:N
        for(i in 1:N){
        
          X[t,i,] <- as.matrix(sample_twisted_transition(X[t-1, i,], model, psi_pa[t,], 1))
         
          log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
          
          w[t,i] <- w[t-1,i] + eval_twisted_potential(model, 
                   list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
         
        }
      }
      
    }
  
  
  
  t = psi_u
  #resample_time <- c(resample_time, t)
  if(compute_ESS_log(w[t-1,]) <= 50){
    re = re + 1
    
    
    ancestors[t,] <- resample(w[t-1,])
    logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
    cat('logZ =', logZ, 're= ', t, "add =", normalise_weights_in_log_space(w[t-1,])[[2]])
    
    for(i in 1:N){
      #filtering particles
      X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],], model, psi_pa[t,], 1)
      log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
      w[t,i] <- eval_twisted_potential(model, list(NA, NA, psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
      
    }
  }else{
 
    ancestors[t,] <- 1:N #?
    for(i in 1:N){
      
      X[t,i,] <- sample_twisted_transition(X[t-1, i,], model, psi_pa[t,], 1)
      log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
      w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(NA, NA, psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
    }
  }
  
  
  
  logZ <- logZ #+ normalise_weights_in_log_space(w[t,])[[2]]
  if(t == Time){
    logZ <- logZ + normalise_weights_in_log_space(w[t,])[[2]]
  }
  w[t,] <- w[t,]  - normalise_weights_in_log_space(w[t,])[[2]]
  cat('logZ =', logZ, 're= ', t, "add =", normalise_weights_in_log_space(w[t,])[[2]])
return(list(X, w, logZ, ancestors))
}

