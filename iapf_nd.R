library(mvtnorm)
library(MASS)
library(FKF)
library(mvnfast)


Num <- function(Z, l, k){
  return(sd(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l])))/mean(exp(Z[max(l-k,1):l]-max(Z[max(l-k,1):l]))))
}

####APF function####
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

iAPF <- function(data, Napf, model){
  N <- vector()
  N[1] <- Napf
  k = 5
  tau <- 0.5
  Z <- vector()
  H <- vector('list', Time)
  
  #H[[1]] <- list(X = X0, logW = w0, logZ = 0)
  
  index = 1
  l = 1  #actually 0
  obs <- data$obs
  d = ncol(obs)
  Time <- nrow(obs)
  X <- array(NA, dim = c(Time, N, d))
  w <- matrix(NA, Time, N)
  kappa <- model$parameters$kappa
  
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  A <- model$tran_mu
  B <- model$tran_cov
  obs_params <- model$obs_params
  log_likelihoods <- matrix(NA, Time, N)
  ancestors <- matrix(NA, Time, N)
  psi_pa <- matrix(NA, Time + 1, 2*d)
  
  X[1, ,] <- rnorm(N[l]*d)  #particles
  for(i in 1:N[l]){
    w[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
     # log(g(obs[1,], X[1,i,]))  #weights
  }
  
  ancestors[1,] <- seq(1:N)
  log_likelihoods[1,] <- w[1,]
  Z[l] <- 0
  
  H[[1]] <- list(X = X[1, ,], logW =  w[1,], log_li = log_likelihoods[1,])
  #t=2:T
  #2. conditional sample
  #re = 0
  if(Time >= 2){
    for(t in 2:Time){
      #print(t)
      
      #a)
      if(compute_ESS_log(w[t-1,]) <= kappa*N[l]){
        ancestors[t,] <- resample(w[t-1,])
        Z[l] = Z[l] + normalise_weights_in_log_space(w[t-1,])[[2]]
        add <- rep(0, N)
        
      }else{
        ancestors[t,] <- c(1:N[l])
        add <- w[t-1,]
      }
      
      for(i in 1:N[l]){
        X[t, i,] <- mvnfast::rmvn(1, ini_mu + 
                                    A%*%(as.vector(X[t-1, ancestors[t,i],]) - ini_mu), B)
        log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], 
                                                      obs[t,, drop = FALSE], obs_params)
      }
      
      w[t,] <- log_likelihoods[t,] + add
      
      H[[t]] <- list(X = X[t, ,], logW =  w[t,], log_li = log_likelihoods[t,])
    }
    
    Z[l] = Z[l] + normalise_weights_in_log_space(w[t,])[[2]]
  }
  

  #print(paste0('re=',re))
  
  ####iapf####
  while(index){
    
    print(l)
    #a)
    output <- list()
    
    if(l != 1){
      output <- APF(psi_pa, data, N[l], model)
      X <- output$X
      w <- output$w
      Z[l] <- output$logZ
      log_likelihoods <- output$log_likelihoods
      
      for(t in 1:Time){
        H[[t]] <- list(X = X[t, ,], logW =  w[t,], log_li = log_likelihoods[t,])
      }
    }
    
    #b)
    print(Z[l])
    if(l <= k | (Num(Z, l, k) >= tau)){
      #psi^{l+1}
      
      for(s in Time:1){
        psi_pa[s,] <- learn_psi(s, psi_pa[s+1,, drop = FALSE ], H[[s]], model)
      }
      
      
      
      if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z[max(l-k,1):l])){
        N[l+1] <- 2*N[l]
        
      }else{
        N[l+1] <- N[l]
      }
      
      #print(paste0('Z[l]=',Z[l]))
      
      #print(paste0('Z=',fkf.obj))
      
      l <- l+1
    }else break
    
  }
  
  #3.
  #output <- APF(psi_pa, data, N[l], model)
  Z[l] <- output$logZ
  #psi_pa <- output[[5]]
  Z_appro <- Z[l]
  print( Z_appro)
 return(list(Z = Z_appro, psi = psi_pa))
}
