#' Function of psi_APF
#'
#' This function performs sampling and resampling procedures of
#' psi-Auxiliary Particle Filter with $/kappa$-adaptive resampling
#'
#' @param model List containing model parameters
#' @param data List containing time series on which to run the filter
#' @param N Number of particles to use
#' @param psi_pa Parameters of the twisting function
#' @param init Index to decide whether the algorithm is at its initialization state
#'
#' @return A list containing:
#' x is the particle set generated during the sampling and resampling approach
#' w is the weights of particles
#' ancestors is indices associated with every resampling event
#' logZ is normalizing constant estimate
#' resample_time is the specific time that resample happens
#'
#' @export
#'
run_psi_APF <- function(model, data, N, psi_pa, init, jump_ini = FALSE){
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  A <- model$tran_mu
  B <- model$tran_cov
  obs_params <- model$obs_params
  re = 0
  
  obs <- as.matrix(data[[1]])
  breaks <- data[[2]]
  w_previous <- data[[3]]
  X_previous <- data[[4]]
  
 
  
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
  if(init){
    if(breaks[1] == 1){
      #the first block. break controls which block the algorithm is running

      X[1, ,] <- mvnfast::rmvn(N , ini_mu, ini_cov)
      
      for(i in 1:N){
        w[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
#print(w[1,i])
      }
      

    }else{

      s <- resample(w_previous, mode = 'multi')
      for(i in 1:N){
        X[1,i,] <- mvnfast::rmvn(1, ini_mu + A%*%(as.vector(X_previous[s[i],, drop = FALSE]) - ini_mu), B)
        w[1, i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
      }

    }

    ancestors[1,] <- seq(1:N)
    log_likelihoods[1,] <- w[1,]

    for(t in 2:Time){

      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        re = re + 1

        ancestors[t,] <- resample(w[t-1,], mode = 'multi')
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        resample_time <- c(resample_time, t-1)
#print(normalise_weights_in_log_space(w[t-1,])[[2]])


        for(i in 1:N){
          X[t,i,] <- mvnfast::rmvn(1, ini_mu + A%*%(X[t-1, ancestors[t,i],] - ini_mu), B)
            #mvnfast::rmvn(1, A%*%X[t-1, ancestors[t,i],], B)
          w[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
        }

        log_likelihoods[t,] <- w[t,]

      }else{
        ancestors[t,] <- 1:N


        for(i in 1:N){
          X[t,i,] <-  mvnfast::rmvn(1, ini_mu + A%*%(X[t-1, i,] - ini_mu), B)
            #mvnfast::rmvn(1, A%*%X[t-1, i,], B)
          log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
        }
        w[t,] <- w[t-1,] + log_likelihoods[t,]
      }
      #print(t)
      #print(compute_ESS_log(w[t-1,]))
    }

    #print(re)
    resample_time <- c(resample_time, Time)


    logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]

  }else{
   if(jump_ini == FALSE){
     if(breaks[1] == 1){ #new
       #if(breaks[1] == lag){
       
       
       X[1,,] <- sample_twisted_initial(list(mean = ini_mu, cov = as.matrix(ini_cov)[1,1]), psi_pa[1,], N)
       
       for(i in 1:N){
         log_likelihoods[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
         w[1,i] <- eval_twisted_potential(model, list(psi_pa[1,], psi_pa[2,], psi_pa[1,]), X[1,i,],  log_likelihoods[1,i])
         #}
       } #new 
       
     }else{
       #print(1)
       w_adj <- vector()
       
       for(i in 1:N){
         X_previous[i,] <- as.vector(X_previous[i,, drop = FALSE])
         w_adj[i] <- w_previous[i]*exp(1/2*(t(ini_mu + A%*%(X_previous[i,]) - ini_mu) + t(psi_pa[1,1:d])%*%
                                              diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d))%*%diag((psi_pa[1, (d+1):(d+d)]^(-2) + 1)^(-1), nrow=d,ncol=d)%*%
                                         (ini_mu + A%*%(X_previous[i, ] - ini_mu) + diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d)%*%psi_pa[1,1:d]) -
                                         1/2*(t(ini_mu + A%*%(X_previous[i, ]- ini_mu))%*%(ini_mu + A%*%X_previous[i,] - ini_mu) +
                                                t(psi_pa[1,1:d])%*%diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d)%*%psi_pa[1,1:d]))
       }
       
       s <- resample(w_adj, mode = 'multi')
       
       for (i in 1:N){
         X[1, i, ] <- sample_twisted_transition(as.vector(X_previous[s[i],, drop = FALSE]), model, psi_pa[1,], 1)
         log_likelihoods[1,i] <- model$eval_likelihood(X[1,i,], obs[1,, drop = FALSE], obs_params)
         w[1, i] <- eval_twisted_potential(model, list(psi_pa[1,], psi_pa[2,], psi_pa[1,]), X[1,i,], log_likelihoods[1,i])
       }
       
     }
   }
    
    if(jump_ini){
      psi_l = breaks[1]
      Time = breaks[2]
      X[psi_l - 1,,] <- X_previous
      w[psi_l - 1,] <- w_previous
    }else{
      psi_l = 2
    }
    
    ancestors[1,] <- seq(1:N)

    if(Time > 2){ 
      for(t in psi_l:(Time - 1)){
        
      #print(Time)
      #cat('t=', t)
      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        re = re + 1
        
        ancestors[t,] <- resample(w[t-1,], mode = 'multi')
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        resample_time <- c(resample_time, t-1)
        
        for(i in 1:N){
          #print(psi_pa[t+1,])
          X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],], model, psi_pa[t,], 1)
          log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
          w[t,i] <- eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
          
        }
        
        
      }else{
       
        ancestors[t,] <- 1:N
        for(i in 1:N){
          #print(psi_pa[t+1,])
          X[t,i,] <- sample_twisted_transition(X[t-1, i,], model, psi_pa[t,], 1)
          #print(obs_params$obs_mean)
          #print(matrix(X[t,i,], 1))
          #print(obs_params$obs_mean%*%matrix(X[t,i,], 1))
          log_likelihoods[t,i] <- model$eval_likelihood(X[t,i,], obs[t,, drop = FALSE], obs_params)
          
          w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], log_likelihoods[t,i])
        }
      }
      
      }
    }
  

    t = Time
    resample_time <- c(resample_time, t)

    if(compute_ESS_log(w[t-1,]) <= kappa*N){
      re = re + 1

      ancestors[t,] <- resample(w[t-1,], mode = 'multi')
      logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]

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

    #print(re)

    logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]
  }

  for (t in 1:Time){
   avg[, t] <- length(unique(X[t,,,drop = FALSE]))/d
  }


  return(list(X, w, logZ, ancestors, normalise_weights_in_log_space(w[Time,])[[2]], log_likelihoods, avg))
}
#' @import mvnfast
#' @import stats
