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
  tran_mu <- rep(as.matrix(model$tran_mu)[1,1], d)
  tran_cov <- model$tran_c
  
  X[1,,] <- mvnfast::rmvn(Napf, ini, ini_c)
  for(i in 1:Napf){
    w[1,i] <- model$eval_likelihood(X[1,i,], obs[1,], model$obs_params)
  }
  logZ = 0
  
  logZ <- logZ + normalise_weights_in_log_space(w[1,])[[2]]
  mix <- resample(w[1,], mode = 'multi')
  X[1,,] <- X[1,mix,]
  ancestors[1,] <- seq(1:N)
  for (t in 2:Time) {
    
    for(i in 1:Napf){
     
      X[t,i,] <-  sample_transition(X[t-1,i,], model, 1)
      w[t,i] <-  model$eval_likelihood(X[t,i,], obs[t,], model$obs_params)
    }
    
    
    logZ <- logZ + normalise_weights_in_log_space(w[t,])[[2]]
    mix <- resample(w[t,], mode = 'multi')
    X[t,,] <- X[t,mix,]
    ancestors[t,] <- mix
    #V[t,] <- V[t-1,]
    #V[t, -unique(I[t-1,mix])] <- t
    
    
  }
  #logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]
  
  return(list(X = X, w = w, logZ = logZ))
}

library(mvnfast)
start = proc.time()
set.seed(123)
library(FKF)
d=1
Num <- 200
Time = 200
alpha = 0.42
B = C = D = diag(1, nrow = d, ncol = d)
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}

dt <- ct <- matrix(0,d,1)
Tt <- A
P0 <- Zt <- Ht <- Gt <- diag(1,d,d)
a0 <- rep(0,d)

X_true <- matrix(NA, Time, d)
final1 <- vector()

#avg <- matrix(nrow=50, ncol=Time)
Obs <- function(){
  #set.seed(seed)
  X_true[1,] <- rnorm(d) 
  for(t in 2:Time){  #observations
    #set.seed(seed)
    X_true[t,] <- rnorm(d) + as.vector(A%*%X_true[t-1,])  #t(rnorm(d) + A%*%x)
  }
  #set.seed(seed)
  return(matrix(rnorm(Time*d, X_true, 1), ncol = d))
}
