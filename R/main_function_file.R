#' @examples
#' \dontrun{
#' The following parameters are provided by users
Napf = N = 200
lag = 50
Time = 50
d_ = 1
alpha = 0.42
tran_m <- matrix(nrow = d_, ncol = d_)
for (i in 1:d_){
  for (j in 1:d_){
    tran_m[i,j] = alpha^(abs(i-j) + 1)
  }
}
ini <- rep(0, d_)
ini_c  <- alpha^2 / (1-alpha)^2
tran_c = den_mean = den_cov = diag(1, nrow = d_, ncol = d_)
parameters_ <- list(k = 5, tau = 0.5, kappa = 0.5)

output <- generate_blocks(lag, Time)
breaks_ <- output[[1]]
psi_index_ <- output[[2]]

obs_ <- sample_obs(tran_m, tran_c, den_mean, den_cov, Time, dist = 'lg') #provided by users

dt_ <- ct_ <- matrix(0, d_, 1)
Tt_ <- a
P0_ <- Zt_ <- Ht_ <- Gt_ <- diag(1, d_, d_)
a0_ <- rep(0, d_)
params <- list(dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_, Zt = Zt_,
               Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)

output <- compute_filtering(params, obs_)
fkf.obj_ <- output[[1]]
fks.obj_ <- output[[2]]

model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m, tran_cov = tran_c,
eval_likelihood = evaluate_likelihood(), parameters = parameters_, dist = 'lg')

data <- list(obs = obs_, breaks = breaks_, psi_index = psi_index_)

kalman <- list(fkf.obj = fkf.obj_, fks.obj  = fks.obj_ ) #provided by users

#run the algorithm
output <- run_quasi_online_pf(model, data, lag, Napf, N, kalman)
#' }

