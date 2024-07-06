set.seed(1234)
library(mvnfast)
library(FKF)
library(roxygen2)
Napf = N = 200
lag = 10
Time = 200
d_ = 2
alpha = 0.42
a <- matrix(nrow = d_, ncol = d_)
for (i in 1:d_){
  for (j in 1:d_){
    a[i,j] = alpha^(abs(i-j) + 1)
  }
}
ini <- rep(0, d_)
b = c = D_ = ini_c = diag(1, nrow = d_, ncol = d_)
k_ <- 5
tau_ <- 0.5
kappa_ = 0.5
dt_ <- ct_ <- matrix(0, d_, 1)
Tt_ <- a
P0_ <- Zt_ <- Ht_ <- Gt_ <- diag(1, d_, d_)
a0_ <- rep(0, d_)
model <- list(ini_mu = ini, ini_cov = ini_c, A = a, B = b, C = c, D = D_, k = k_,
              tau = tau_, kappa = kappa_, dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_, Zt = Zt_,
              Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)

obs_ <- sample_obs(model, Time)

output <- generate_blocks(lag, Time)
breaks_ <- output[[1]]
psi_index_ <- output[[2]]

data <- list(obs = obs_, breaks = breaks_, psi_index = psi_index_)

output <- compute_filtering(model, data)
fkf.obj_ <- output[[1]]
fks.obj_ <- output[[2]]

kalman <- list(fkf.obj = fkf.obj_, fks.obj  = fks.obj_ )

ratio <- vector()
dist <- vector()
for (qq in 1:1) {
  set.seed(qq*2)
  output <- run_quasi_online_pf(model, data, lag, Napf, N, kalman)
  ratio[qq] <- output[[5]]
  dist[qq] <- output[[6]]
}

