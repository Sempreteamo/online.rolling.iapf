seed = 1

Napf = N = 100

Time = 50
lag = 20
d_ = 4
# 2^((seed - 1) %/% 40)

#ini <- rep(0, d_)

tran_m = diag(1, nrow = d_, ncol = d_)
tran_c = 10^(-3) * matrix(
  c(9.65,  11.42,  3.97,  12.07,
    11.42, 20.43,  5.44,  21.08,
    3.97,  5.44,   5.45,  7.08,
    12.07, 21.08,  7.08,  22.31),
  nrow = 4,
  ncol = 4,
  byrow = TRUE  #
)

#ini_c <- diag(1, d_, d_)
#ini_c = diag(1, nrow = d_, ncol = d_)
#obs_m = diag(1, nrow = d_, ncol = d_)
#obs_c = diag(1, nrow = d_, ncol = d_)
parameters_ <- list(k = 6, tau = 0.5, kappa = 0.5)

obs_p <- list(c = -1.27, cov = matrix(c(1.000, 0.404, 0.280, 0.347,
                                        0.404, 1.000, 0.400, 0.541,
                                        0.280, 0.400, 1.000, 0.622,
                                        0.347, 0.541, 0.622, 1.000),
                                      nrow = d_, byrow = TRUE))

model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m, tran_cov = tran_c, obs_params = obs_p,
              eval_likelihood = evaluate_likelihood_msv, simu_observation = simulate_observation_msv,
              parameters = parameters_)

obs_ <- final_df[1:50, ,drop = FALSE] # provided by users
#obs_ <- sample_obs(model, Time, d_)

a0_ = ini     # Initial state mean
P0_ = ini_c    # Initial state covariance
   # Observation matrix (C)
Ht_ = tran_c    # Observation noise covariance (R)
Gt_ = obs_p[[2]]   # Process noise covariance (Q)

dt_ <- matrix(0, nrow = d_, ncol = 1)
ct_ <- matrix(rep(-1.27, d_), nrow = d_, ncol = 1)

Tt_ <- diag(d_)
Zt_ <- diag(d_)

params <- list(dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_, Zt = Zt_,
               Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)

filter <- compute_fkf(params, obs_)
filtering <- filter[[1]]
smoothing <- filter[[2]]

data_ =  data <- list(obs = obs_)
#data <- list(obs = obs_, breaks = breaks_,
#            psi_index = psi_index_)


filtering_estimates <- 0

logZ_matrix_rolling <- vector()
logZ_matrix_rolling_iapf <- vector()

set.seed(seed)

time_iapf <- system.time({
  output_iapf <- iAPF(data, Napf, model)
  output <-  Orc_SMC(lag, data, model, Napf)
})[["elapsed"]]
print(time_iapf)

H <- output$H_forward
logZ_matrix_rolling_iapf <- output$logZ


#filtering_estimates <- output$f_means

log_ratio_rolling_iapf <- compute_log_ratio(logZ_matrix_rolling_iapf[Time], filtering)
print(log_ratio_rolling_iapf)

smoothed_particles <- compute_backward_smoothing(H)

fkf_smoothed_means <- t(filtering$att)
fkf_smoothed_covs  <- filtering$Ptt
smooth_means <- t(smoothing$ahatt)
smooth_covs <- smoothing$Vt

plot_data_list <- list()

smoothed_means <- apply(smoothed_particles, c(2, 3), mean)
smoothed_quantiles_lower <- apply(smoothed_particles, c(2, 3), function(x) quantile(x, 0.05))
smoothed_quantiles_upper <- apply(smoothed_particles, c(2, 3), function(x) quantile(x, 0.95))


for (i in 1:d_) {
  #
  true_state_df <- data.frame(
    Time = 1:Time,
    Value = true_states[, i],
    Type = "True State",
    Coordinate = paste0("Coordinate ", i)
  )
  plot_data_list[[paste0("true_state_", i)]] <- true_state_df

  # 粒子平滑均值
  pf_smooth_mean_df <- data.frame(
    Time = 1:Time,
    Value = smoothed_means[, i],
    Type = "PF Smooth Mean",
    Coordinate = paste0("Coordinate ", i)
  )
  plot_data_list[[paste0("pf_smooth_mean_", i)]] <- pf_smooth_mean_df

  #
  pf_smooth_ci_df <- data.frame(
    Time = 1:Time,
    Lower = smoothed_quantiles_lower[, i],
    Upper = smoothed_quantiles_upper[, i],
    Type = "PF Smooth CI",
    Coordinate = paste0("Coordinate ", i)
  )
  plot_data_list[[paste0("pf_smooth_ci_", i)]] <- pf_smooth_ci_df

  # FKF 平滑均值
  fkf_smooth_mean_df <- data.frame(
    Time = 1:Time,
    Value = fkf_smoothed_means[i, ],
    Type = "FKF Smooth Mean",
    Coordinate = paste0("Coordinate ", i)
  )
  plot_data_list[[paste0("fkf_smooth_mean_", i)]] <- fkf_smooth_mean_df

  # FKF 平滑置信区间 (例如 95% CI，假设是高斯分布)
  # 计算标准差 sqrt(diag(covariance))
  fkf_smooth_sd <- sqrt(sapply(1:Time, function(t) fkf_smoothed_covs[i, i, t]))
  fkf_smooth_ci_df <- data.frame(
    Time = 1:Time,
    Lower = fkf_smoothed_means[i, ] - 1.96 * fkf_smooth_sd,
    Upper = fkf_smoothed_means[i, ] + 1.96 * fkf_smooth_sd,
    Type = "FKF Smooth CI",
    Coordinate = paste0("Coordinate ", i)
  )
  plot_data_list[[paste0("fkf_smooth_ci_", i)]] <- fkf_smooth_ci_df
}
