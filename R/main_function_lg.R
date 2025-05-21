#' @examples
#' \dontrun{
#' The following parameters are provided by users
#' library(mvnfast)
#' library(FKF)
#' Napf = N = 200
#' lag = 4
#' Time = 10
#' d_ = 2
#'
#' alpha = 0.42
#' tran_m <- matrix(nrow = d_, ncol = d_)
#' for (i in 1:d_){
#'   for (j in 1:d_){
#'       tran_m[i,j] = alpha^(abs(i-j) + 1)
#'   }
#' }
#' ini <- rep(0, d_)
#'
#' #
#' tran_m = diag(1, nrow = d_, ncol = d_)
#' tran_c =  diag(1, nrow = d_, ncol = d_)
#' ini_c = diag(1, nrow = d_, ncol = d_)
#' obs_m = diag(1, nrow = d_, ncol = d_)
#' obs_c = diag(1, nrow = d_, ncol = d_)
#' parameters_ <- list(k = 7, tau = 0.5, kappa = 0.5)
#' obs_p <- list(obs_mean = obs_m, obs_cov = obs_c)
#'
#' #output <- generate_blocks(lag, Time)
#' #output <- generate_blocks_half(lag, Time)
#' #breaks_ <- output[[1]]
#' #psi_index_ <- output[[2]]
#'
#' model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m, tran_cov = tran_c, obs_params = obs_p,
#'  eval_likelihood = evaluate_likelihood_lg, simu_observation = simulate_observation_lg,
#'  parameters = parameters_, dist = 'lg')
#'
#' obs_ <- sample_obs(model, Time, d_) #provided by users
#' #obs_ <- as.matrix(read.csv('C:/Users/ip21972/Downloads/quasi.online.iAPF.package-main/quasi.online.iAPF.package-main/try_data.csv'))
#'
#'a0_ = ini     # Initial state mean
#'P0_ = ini_c    # Initial state covariance
#'Zt_ = obs_m    # Observation matrix (C)
#'Ht_ = tran_c    # Observation noise covariance (R)
#'Gt_ = obs_c   # Process noise covariance (Q)
#'
#' dt_ <- ct_ <- matrix(0, d_, 1)
#' Tt_ <- as.matrix(tran_m)
#' params <- list(dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_ , Zt = Zt_,
#'                Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)
#'
#' filter <- compute_fkf(params, obs_)
#' fkf_obj <- filter[[1]]
#' fks_obj <- filter[[2]]
#'
#'filt_means <- t(fkf_obj$att)
#'filt_covs  <- fkf_obj$Ptt
#'smooth_means <- t(fks_obj$ahatt)
#'smooth_covs <- fks_obj$Vt
#'
#'data_ =  data <- list(obs = obs_)
#'#data <- list(obs = obs_, breaks = breaks_,
#'#psi_index = psi_index_)
#'
#'
#' log_ratio <- vector()
#' log_ratio_rolling <- vector()
#' log_ratio_iapf <- vector()
#' log_ratio_apf <- vector()
#' avg <- matrix(nrow = 1, ncol = Time)
#' fkf_obj_estimates <- 0
#'
#'
#' num_runs <- 1
#' logZ_matrix_rolling <- matrix(NA, nrow = num_runs, ncol = Time)
#'
#' for(i in 1:num_runs){
#' set.seed(i*2)
#'
#' time_info <- system.time({
#' output <- Orc_SMC(lag, data, model, Napf)
#' })
#' print(time_info["elapsed"])
#'
#'
#' logZ_matrix_rolling[i, ] <- output$logZ
#' #fkf_obj_estimates <- output$f_means
#' filter <- compute_fkf(params, obs_[1:Time,,drop = FALSE])
#' fkf_obj <- filter[[1]]
#'
#' log_ratio_rolling[i] <- compute_ratio(logZ_matrix_rolling[i,Time], fkf_obj)
#' print(log_ratio_rolling[i] )
#' }
#'
#'#The Backward smoothing and L1 error
#' H_list <- output$H_forward
#' smooth_particles <- backward_smoothing(H_list)
#' L1_error <- compute_l1_at_t(Time, data, smooth_particles)
#'
#' logZ_matrix_iapf <- matrix(NA, nrow = num_runs, ncol = Time)
#' for(i in 1:num_runs){
#' set.seed(i)
#' output <- iAPF(data, Napf, model)
#'
#' psi_pa <- output$psi
#' logZ_matrix_iapf[i, ] <- output$Z
#' #fkf_obj_estimates <- output$f_means
#' log_ratio_iapf[i] <- compute_ratio(logZ_matrix_iapf[i,Time], fkf_obj)
#' print(log_ratio_iapf[i] )
#' }
#'
#'
#'
#' log_ratio_rolling_vec <- matrix(NA, nrow = num_runs, ncol = Time)
#' for(n in 1:num_runs){
#' for(i in 1:1){
#' filter <- compute_fkf(params, obs_[1:i,, drop = FALSE])
#' fkf_obj <- filter[[1]]
#' log_ratio_rolling_vec[n, i] <- compute_log_ratio(logZ_matrix_rolling[n, i], fkf_obj)
#' }
#' }
#'
#'
#' fks_obj <- filter[[2]]
#'
#'for(i in 1:50){
#' set.seed(i*2)
#' output_apf <- run_bpf(model, data, lag, Napf)
#' logZ <- output_apf[[3]]
#' log_ratio_apf[i] <- compute_log_ratio(logZ, fkf_obj)
#' print(log_ratio_apf[i] )
#'}
#'
#'start_time <- Sys.time()
#'for(i in 1:50){
#' set.seed(i*2)
#' output_apf <- run_bpf(model, data, lag, Napf)
#' logZ <- output_apf[[3]]
#' log_ratio_apf[i] <- compute_log_ratio(logZ, fkf_obj)
#' print(log_ratio_apf[i] )
#'}
#' end_time <- Sys.time()
#' bpf_time <- end_time - start_time
#'
#'start_time <- Sys.time()
#'for(i in 1:1){
#'set.seed(i+2)
#' #run the algorithm
#' output <- run_quasi_online_pf(model, data, Napf, N)
#'
#' X<- output[[1]]
#' w<- output[[2]]
#' logZ <- output[[3]]
#' psi_final <- output[[4]]
#' X_apf <- output[[5]]
#' w_apf <- output[[6]]
#' ancestors <- output[[7]]
#'
#' log_ratio[i] <- compute_log_ratio(logZ, fkf_obj)
#' print(log_ratio[i] )
#' dist <- compute_dKS(X, w, fks_obj)
#'
#' for (t in 1:Time){
#' avg[, t] <- length(unique(X[t,,,drop = FALSE]))/d_
#' }
#'
#' plot(x = c(1:Time), y = avg[1,])
#' }
#' end_time <- Sys.time()
#' iapf_time <- end_time - start_time
#'
#' specific_time = 105
#'
#' #update observations:
#'
#' if(nrow(obs_) < specific_time){
#' data$obs <- rbind(obs_, sample_obs(model, specific_time - nrow(obs_), d_))
#' }
#'
#' output <- generate_blocks(lag, specific_time)
#' data$breaks <- output[[1]]
#' #data$psi_index <- output[[2]]
#'
#' previous_info <- list(previous_time = 100,
#' X_pre = X_apf, w_pre = w_apf, psi_final = psi_final, X = X, w = w)
#' output1 <- run_quasi_online_pf(model, data, Napf, N, previous_info)
#' logZ <- output1[[3]]
#'
#' filter_t <- compute_fkf_fkf_obj(params, data$obs)
#' fkf_obj_t <- filter_t[[1]]
#' log_ratio_t <- compute_log_ratio(logZ, fkf_obj_t)
#' }

