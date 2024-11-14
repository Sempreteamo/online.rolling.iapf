#' @examples
#' \dontrun{
#' The following parameters are provided by users
#' library(mvnfast)
#' library(FKF)
#' Napf = N = 200
#' lag = 10
#' Time = 51
#' d_ = 1
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
#' ini_c = tran_c = obs_m = obs_c = diag(1, nrow = d_, ncol = d_)
#' parameters_ <- list(k = 5, tau = 0.5, kappa = 0.5)
#' obs_p <- list(obs_mean = obs_m, obs_cov = obs_c)
#'
#' output <- generate_blocks(lag, Time)
#' breaks_ <- output[[1]]
#' psi_index_ <- output[[2]]
#'
#' model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m, tran_cov = tran_c, obs_params = obs_p,
#'  eval_likelihood = evaluate_likelihood_lg, simu_observation = simulate_observation_lg,
#'  parameters = parameters_, dist = 'lg')
#'
#' obs_ <- sample_obs(model, Time, d_) #provided by users
#'
#' dt_ <- ct_ <- matrix(0, d_, 1)
#' Tt_ <- as.matrix(tran_m)
#' P0_ <- Zt_ <- Ht_ <- Gt_ <- diag(1, d_, d_)
#' a0_ <- rep(0, d_)
#' params <- list(dt = dt_, ct = ct_, Tt = Tt_, P0 = P0_, Zt = Zt_,
#'                Ht = Ht_, Gt = Gt_, a0 = a0_, d = d_)
#'
#' filter <- compute_fkf_filtering(params, obs_)
#' filtering <- filter[[1]]
#' smoothing <- filter[[2]]
#'
#' data <- list(obs = obs_, breaks = breaks_, 
#' psi_index = psi_index_)
#'
#' kalman <- list(fkf.obj = filtering, fks.obj  = smoothing ) #provided by users
#'
#' log_ratio <- vector()
#'
#'for(i in 1:1){
#'set.seed(i*2)
#' #run the algorithm
#' output <- run_quasi_online_pf(model, data, lag, Napf, N)
#' X<- output[[1]]
#' w<- output[[2]]
#' logZ <- output[[3]]
#' psi <- output[[4]]
#' #avg <- output[[5]]
#'
#' log_ratio[i] <- compute_log_ratio(logZ, filtering)
#' print(log_ratio[i] )
#' #dist <- compute_dKS(X, w, smoothing)
#'
#' #plot(x = c(1:Time), y = avg[1,])
#' }
#' 
#' specific_time = 60
#' 
#' #update observations:
#' 
#' if(nrow(obs_) < specific_time){
#' obs_ <- rbind(obs_, sample_obs(model, specific_time - nrow(obs_), d_)) 
#' 
#' #find and extend the nearest block to t
#' combined_values <- combined_values[-c(length(breaks[[1]]), length(breaks[[1]]) + length(breaks[[2]]))]
#' 
#' nearest <- max(combined_values[combined_values < specific_time])
#' }
#' 
#' 
#' output1 <- position_in_list(breaks_, nearest, combined_values)
#' index <- output1[[1]]
#' position <- output1[[2]] + 1
#' 
#' #update data
#' data$breaks[[index]][position] <- specific_time + 1
#' data$run_block <- c(index, position)
#' data$past <- list(X[nearest - 1,,], w[nearest - 1,])
#' data$obs <- obs_
#' 
#' #psi updated at times min(\Delta_{k, t}^m) ... t-1
#' output_t <- run_iAPF(model, data, Napf)
#' psi_update <- output_t[[3]] 
#' 
#' psi_t <- rbind(psi[1:nearest - 1,], psi_update)
#' 
#' #run psi-APF over times min(\Delta_{k, t}^m) ... t-1
#' output2 <- run_psi_APF(model, list(obs, breaks[[1]][1], 0, 0), Napf, psi_t, init = FALSE)
#' X <- output2[[1]]
#' w <- output2[[2]]
#' logZ_update <- output2[[3]]
#' 

#' filter_t <- compute_fkf_filtering(params, obs_[1:specific_time,])
#' filtering_t <- filter_t[[1]]
#' log_ratio <- compute_log_ratio(logZ_update, filtering_t)
#' }

