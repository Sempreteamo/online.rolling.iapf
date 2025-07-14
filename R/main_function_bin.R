#' @examples
#' \dontrun{
#' The following parameters are provided by users
#' library(mvnfast)
#' library(FKF)
#' library(quasi.online.iAPF)
#' Napf = N = 1000
#' lag = 64
#' Time = 100
#' M = 50
#' d_ = 1
#'
#' parameters_ <- list(k = 5, tau = 0.5, kappa = 0.5)
#'
#' ini <- rep(0, d_)
#' ini_c = diag(1, nrow = d_, ncol = d_)
#' tran_m =  diag(0.99, nrow = d_, ncol = d_)
#' tran_c =  diag(0.11, nrow = d_, ncol = d_)
#'
#' obs_p <- M
#'
#'
#' model <- list(ini_mu = ini, ini_cov = ini_c, tran_mu = tran_m,
#' tran_cov = tran_c, obs_params = obs_p,
#' eval_likelihood = evaluate_likelihood_bin, simu_observation = simulate_observation_bin,
#' parameters = parameters_)
#'
#' set.seed(123)
#' obs_ <- sample_obs(model, Time, d_)
#' data <- list(obs = obs_)
#'
#'
#' logZ <- vector()
#'
#' num_runs <- 1
#'
#'for(i in 1:num_runs){
#'set.seed(i*2)
#' #run the algorithm
#' output <- Orc_SMC(lag, data, model, Napf)
#' #output <- run_bpf(model, data, lag, Napf)
#'
#' logZ <- output$logZ
#' #avg <- output[[5]]
#'
#' #log_ratio[i] <- compute_log_ratio(logZ, filtering)
#' print(logZ)
#' #dist <- compute_dKS(X, w, smoothing)
#'
#' #plot(x = c(1:Time), y = avg[1,])
#' 
#' ess_l64 <- output$ess_history
#' plot(ess)
#' }
#' }

