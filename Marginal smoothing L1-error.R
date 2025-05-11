#for each time t, coor j, L1 error

compute_l1_at_t <- function(t0, d_) {
  sum = 0
  for(j in 1:d_){
    μ_true    <- smooth_means[t0, j]
    σ2_true   <- smooth_covs[j, j, t0]
    xs        <- seq(μ_true - 4*sqrt(σ2_true),
                     μ_true + 4*sqrt(σ2_true),
                     length.out = 200)
    #cdf of the real kalman distribution
    cdf_true <- pnorm(xs, mean = μ_true, sd = sqrt(σ2_true))
    
    #empirical cdf
    particles_k <- smooth_particles[, t0, j]
    cdf_emp <- sapply(xs, function(x) mean(particles_k <= x))
    
    #L1 distance between cdfs
    dx <- xs[2] - xs[1]
    l1_k <- sum(abs(cdf_emp - cdf_true)) * dx
    
    sum = sum + l1_k
  }
  return(sum/d_)
}

L1_time <- numeric(Time)
for(t in 1:Time){
  L1_time[t] <- compute_l1_at_t(t, d_)
}

df <- data.frame(time = 1:T, L1 = L1_time)
ggplot(df, aes(x = time, y = L1)) +
  geom_line() +
  labs(x = "Time t", y = "Average marginal L1 error",
       title = "Marginal smoothing L1-error over time for N500T1000, d=2") +
  theme_minimal()
