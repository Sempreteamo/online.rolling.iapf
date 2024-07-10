#ess resampling

#logW_un = w[t-1,]

compute_ESS_log <- function(logW_un){
  #print(log_un)
    logW <- normalise_weights_in_log_space(logW_un)[[1]]
    ess <- 1/sum(logW^2)
    return(ess)  
}
