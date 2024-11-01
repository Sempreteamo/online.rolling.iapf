
perform_online_setting <- function(specific_time, w_record, X_record, Napf){
  breaks_ <- data$breaks
  obs <- data$obs
  t = specific_time
  nearest_start_time <-  max(unlist(breaks_)[unlist(breaks_) <= t])

  #find the nearest start time
  Z_apf <- vector()
  l = 1
  k = 5

  output1 <- run_psi_APF(model, 
                         list(obs[nearest_start_time:t,], nearest_start_time:t, w_record[nearest_start_time-1,],
                              X_record[nearest_start_time-1,,]), 
                         Napf, psi_pa = 0, init = TRUE)
  X_apf <- output1[[1]]
  w_apf <- output1[[2]]
  Z_apf[l] <- output1[[3]]
  log_likelihoods <- output1[[6]]
  while(TRUE){
    
    output <- list()
    
    if(l != 1){
      #print(l)
      #generate filtering particles X_apf for psi the next iteration
      #APF outputs filtering X_apf for the next psi, and smoothing X_apf_s
      #for the final calculation
      #print(psi_pa)
      output <- run_psi_APF(model, list(obs[nearest_start_time:t,],
                                        nearest_start_time:t, w_apf[nrow(w_apf),], X_apf[nrow(X_apf),,]),
                            Napf, psi_pa, init = FALSE)
      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf[l] <- output[[3]]
      log_likelihoods <- output[[6]]
      
    }
    
    #to speed up the algorithm, I just fix the number of iterations to be k.
    #Here k = 5
    
    if(l <= k ){
      #cat('l=',l)
      #receive filtering particles X_apf for psi
      psi_pa <- learn_psi(X_apf, obs[nearest_start_time:t,],
                          model, log_likelihoods)
      
      if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z_apf[max(l-k,1):l])){
        N[l+1] <- 2*N[l]
        
      }else{
        N[l+1] <- N[l]
      }
      
      l <- l+1
      #print('finish l')
    }else break
  }
 
  psi_pa_final <- rbind(psi_final1[1:(nearest_start_time-1),],  psi_pa)
  
  output1 <- run_psi_APF(model, list(obs[1:t,], 
                                     breaks[[1]][1], 0, 0), Napf, psi_pa_final, init = FALSE)
  X <- output1[[1]]
  w <- output1[[2]]
  logZ <- output1[[3]]
  
  return(X = X, w = w, logZ = logZ, psi_final = psi_pa_final)
  
}

