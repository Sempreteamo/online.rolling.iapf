
perform_online_setting <- function(data, specific_time, w_record, X_record, Napf, psi_final){
  breaks_ <- data$breaks
  obs <- data$obs
  psi_index <- data$psi_index
  run_block <- data$run_block
  t = specific_time
  Time <- nrow(obs)
  pre_t <- nrow(psi_final)
  
  if(pre_t < Time){
    psi_final <- rbind(psi_final, matrix(NA, Time - pre_t, dim(psi_final)[2]))
  }
  
  output <- run_iAPF(model, data, Napf, X_record, w_record)
  
  X_record <- output[[1]]
  w_record <- output[[2]]
  psi <- output[[3]]
  
  M = max(breaks_[[1]][run_block[1] - 1], breaks[[2]][run_block[3] - 1])
  

  psi_final[M:specific_time,] <- combine_psi(psi, psi_index[M:specific_time])
  
  output1 <- run_psi_APF(model, list(obs[1:t,], 
                                     breaks[[1]][1], 0, 0), Napf, psi_final, init = FALSE)
  #X <- output1[[1]]
  #w <- output1[[2]]
  logZ <- output1[[3]]
  
  return(list(X = X_record, w = w_record, logZ = logZ, psi_final = psi_final))
  
}

