
perform_online_setting <- function(data, specific_time, X_apf, w_apf, psi_final, logZ){
  breaks <- data$breaks
  obs <- data$obs
  Time <- nrow(obs)
  
  #find and extend the nearest block to t
  nearest <- unlist(breaks)[which.max(unlist(breaks) < specific_time)]
  
  if(nearest %in% breaks[[1]] ){
    index <- 1
  }else{
    index <- 2
  }

  
  data$run_block <- c(index, nearest + 1, specific_time)
  data$past <- list(X_apf[[index]], w_apf[[index]]) 
  
  #Run iAPF to get psi and X, w at final time point
  output <- run_iAPF(model, data, Napf)
  X_apf[[index]] <- output[[1]]
  w_apf[[index]] <- output[[2]]
  psi_pa <- output[[3]]
  
  #psi_pa <- rbind(psi_final, psi_pa)
  
  #Update psi-APF
  output1 <- run_psi_APF(model, list(obs[(nearest + 1):specific_time,], breaks[[1]][1], 0, 0), 
                         Napf, psi_pa, init = FALSE)
  logZ <- logZ + output1[[3]]
  
  return(list(X = X_apf, w = w_apf, logZ = logZ, psi_final = psi_final))
  
}

