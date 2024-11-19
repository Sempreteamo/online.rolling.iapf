#' Function to run the whole iAPF processes and generate final results
#'
#' @param model Model information
#' @param data Observations
#' @param lag Length of time intervals
#' @param Napf Number of particles used in the smoothing process
#' @param N Number of particles used in the filtering process
#'
#' @return A list contains
#' X is the particles
#' w is the weights of particles
#' logZ is the estimated normalising constant
#' Xs is the smoothing particles
#' log_ratio is the log_ratio between estimates and real normalising constant
#' dist is the distance between real estimates and smoothing distribution
#'
#' @export
#'
run_quasi_online_pf <- function(model, data, lag, Napf, N){
  breaks <- data$breaks
  psi_index <- data$psi_index
  obs <- data$obs
  Time <- nrow(obs)
  Xs <- array(NA, dim = c(nrow(obs), N, d))
  d = ncol(obs)
  psi_final <- matrix(NA, Time, 2*d)
  
  X_apf <- list(NULL, NULL)
  w_apf <- list(NULL, NULL)
  
  logZ = 0
  psi_l = 1
  
  #step 1: Decide whether t is the upper boundary of a block
  for(t in 1:Time){
    
    if(t %in% unlist(breaks)){

      indices <- which(sapply(breaks, function(b) t %in% b)) 
      
      for (index in indices) {
        
        if(t == Time){
          b_s <- if (length(breaks[[index]]) > 1) (breaks[[index]][-2] + 1) else 1
          
        }else{
          b_s <- max((t - lag + 1), 1)
        }

        data$run_block <- c(index, b_s, t) #which layer, which block iAPF runs
        data$past <- list(X_apf[[index]], w_apf[[index]]) #initialize distribution
        
        #Step 2: Run iAPF to get psi and X, w at final time point
        output <- run_iAPF(model, data, Napf)
        X_apf[[index]] <- output[[1]]
        w_apf[[index]] <- output[[2]]
        psi_pa <- output[[3]]
        
        #Step 3: Update psi
        a = 1
        for(i in b_s:t){
          
          if(psi_index[i] == index){
            psi_final[i,] <- psi_pa[a,]
            
          }
          a = a + 1
        }
        
        #Step 4: Update psi-APF
        
        psi_u <- ifelse(is.na(which(is.na(psi_final))[1] - 1), Time, 
               which(is.na(psi_final))[1] - 1) #psi updated to where
        
        if( psi_u > 0){
          output1 <- run_psi_APF(model, list(obs, breaks[[1]][1], 0, 0), 
                                 Napf, psi_final, init = FALSE)
          logZ <- logZ + output1[[3]]
          
        }
        
        psi_l = psi_u + 1
        
      }
      
    }

  }

  return(list(X = X, w = w, logZ = logZ, psi_final = psi_final1))
}



