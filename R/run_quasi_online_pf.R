#' Function to run the whole quasi-online-iAPF processes and generate final results
#'
#' @param model Model information
#' @param data Observations
#' @param Napf Number of particles used in the smoothing process
#' @param N Number of particles used in the filtering process
#' @param previous_info Transfer information of the old results when the function is used online, 
#'
#' @return A list contains
#' X is the particles generated within the psi-APF
#' w is the weights of particles generated within the psi-APF
#' logZ is the estimated normalising constant
#' psi_final is psi function generated
#' X_apf is the particles generated within the iAPF
#' w_apf is the weights of particles generated within the iAPF
#'
#' @export
#'
run_quasi_online_pf <- function(model, data, Napf, N, previous_info = NULL){
  breaks <- data$breaks
  psi_index <- data$psi_index
  obs <- data$obs
  Time <- nrow(obs)
  d = ncol(obs)
  #psi_final[previous_time, ] <- data$
  block <- list(0, 0)
  ancestors <- matrix(NA, Time, N)
  
  if(is.null(previous_info)){
    previous_time <- 0
    X_apf <- list(vector(), vector())
    w_apf <- list(vector(), vector()) 
    psi_final <- matrix(NA, Time, 2*d)
    psi_l = 1
    X <- array(NA, c(Time, N, d))
    w <- matrix(NA, Time, N)
  }else{
    previous_time <- previous_info[[1]]
    X_apf <- previous_info[[2]]
    w_apf <- previous_info[[3]]
    psi_final <- rbind(previous_info[[4]], matrix(NA, Time - previous_time, 2*d))
    X <- previous_info[[5]]
    w <- previous_info[[6]]
    psi_l <- previous_time 
  }
  
  logZ = 0
  #b = 0
  
  #step 1: Decide whether t is the upper boundary of a block
  for(t in (previous_time + 1):Time){
    print(t)
    
    if(t %in% unlist(breaks)){

      if(previous_time != 0){
       indices <- which(sapply(breaks, function(b) 
          max(unlist(breaks)[unlist(breaks) < previous_time] ) %in% b)) 
       psi_index <- c(psi_index, rep(indices, Time - previous_time))
       
      }else{
        indices <- which(sapply(breaks, function(b) t %in% b)) 
      }
      
      for (index in indices){
        
        if(t == Time){
          b_s <- if (length(breaks[[index]]) > 1) 
            breaks[[index]][length(breaks[[index]]) - 1] + 1 else 1

          
        }else{
          b_s <- max(breaks[[index]][block[[index]]], 1)
        }
        
        block[[index]] <- block[[index]] + 1

        data$run_block <- c(index, b_s, t) #which layer, which block iAPF runs
        data$past <- list(X_apf[[index]], w_apf[[index]]) #initialize distribution
        
        
        #Step 2: Run iAPF to get psi and X, w at final time point
        #cat('no')
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
                          which(is.na(psi_final))[1] - 2) #psi updated to where
        
       
        if( psi_u > 0){
          print(c(psi_l, psi_u))
          #output1 <- APF_with_resampling_logZ(model, obs, w, X, Napf, psi_final, psi_l , psi_u, logZ, ancestors)
          if(psi_l != 1){
            if(psi_u != Time){
              output1 <- run_psi_APF(model, list(obs[psi_l:psi_u,], c(psi_l,psi_u - 1), w[psi_l - 1,], as.matrix(X[psi_l - 1,,])), 
                                     Napf, psi_final[psi_l:(psi_u+1),], init = FALSE, jump_ini = TRUE, jump_last = TRUE)
            }else{
              output1 <- run_psi_APF(model, list(obs[psi_l:psi_u,], c(psi_l,psi_u), w[psi_l - 1,], as.matrix(X[psi_l - 1,,])), 
                                     Napf, psi_final[psi_l:psi_u,], init = FALSE, jump_ini = TRUE)
            }
            
          }else{
            output1 <- run_psi_APF(model, list(obs[psi_l:psi_u,], c(psi_l,psi_u),0,0), 
                                   Napf, psi_final, init = FALSE, jump_last = TRUE)
          }
          
          if(is.null(previous_info)){
            X[psi_l:psi_u,,] <- output1[[1]]
            w[psi_l:psi_u,] <- output1[[2]]
          }
          
          logZ <- output1[[3]]+logZ
          print(logZ)
          
          #if(psi_u != Time){
           # w[psi_u,] <- w[psi_u,] - normalise_weights_in_log_space(w[psi_u,])[[2]]  #+ logZ
          #ancestors[psi_l:psi_u,] <- output1[[4]][psi_l:psi_u,]
          #}
        
        psi_l = max(psi_u + 1, 1)
        
      }
      
    }
    

  }
}
  return(list(X = X, w = w, logZ = logZ, psi_final = psi_final, X_apf = X_apf, w_apf = w_apf))
}



