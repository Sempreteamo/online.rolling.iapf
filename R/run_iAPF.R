run_iAPF <- function(model, data, Napf){
  breaks <- data$breaks
  psi_index <- data$psi_index
  obs <- data$obs
  k <- model$k
  Time <- nrow(obs)
  d = ncol(obs)
  psi_final <- list()

  for(index in 1:2){
    combined <- array(NA, dim = c(Napf, Time, d))
    #cat('index=', index)
    psi_pa1 = NULL
    #ancestors1 = NULL
    #w_apf1 = NULL

    for(b in 2:length(breaks[[index]])){
      #print(b)
      l = 1
      Z_apf <- vector()
      N[l] = Napf
      if(b == 2){
        output <- run_psi_APF(model, list(obs[breaks[[index]][(b-1)]:(breaks[[index]][b]-1),],
                   breaks[[index]][(b-1):b], 0, 0), N[l], psi_pa = 0, init = TRUE)
        #print('done b=2')
      }else{
        output <- run_psi_APF(model, list(obs[breaks[[index]][(b-1)]:(breaks[[index]][b]-1),],
                  breaks[[index]][(b-1):b], w_apf[nrow(w_apf),], X_apf[nrow(X_apf),,]), N[l], psi_pa = 0, init = TRUE)
      }

      X_apf <- output[[1]]
      w_apf <- output[[2]]
      Z_apf[l] <- output[[3]]
      ancestors <- output[[4]]

      while(TRUE){

        output <- list()

        if(l != 1){
          #print(l)
          #generate filtering particles X_apf for psi the next iteration
          #APF outputs filtering X_apf for the next psi, and smoothing X_apf_s
          #for the final calculation

          output <- run_psi_APF(model, list(obs[breaks[[index]][(b-1)]:(breaks[[index]][b]-1),],
                                            breaks[[index]][(b-1):b], w_apf[nrow(w_apf),], X_apf[nrow(X_apf),,]),
                                N[l], psi_pa, init = FALSE)
          X_apf <- output[[1]]
          w_apf <- output[[2]]
          Z_apf[l] <- output[[3]]
          ancestors <- output[[4]]

        }

        #to speed up the algorithm, I just fix the number of iterations to be k.
        #Here k = 5

        if(l <= k ){
          #cat('l=',l)
          #receive filtering particles X_apf for psi
          psi_pa <- learn_psi(X_apf, obs[breaks[[index]][(b-1)]:(breaks[[index]][b]-1),], model)

          if(l > k & N[max(l-k,1)] == N[l] & is.unsorted(Z_apf[max(l-k,1):l])){
            N[l+1] <- 2*N[l]

          }else{
            N[l+1] <- N[l]
          }

          l <- l+1
          #print('finish l')
        }else break
      }

      psi_pa1 <- rbind(psi_pa1, psi_pa)
      #ancestors1 <- rbind(ancestors1, ancestors)
      #w_apf1 <- rbind(w_apf1, w_apf)
      #combined[breaks[[index]][(b-1)]:(breaks[[index]][b]-1),,] <- X_apf

    }

    psi_final[[index]] <- psi_pa1
  }
  #output psi
  return(list(X_apf, w_apf, psi_final, Z_apf[l], ancestors))
}
