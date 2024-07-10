
psi1 <- matrix(1:10, 5, 2)
psi2 <- matrix(11:20, 5, 2)
psi <- list(psi1, psi2)
index <- sample(1:2, 5, replace = TRUE)
test_matrix <- rbind(psi[[index[1]]][1,], psi[[index[2]]][2,],
                        psi[[index[3]]][3,], psi[[index[4]]][4,],
                        psi[[index[5]]][5,])


test_combine_psi <- function(psi, index, test_matrix){
  output <- combine_psi(psi, index)
  
  if (!identical(output, test_matrix)) {
    stop("Test failed: Incorrect output from optimization function.")
  }
  
}

test_combine_psi(psi, index, test_matrix)