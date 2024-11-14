#find which block we will run iAPF within
position_in_list <- function(breaks_, nearest, combined_values) {
  if (nearest %in% breaks_[[1]]) {
    list_number <- 1
    position_in_list <- which(breaks_[[1]] == nearest)
    } else {
      list_number <- 2
      position_in_list <- which(breaks_[[2]] == nearest)
      }
  return(c(list_number, position_in_list))
  }