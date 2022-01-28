neigh <- function(dir) {
  mov <- c(0,0)
  switch(dir,
         {mov[1] <- 1; mov[2] <- 0},   # case 1 (E)
         {mov[1] <- 1; mov[2] <- -1},   # case 2 (NE) 
         {mov[1] <- 0; mov[2] <- -1},   # case 3 (N)
         {mov[1] <- -1; mov[2] <- -1},  # case 4 (NW)
         {mov[1] <- -1; mov[2] <- 0},  # case 5 (W)
         {mov[1] <- -1; mov[2] <- 1}, # case 6 (SW)
         {mov[1] <- 0; mov[2] <- 1},  # case 7 (S)
         {mov[1] <- 1; mov[2] <- 1})  # case 8 (SE)
  return(mov)
}



