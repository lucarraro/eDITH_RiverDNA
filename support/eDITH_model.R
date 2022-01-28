# fast eDITH implementation
eDITH_model <- function(p,tau,width,leng,velocity,Q,nNodes,perm,downNode){
  Conc <- p*width*leng*exp(-leng/tau/velocity)/Q # conc from local eDNA production
  for (i in 1:(nNodes-1)){
    k <- perm[i] # perm allows exploring nodes in the downstream direction
    j <- downNode[k]
    Conc[j] <- Conc[j] + Conc[k]*Q[k]*exp(-leng[j]/velocity[j]/tau)/Q[j]
  }
  return(Conc)
}

