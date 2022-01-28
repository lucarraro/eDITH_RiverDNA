extract_unconnected_pairs <- function(catch, alpha=NULL){
  
  # if alpha is provided, only pairs where at least one site has alpha>0 are picked
 
  vv <- which(catch$AG$A <= median(catch$AG$A))
  ww <- as.matrix(catch$AG$downstreamLengthUnconnected[vv,vv])
  ww[upper.tri(ww)] <- 0
  zz <- which(ww>0,arr.ind = T)
  yy <- matrix(c(vv[zz[,1]],vv[zz[,2]]),ncol=2) # all pairs of flow-unconnected sites with A <= median
  yy_tmp <- yy
  
  if (!is.null(alpha)){
    alpha_mat <- cbind(alpha[yy[,1]], alpha[yy[,2]])
    foo <- which(rowSums(alpha_mat)==0)
    if (length(foo)>0){
    yy_tmp <- yy_tmp[-foo,] # remove pairs where both sites have alpha=0
  }}
  
  HeadPairs <- matrix(0,1000,2)
  i <- 1
  while (length(yy_tmp)>0){
    if (length(yy_tmp)>2){
      sampled <- yy_tmp[sample(1:(length(yy_tmp)/2),1),]
      rowsToRemove <- unique(which(yy_tmp[,1]==sampled[1] | yy_tmp[,2]==sampled[1] | yy_tmp[,1]==sampled[2] | yy_tmp[,2]==sampled[2]))
      yy_tmp <- yy_tmp[-rowsToRemove,]
    } else {
      sampled <- yy_tmp
      yy_tmp <- numeric(0)}
    HeadPairs[i,] <- sampled
    i <- i+1
  }
  HeadPairs <- HeadPairs[1:(i-1),]
  headIndices <- HeadPairs[,1] + (HeadPairs[,2] - 1)*catch$AG$nNodes
  
  vv <- which(catch$AG$A > median(catch$AG$A))
  ww <- as.matrix(catch$AG$downstreamLengthUnconnected[vv,vv])
  ww[upper.tri(ww)] <- 0
  zz <- which(ww>0,arr.ind = T)
  yy <- matrix(c(vv[zz[,1]],vv[zz[,2]]),ncol=2) # all pairs of flow-unconnected sites with A > median
  yy_tmp <- yy
  
  if (!is.null(alpha)){
    alpha_mat <- cbind(alpha[yy[,1]], alpha[yy[,2]])
    foo <- which(rowSums(alpha_mat)==0)
    if (length(foo)>0){
      yy_tmp <- yy_tmp[-foo,] # remove pairs where both sites have alpha=0
    }}
  
  DownPairs <- matrix(0,1000,2)
  i <- 1
  while (length(yy_tmp)>0){
    if (length(yy_tmp)>2){
      sampled <- yy_tmp[sample(1:(length(yy_tmp)/2),1),]
      rowsToRemove <- unique(which(yy_tmp[,1]==sampled[1] | yy_tmp[,2]==sampled[1] | yy_tmp[,1]==sampled[2] | yy_tmp[,2]==sampled[2]))
      yy_tmp <- yy_tmp[-rowsToRemove,]
    } else {
      sampled <- yy_tmp
      yy_tmp <- numeric(0)}
    DownPairs[i,] <- sampled
    i <- i+1
  }
  DownPairs <- DownPairs[1:(i-1),]
  downIndices <- DownPairs[,1] + (DownPairs[,2] - 1)*catch$AG$nNodes
  
  ls <- list(headIndices=headIndices, downIndices=downIndices)
  return(ls)
}