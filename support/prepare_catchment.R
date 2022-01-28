prepare_catchment <- function(z, thrA, x_outlet, y_outlet, maxReachLength){
  
  cellsizeX <- (z@extent@xmax - z@extent@xmin)/z@ncols
  cellsizeY <- (z@extent@ymax - z@extent@ymin)/z@nrows
  
  # Pitremove
  system("mpiexec -n 8 pitremove -z support/DTM.asc -fel DTMfel.tif")
  fel=raster("DTMfel.tif")
  #plot(fel)
  
  
  # D8 flow directions
  system("mpiexec -n 8 D8Flowdir -p DTMp.tif -sd8 DTMsd8.tif -fel DTMfel.tif",show.output.on.console=F,invisible=F)
  p=raster("DTMp.tif")
  #plot(p)
  sd8=raster("DTMsd8.tif")
  #plot(sd8)
  
  # Contributing area
  system("mpiexec -n 8 AreaD8 -p DTMp.tif -ad8 DTMad8.tif")
  ad8=raster("DTMad8.tif")
  #plot(log(ad8))
  #zoom(log(ad8))
  
  # Threshold
  system(sprintf("mpiexec -n 8 Threshold -ssa DTMad8.tif -src DTMsrc.tif -thresh %.2f",thrA/cellsizeX/cellsizeY))
  src=raster("DTMsrc.tif")
  #plot(src)
  #zoom(src)
  
  shp.point(x_outlet,y_outlet,"ApproxOutlet")
  
  #makeshape.r("ApproxOutlets")
  
  # Move Outlets
  system("mpiexec -n 8 moveoutletstostreams -p DTMp.tif -src DTMsrc.tif -o ApproxOutlet.shp -om Outlet.shp")
  
  # Contributing area upstream of outlet
  system("mpiexec -n 8 Aread8 -p DTMp.tif -o Outlet.shp -ad8 DTMssa.tif")
  ssa=raster("DTMssa.tif")
  #plot(ssa) 
  
  # new stuff ####
  tmp <- coordinates(ssa)
  X_FD <- tmp[,1]
  Y_FD <- tmp[,2]
  rm(tmp)
  Z_FD <- values(z)
  
  #xllcorner <- z@extent@xmin
  #yllcorner <- z@extent@ymin
  ncols <- z@ncols
  nrows <- z@nrows
  cellsize <- sqrt(cellsizeX*cellsizeY)
  
  cont <- rasterToContour(reclassify(ssa,cbind(NA,0)),levels=0.5)
  
  indexRNNodes <- which(values(ssa)>=thrA/cellsize^2)
  tmp <- xyFromCell(ssa,indexRNNodes)
  X_RN <- tmp[,1]
  Y_RN <- tmp[,2]
  
  tmp <- values(ssa)
  A_RN <- tmp[which(tmp>=thrA/cellsize^2)]
  tmp <- values(fel)
  Z_RN <- tmp[which(values(ssa)>=thrA/cellsize^2)]
  tmp <- values(p)
  Length_RN <- cellsize*(1 + as.numeric(tmp[which(values(ssa)>=thrA/cellsize^2)] %% 2 ==0)*(sqrt(2)-1)) 
  rm(tmp)
  
  ## W, downNode at RN level
  
  nNodes_RN <- length(indexRNNodes)
  W_RN <- spam(0,nNodes_RN,nNodes_RN)
  ind <- matrix(0,nNodes_RN,2)
  downNode_RN <- numeric(nNodes_RN)
  flowDir <- values(p)
  Slope_RN <-  numeric(nNodes_RN)
  
  k <- 1
  for (i in 1:nNodes_RN){
    mov <- neigh(flowDir[indexRNNodes[i]])
    d <- which(indexRNNodes==(indexRNNodes[i]+mov[1]+mov[2]*ncols)) # indices from top-left corner to the right, then next row...
    if (length(d)!=0){
      ind[k, ] <- c(i,d)
      k <- k + 1
      Slope_RN[i] <- (Z_RN[i]-Z_RN[d])/Length_RN[i]
    } 
  }
  ind <- ind[-k, ]
  downNode_RN[ind[,1]] <- ind[,2]
  W_RN[ind] <- 1
  rm(ind)
  Outlet_RN <- which(downNode_RN==0)
  
  # find AG nodes ####
  DegreeIn <- colSums(W_RN)
  DegreeOut <- rowSums(W_RN)
  Confluence <- DegreeIn>1
  Source <- DegreeIn==0
  SourceOrConfluence <- Source|Confluence
  ConfluenceNotOutlet <- Confluence&(downNode_RN!=0)
  ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet
  
  OutletNotChannelHead <- (downNode_RN==0)&(!ChannelHeads)
  IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
  whichNodeAG <- which(IsNodeAG)
  
  nNodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(nNodes_AG)
  RN_to_AG <- numeric(nNodes_RN)
  reachID <- 1
  X_AG <- NaN*numeric(nNodes_AG)
  Y_AG <- NaN*numeric(nNodes_AG)
  Z_AG <- NaN*numeric(nNodes_AG)
  A_AG <- NaN*numeric(nNodes_AG)
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID 
    j <- downNode_RN[i] 
    X_AG[reachID] <- X_RN[i]
    Y_AG[reachID] <- Y_RN[i]
    Z_AG[reachID] <- Z_RN[i]
    A_AG[reachID] <- A_RN[i]
    Length_AG[reachID] <- Length_RN[i]
    tmp_length <- Length_RN[i]
    tmp <- NULL
    j0 <- j
    while (!IsNodeAG[j] && j!=0) {
      tmp <- c(tmp, j)
      tmp_length <-  tmp_length + Length_RN[j]
      j_old <- j
      j <- downNode_RN[j]} 
    
    if (tmp_length > maxReachLength){
      n_splits <- ceiling(tmp_length/maxReachLength)
      new_maxLength <- tmp_length/n_splits
      j <- j0
      while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= new_maxLength) {
        RN_to_AG[j] <- reachID 
        Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
        j_old <- j
        j <- downNode_RN[j]}
      if (Length_AG[reachID] > new_maxLength){
        j <- j_old
        Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
        ChannelHeads[j] <- 1
        whichNodeAG <- c(whichNodeAG,j)}
      
    } else {
      RN_to_AG[tmp] <- reachID
      Length_AG[reachID] <- tmp_length
    }
    
    reachID <- reachID + 1
    whichNodeAG <- whichNodeAG[-1]
  }
  nNodes_AG <- length(X_AG)
  
  # W, downNode at AG level ####
  
  downNode_AG <- numeric(nNodes_AG)
  W_AG <- spam(0,nNodes_AG,nNodes_AG)
  ind <- matrix(0,nNodes_AG,2)
  reachID <- sum(ChannelHeads) + 1
  for (i in 1:nNodes_RN){ 
    if (downNode_RN[i] != 0 && RN_to_AG[downNode_RN[i]] != RN_to_AG[i]) {
      downNode_AG[RN_to_AG[i]] <- RN_to_AG[downNode_RN[i]]
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],downNode_AG[RN_to_AG[i]])
    }
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  Outlet_AG <- RN_to_AG[Outlet_RN]
  
  AG_to_RN <- vector("list", nNodes_AG)
  for(i in 1:nNodes_AG) { # attribute river network pixels to fields of the AG_to_FD list 
    AG_to_RN[[i]] <- which(RN_to_AG==i) 
  }
  
  FD_to_SC <- NA*numeric(length(flowDir))
  SC_to_FD <- vector("list",nNodes_AG)
  FD_to_SC[indexRNNodes] <- RN_to_AG
  for (i in 1:nNodes_AG){
    SC_to_FD[[i]] <- indexRNNodes[which(RN_to_AG==i)]
  }
    
  
  # find FD_to_SC
  drainageArea <- values(ssa)
  indexFDNodes <- which(drainageArea>0)
  ind_head <- indexFDNodes[which(drainageArea[indexFDNodes]==1)]
  
  for (i in 1:length(ind_head)){
    d <- ind_head[i]
    k <- NA; d_new <- d; sub_d <- numeric(0)
    while (is.na(k)){
      k <- FD_to_SC[d_new]
      if (is.na(k)){
        sub_d <- c(sub_d, d_new)
        mov <- neigh(flowDir[d_new])
        d_new <- d_new + mov[1] + mov[2]*ncols # indices from top-left corner to the right, then next row...
      }
    }
    FD_to_SC[sub_d] <- k
    SC_to_FD[[k]] <- c(SC_to_FD[[k]], sub_d)
  }
  
  
  # Upstream_RN : list containing IDs of all reaches upstream of each reach (plus reach itself)
  # Upstream_RN <- vector("list",nNodes_RN)
  # nUpstream_RN <- numeric(nNodes_RN)
  # for (i in 1:nNodes_RN){
  #   UpOneLevel <- which(downNode_RN==i) # find reaches at one level upstream
  #   Upstream_RN[[i]] <- UpOneLevel      # add them to the list
  #   while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
  #     ContinuePath <- UpOneLevel # jump 1 level above
  #     UpOneLevel <- which(downNode_RN %in% ContinuePath) # find reaches at one level upstream
  #     Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
  #   }
  #   Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
  #   nUpstream_RN[i] <- length(Upstream_RN[[i]])
  #   if ((i %% 100)==0){
  #     cat(sprintf("%.2f%% done\n",100*i/nNodes_RN))
  #   }
  # }
  
  # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_AG <- vector("list",nNodes_AG)
  nUpstream_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    UpOneLevel <- which(downNode_AG==i) # find reaches at one level upstream
    Upstream_AG[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(downNode_AG %in% ContinuePath) # find reaches at one level upstream
      Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
    }
    Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
    nUpstream_AG[i] <- length(Upstream_AG[[i]])
  }
  
  
  # calculate Strahler stream order
  StreamOrder_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    j <- order(nUpstream_AG)[i] # index that explores reaches in a downstream direction
    tmp <- which(downNode_AG==j) # set of reaches draining into j
    if (length(tmp)>0){
      IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
      if (IncreaseOrder > 1) {
        StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
      } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
    } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
  }
  
  
  #print('Length and slope at AG level...',quote=FALSE) 
  # Calculate length and slopes of reaches
  #Length_AG <- rep(0,Nnodes_AG)
  Slope_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    #Length_AG[i] <- sum(OCN$FD$leng[AG_to_FD[[i]]])
    Slope_AG[i] <- (Slope_RN[RN_to_AG==i] %*% Length_RN[RN_to_AG==i])/Length_AG[i] # scalar product between vector of slopes and lengths of nodes at RN level belonging to reach i 
  }
  
  
  nNodes_SC <- nNodes_AG
  Z_SC <- numeric(nNodes_SC)
  Alocal_SC <- numeric(nNodes_SC)
  for (i in 1:nNodes_SC) {
    Z_SC[i] <- mean(Z_FD[which(FD_to_SC==i)])
    Alocal_SC[i] <- sum(FD_to_SC==i,na.rm=T)
  }
  
  Areach_AG <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG) {
    Areach_AG[i] <- sum(Alocal_SC[Upstream_AG[[i]]])  
  }
  
  # coordinates of AG nodes considered at the downstream end of the respective edge
  XReach <- numeric(nNodes_AG)
  YReach <- numeric(nNodes_AG)
  ZReach <- numeric(nNodes_AG)
  for (i in 1:nNodes_AG){
    tmp <- AG_to_RN[[i]]
    ind <- which(A_RN[tmp]==max(A_RN[tmp]))
    node <- tmp[ind]
    XReach[i] <- X_RN[node]
    YReach[i] <- Y_RN[node]
    ZReach[i] <- Z_RN[node]
  }
  XReach[Outlet_AG] <- NaN
  YReach[Outlet_AG] <- NaN
  ZReach[Outlet_AG] <- NaN
  
  # # X,Y of subcatchment centroids
  # X_SC <- numeric(nNodes_SC)
  # Y_SC <- numeric(nNodes_SC)
  # for (i in 1:nNodes_SC){
  #   X_SC[i] <- mean(OCN$FD$X[FD_to_SC[[i]]])
  #   Y_SC[i] <- mean(OCN$FD$Y[FD_to_SC[[i]]])
  # }
  
  
  # RN level
  # RN_DownstreamPathLength <- spam(0,nNodes_RN,nNodes_RN)
  # RN_DownstreamPath <- vector("list",nNodes_RN)
  # indices <- matrix(0,1000*nNodes_RN,2)
  # values <- numeric(1000*nNodes_RN)
  # counter <- 1
  # for (i in 1:nNodes_RN){RN_DownstreamPath[[i]] <- vector("list",nNodes_RN)}
  # for (i in 1:nNodes_RN){
  #   Ups <- Upstream_RN[[i]]
  #   for (j in 1:length(Ups)){
  #     k <- Ups[j]
  #     Path <- k
  #     while (k != i) {
  #       k <- downNode_RN[k]
  #       Path <- c(Path, k)
  #     }
  #     RN_DownstreamPath[[Ups[j]]][[i]] <- Path
  #     #AG_DownstreamPathLength[Ups[j],i] <- sum(OCN$AG$leng[Path])
  #     indices[counter, ] <- c(Ups[j],i)
  #     values[counter] <- sum(Length_RN[Path]) 
  #     #if (includeDownstreamNode==FALSE){
  #     values[counter] <- values[counter]  - Length_RN[k]
  #     #}
  #     counter <- counter + 1
  #   }
  #   if ((i %% 100)==0){
  #     cat(sprintf("%.2f%% done\n",100*i/nNodes_RN))
  #   }
  #   
  # }
  # indices <- indices[1:(counter-1), ]
  # values <- values[1:(counter-1)]
  # RN_DownstreamPathLength[indices] <- values
  
  # AG level
  AG_DownstreamPathLength <- spam(0,nNodes_AG,nNodes_AG)
  AG_DownstreamPath <- vector("list",nNodes_AG)
  indices <- matrix(0,1000*nNodes_AG,2)
  values <- numeric(1000*nNodes_AG)
  counter <- 1
  for (i in 1:nNodes_AG){AG_DownstreamPath[[i]] <- vector("list",nNodes_AG)}
  for (i in 1:nNodes_AG){
    Ups <- Upstream_AG[[i]]
    for (j in 1:length(Ups)){
      k <- Ups[j]
      Path <- k
      while (k != i) {
        k <- downNode_AG[k]
        Path <- c(Path, k)
      }
      AG_DownstreamPath[[Ups[j]]][[i]] <- Path
      #AG_DownstreamPathLength[Ups[j],i] <- sum(OCN$AG$leng[Path])
      indices[counter, ] <- c(Ups[j],i)
      values[counter] <- sum(Length_AG[Path]) 
      #if (includeDownstreamNode==FALSE){
      #values[counter] <- values[counter]  - Length_AG[k]
      #}
      counter <- counter + 1
    }
  }
  indices <- indices[1:(counter-1), ]
  values <- values[1:(counter-1)]
  AG_DownstreamPathLength[indices] <- values
  
  
  
  AG_DwnstrLength_unconnected <- spam(0,nNodes_AG,nNodes_AG)
  indices <- matrix(0,nNodes_AG^2,2)
  values <- numeric(nNodes_AG^2)
  counter <- 1
  for (i in 1:nNodes_AG){
    for (j in 1:nNodes_AG){
        if (AG_DownstreamPathLength[i,j]==0 & AG_DownstreamPathLength[j,i]==0){
        k  <- intersect(AG_DownstreamPath[[i]][[Outlet_AG]],AG_DownstreamPath[[j]][[Outlet_AG]])[1]
        #AG_DwnstrLength_unconnected[i,j] <- AG_DownstreamPathLength[i,k]
        indices[counter,] <- c(i,j)
        values[counter] <- AG_DownstreamPathLength[i,k]
        counter <- counter + 1
        #CountPaths[k] <- CountPaths[k]  + 1 
        }
    }
    cat(sprintf(" %d \n",i))
  }
  indices <- indices[1:(counter-1), ]
  values <- values[1:(counter-1)]
  AG_DwnstrLength_unconnected[indices] <- values
  
  
  catch <- vector("list")
  
  catch$FD[["X"]] <- X_FD
  catch$FD[["Y"]] <- Y_FD
  catch$FD[["Z"]] <- Z_FD
  catch$FD[["toSC"]] <- FD_to_SC
  
  catch$RN[["A"]] <- A_RN*cellsize^2
  catch$RN[["downNode"]] <- downNode_RN
  # catch$RN[["Upstream"]] <- Upstream_RN
  # catch$AG[["downstreamPath"]] <- AG_DownstreamPath
  # catch$AG[["downstreamPathLength"]] <- AG_DownstreamPathLength
  catch$RN[["X"]] <- X_RN
  catch$RN[["Y"]] <- Y_RN
  catch$RN[["toAGReach"]] <- RN_to_AG
  catch$RN[["nNodes"]] <- nNodes_RN
  
  catch$AG[["A"]] <- A_AG*cellsize^2
  catch$AG[["AReach"]] <- Areach_AG*cellsize^2
  catch$AG[["downNode"]] <- downNode_AG
  catch$AG[["downstreamPath"]] <- AG_DownstreamPath
  catch$AG[["downstreamPathLength"]] <- AG_DownstreamPathLength
  catch$AG[["downstreamLengthUnconnected"]] <- AG_DwnstrLength_unconnected
  catch$AG[["leng"]] <- Length_AG
  catch$AG[["outlet"]] <- Outlet_AG
  catch$AG[["slope"]] <- Slope_AG
  catch$AG[["streamOrder"]] <- StreamOrder_AG
  catch$AG[["Upstream"]] <- Upstream_AG
  catch$AG[["X"]] <- X_AG
  catch$AG[["Y"]] <- Y_AG
  catch$AG[["Z"]] <- Z_AG
  catch$AG[["W"]] <- W_AG
  catch$AG[["nNodes"]] <- nNodes_AG
  pl <- initial_permutation(downNode_AG)
  catch$AG[["perm"]] <- pl$perm
  
  catch$SC[["toFD"]] <- SC_to_FD
  catch$SC[["A"]] <- Alocal_SC*cellsize^2
  
  catch$CM[["XContour"]] <-  cont@lines[[1]]@Lines[[1]]@coords[,1]
  catch$CM[["YContour"]] <-  cont@lines[[1]]@Lines[[1]]@coords[,2]
  catch$CM[["A"]] <- max(A_RN)*cellsize^2
  
  catch[["cellsize"]] <- cellsize
  catch[["thrA"]] <- thrA
  
  
  invisible(catch)
}

# AUXILIARY FUNCTIONS ####

# function to create point shapefile given coordinates
shp.point <- function(x,y,sname="shape"){
  n <- length(x)
  dd <- data.frame(Id=1:n,X=x,Y=y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}

initial_permutation <- function(DownNode){
  
  Outlet <- which(DownNode==0)
  NodesToExplore <- Outlet # start from outlets
  reverse_perm <- numeric(length(DownNode)) # build permutation vector from outlets to headwaters, then flip it
  
  k <- 0
  while (length(NodesToExplore)>0){ # continue until all the network has been explored
    k <- k + 1
    node <- NodesToExplore[1] # explore a node
    reverse_perm[k] <- node # assign position in the permutation vector
    NodesToExplore <- NodesToExplore[-1] # remove explored node
    UpNodes <- which(DownNode==node) # find nodes upstream of node
    while (length(UpNodes)>0){ # continue upstream until a headwater is found
      k <- k + 1
      node <- UpNodes[1] # explore first upstream node
      reverse_perm[k] <- node
      if (length(UpNodes)>1){ # if there is a bifurcation upstream, add the other upstream connections at the top of NodesToExplore
        NodesToExplore <- c(UpNodes[2:length(UpNodes)],NodesToExplore)
      }
      UpNodes <- which(DownNode==node)
    }
  }
  
  perm <- reverse_perm[length(DownNode):1] # flip permutation
  
  OutList = list(perm=perm,noDAG=0)
  
  invisible(OutList)
}
