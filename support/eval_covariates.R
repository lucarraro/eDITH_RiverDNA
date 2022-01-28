eval_covariates <- function(catch, geology_shp, 
                            landcover_shp,
                            covariates_morphology,
                            covariates_geol,
                            covariates_landcover){
  
  ### read geology raster
  cat("Read geology shapefile \n")
  geol <- shapefile(geology_shp)
  cat("Rasterize geology \n")
  geol$GEOL_F <- as.factor(geol$GEOL_F)
  geol_raster <- rasterize(geol,DEM,field="GEOL_F")
  
  geol <- values(geol_raster)
  geol[is.na(catch$FD$toSC)] <- NA
  geol[geol==1 | geol==11 | geol==13] <- 1001                                   # 1. Alluvial
  geol[geol==138 | geol==152 | geol==153 | geol==156 | geol==157] <- 1002       # 2. Alpine
  geol[geol==10] <- 1003                                                        # 3. Loess
  geol[geol==20 | geol==21 | geol==22 | geol==23 | geol==24 | geol==25] <- 1004 # 4. Molasses
  geol[geol==12] <- 1005                                                        # 5. Moraines
  geol[geol==2] <- 1006                                                         # 6. Peat
  geol[geol==3 | geol==4 | geol==5] <- 1007                                     # 7. Scree
  geol[geol==14] <- 1008                                                        # 8. Water
  geol <- geol - 1000
  
  
  ### read landcover raster
  landcover <- shapefile(landcover_shp)
  landcover$objval <- as.factor(landcover$objval)
  landcover_raster <- rasterize(landcover,DEM,field="objval")
  
  landcover <- values(landcover_raster)
  landcover[is.na(catch$FD$toSC)] <- NA
  landcover[landcover==9] <- 1001                 # 1. Forest
  landcover[landcover==4] <- 1002                 # 2. Lake
  landcover[landcover==3] <- 1003                 # 3. Orchard  
  landcover[landcover==1 | landcover==2] <- 1004  # 4. Rock
  landcover[landcover==8] <- 1005                 # 5. Swamp
  landcover[landcover==5] <- 1006                 # 6. Urban area
  landcover <- landcover - 1000
  
  
  covariates <- data.frame(matrix(nrow=catch$AG$nNodes,ncol=0))
  
  if ("DrainageArea" %in% covariates_morphology){
    covariates[["DrainageArea"]] <- catch$AG$A
  }
  
  if ("StreamOrder" %in% covariates_morphology){
    covariates[["StreamOrder"]] <- catch$AG$streamOrder
  }
  
  if ("SlopeLoc" %in% covariates_morphology){
    covariates[["SlopeLoc"]] <- catch$AG$slope
  }
  
  if ("MeanSlopeUps" %in% covariates_morphology){
    meanSlopeUps <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      tmp <- catch$AG$Upstream[[i]]
      meanSlopeUps[i] <- sum(catch$AG$slope[tmp]*catch$AG$leng[tmp])/sum(catch$AG$leng[tmp])
    }
    covariates[["MeanSlopeUps"]] <- meanSlopeUps
  }
  
  if ("MeanElevLoc" %in% covariates_morphology){
    meanElevLoc <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      meanElevLoc[i] <- mean(catch$FD$Z[catch$SC$toFD[[i]]])
    }
    
    covariates[["MeanElevLoc"]] <- meanElevLoc
  }
  
  if ("MeanElevUps" %in% covariates_morphology){
    meanElevUps <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      tmp <- catch$AG$Upstream[[i]]
      meanElevUps[i] <- sum(meanElevLoc[tmp]*catch$SC$A[tmp])/sum(catch$SC$A[tmp])
    }
    covariates[["MeanElevUps"]] <- meanElevUps
  }
  
  for (cov in 1:length(covariates_geol)){
    # Local covariate
    vec <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      vec[i] <- sum(geol[catch$SC$toFD[[i]]]==cov)/length(catch$SC$toFD[[i]])
    }
    #covariates[[covariates_geol[cov]]] <- vec
    # Upstream covariate
    vecUps <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      tmp <- catch$AG$Upstream[[i]]
      vecUps[i] <- sum(vec[tmp]*catch$SC$A[tmp])/sum(catch$SC$A[tmp])
    }
    covariates[[paste(covariates_geol[cov],"Ups",sep="")]] <- vecUps
  }
  covariates$MolassesUps <- NULL # remove the most abundant class because it's correlated to the others
  
  for (cov in 1:length(covariates_landcover)){
    # Local covariate
    vec <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      vec[i] <- sum(landcover[catch$SC$toFD[[i]]]==cov,na.rm=T)/length(catch$SC$toFD[[i]])
    }
    covariates[[covariates_landcover[cov]]] <- vec
    # Upstream covariate
    vecUps <- numeric(catch$AG$nNodes)
    for (i in 1:catch$AG$nNodes){
      tmp <- catch$AG$Upstream[[i]]
      vecUps[i] <- sum(vec[tmp]*catch$SC$A[tmp])/sum(catch$SC$A[tmp])
    }
    #covariates[[paste(covariates_landcover[cov],"Ups",sep="")]] <- vecUps
  }
  
  
  ## function to locate points and derive geographical covariates
  clusters <- numeric(catch$AG$nNodes)
  
  #xy <- locator(1)
  xy <- vector("list"); xy$x <- 732918.8; xy$y <- 232304.2 # LUT
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 1
  
  xy <- vector("list"); xy$x <- 723246.2; xy$y <- 250333.9 # GON
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 2
  
  xy <- vector("list"); xy$x <- 734156.9; xy$y <- 252345.8 # WIS
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 3
  
  xy <- vector("list"); xy$x <- 731061.7; xy$y <- 242286.3 # NE3
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 4
  
  xy <- vector("list"); xy$x <- 734621.2; xy$y <- 252655.3 # GL2
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 5
  
  xy <- vector("list"); xy$x <- 729600.8; xy$y <- 257712.4 # GL1
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 6
  
  xy <- vector("list"); xy$x <- 728896.7; xy$y <- 257241.7 # GL1
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 6
  
  xy <- vector("list"); xy$x <- 723787.9; xy$y <- 245149.4 # DIE
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 7
  
  xy <- vector("list"); xy$x <- 738412.8; xy$y <- 227970.9 # TH8
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 8
  
  xy <- vector("list"); xy$x <- 732468.1; xy$y <- 232375.4 # TH7
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 9
  
  xy <- vector("list"); xy$x <- 725539.3; xy$y <- 238645.8 # TH6
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 10
  
  xy <- vector("list"); xy$x <- 724020; xy$y <- 244839.8 # TH5
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 11
  
  xy <- vector("list"); xy$x <- 727270; xy$y <- 247316 # NE2
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 12
  
  xy <- vector("list"); xy$x <- 723787.9; xy$y <- 250024.4 # TH4
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 13
  
  xy <- vector("list"); xy$x <- 724252.2; xy$y <- 250411.3 # NE1
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 14
  
  xy <- vector("list"); xy$x <- 723478.4; xy$y <- 257452.9 # TH3
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 15
  
  xy <- vector("list"); xy$x <- 729622.1; xy$y <- 257757.2 # TH2
  tmp <- locate_site(xy$x,xy$y,catch,showPlot = F)
  clusters[((1:catch$AG$nNodes) %in% catch$AG$Upstream[[tmp$AGnode]]) & clusters==0] <- 16
  
  clusters[clusters==0] <- 17
  
  #draw_thematic_catch(clusters,catch,discreteLevels=TRUE,colPalette = c("#000000",rainbow(17)))
  
  covariates_clusters <- c("LUT","GON","WIS","NE3","GL2","GL1","DIE","TH8",
  #                         1     2     3     4     5     6     7     8  
                           "TH7","TH6","TH5","NE2","TH4","NE1","TH3","TH2","TH1")
  #                         9     10    11    12    13    14    15    16    17
  
  for (i in 1:length(covariates_clusters)){
    vec <- numeric(catch$AG$nNodes)
    vec[clusters==i] <- 1
    covariates[[covariates_clusters[i]]] <- vec
  }
  
  # Z-normalize covariates
  Zcovariates <- data.frame(matrix(nrow=catch$AG$nNodes,ncol=0))
  for (i in 1:length(covariates)){
    Zcovariates[[names(covariates)[i]]] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i]) 
  }
  
  explist <- vector("list")
  explist[["covariates"]] <- covariates
  explist[["Zcovariates"]] <- Zcovariates
  
  invisible(explist)
}
