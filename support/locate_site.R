locate_site <- function(X,Y,catch,
                        showPlot=TRUE,
                        title="",
                        zoom_xy_lim=NA){
  
  Xgrid <- catch$FD$X[which(abs(catch$FD$X-X)==min(abs(catch$FD$X-X)))][1]
  Ygrid <- catch$FD$Y[which(abs(catch$FD$Y-Y)==min(abs(catch$FD$Y-Y)))][1]
  
  shp.point(Xgrid,Ygrid,"tmp")
  system("mpiexec -n 8 moveoutletstostreams -p DTMp.tif -src DTMsrc.tif -o tmp.shp -om tmp2.shp")
  tmp2 <- read.shp("tmp2.shp")
  
  Xnew <- as.numeric(tmp2$shp[,2])
  Ynew <- as.numeric(tmp2$shp[,3])
  
  RNnode <- which(catch$RN$X==Xnew & catch$RN$Y==Ynew)
  AGnode <- catch$RN$toAGReach[RNnode]
  
  if (showPlot){
    Xmin <- min(X,Xnew); Xmax <- max(X,Xnew); Ymin <- min(Y,Ynew); Ymax <- max(Y,Ynew)
    
    if (all(is.na(zoom_xy_lim))){
      zoom_xy_lim <- c(Xmin-20*catch$cellsize,Xmax+20*catch$cellsize,Ymin-20*catch$cellsize,Ymax+20*catch$cellsize)
    }
    
    theme <- numeric(catch$AG$nNodes); theme[AGnode] <- 1
    xy_lim <- draw_thematic_catch(catch,theme, discreteLevels = T,colPalette = colorRampPalette(c("blue2","orange")),
                        zoom_xy_lim=zoom_xy_lim)
    title(title)
    legend(x=xy_lim[1]+catch$cellsize, y=xy_lim[4], legend=c("Original site","Relocated site","Assigned reach"),
           col=c("red","black","orange"), pch=c(15,20,NA),lty=c(0,0,1))
    points(X,Y,pch=15,col="red")
    points(Xnew,Ynew,pch=20,col="black")
    
  }
  
  explist <- vector("list")
  explist[["AGnode"]] <- AGnode
  explist[["RNnode"]] <- RNnode
  
  invisible(explist)
}
