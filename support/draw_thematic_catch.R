
draw_thematic_catch <- function(OCN,theme=NA*numeric(OCN$AG$nNodes),
                                chooseAggregation=NULL,
                                discreteLevels=FALSE,
                                colLevels=NULL,
                                cutoff=FALSE,
                                colPalette=colorRampPalette(c("yellow","red","black")),
                                drawNodes=FALSE,
                                nodeType="upstream",
                                cex=2,
                                pch=21,
                                nanColor="#0099FF",
                                riverColor="#0099FF",
                                backgroundColor="#999999",
                                addLegend=TRUE,
                                zoom_xy_lim=NA){
  
  # initialization
  if (discreteLevels == FALSE) {
    if (is.null(colLevels)){
      
      minval <- min(theme[!(is.nan(theme))])
      maxval <- max(theme[!(is.nan(theme))])
      if (is.na(minval) & is.na(maxval)){
        minval <- 0; maxval <- 0;
      }
      N_colLevels <- 1000
      colLevels <- c(minval,maxval,N_colLevels)
    }
    minval <- colLevels[1]; maxval <- colLevels[2]; N_colLevels <- colLevels[3]
    
    if (minval==maxval) {maxval <- minval + 1}
    Breakpoints <- seq(minval,maxval,len = N_colLevels+1)
  } else if (discreteLevels == TRUE) {
    if (is.null(colLevels)){
      N_colLevels <- length(unique(theme[!is.nan(theme)]))
      Breakpoints <- c(sort(unique(theme[!is.nan(theme)])),2*max(theme[!is.nan(theme)]))
    } else {N_colLevels <- length(colLevels) - 1
    Breakpoints <- colLevels}}
  
  if (typeof(colPalette)=="closure") {
    colPalette <- colPalette(N_colLevels)
  } else if (typeof(colPalette)=="character") {
    colPalette <- colPalette[1:N_colLevels] }
  
  if (length(theme)!=OCN$AG$nNodes){
    stop('theme has invalid length')
  }
  
  
  if (length(cex)>1 && length(cex) != length(theme)){
    stop('cex has invalid length')
  }
  
  if (length(pch)>1 && length(pch) != length(theme)){
    stop('pch has invalid length')
  }
  
  
  X <- OCN$RN$X#[which( OCN$FD$toCM %in% chooseCM )]
  Y <- OCN$RN$Y#[which( OCN$FD$toCM %in% chooseCM )]
  Xc <- OCN$CM$XContour
  Yc <- OCN$CM$YContour
  
  
  if (length(cex)==1){
    cex_vec <- cex*rep(1,length(theme))
  } else {cex_vec <- cex}
  
  if (length(pch)==1){
    pch_vec <- pch*rep(1,length(theme))
  } else {pch_vec <- pch}
  
  AvailableNodes <- setdiff(1:OCN$RN$nNodes,OCN$RN$outlet)
  
  if (all(is.na(zoom_xy_lim))){
    zoom_xy_lim <- numeric(4)
    zoom_xy_lim[1] <- min(Xc); zoom_xy_lim[2] <- max(Xc) 
    zoom_xy_lim[3] <- min(Yc); zoom_xy_lim[4] <- max(Yc) 
    showAxes <- FALSE; labX <- ""; labY  <- ""
  } else {showAxes <- TRUE; labX <- "x"; labY <- "y"}
  
  plot(c(zoom_xy_lim[1],zoom_xy_lim[2],zoom_xy_lim[2],zoom_xy_lim[1],zoom_xy_lim[1]), 
       c(zoom_xy_lim[3],zoom_xy_lim[3],zoom_xy_lim[4],zoom_xy_lim[4],zoom_xy_lim[3]),
       type="n",xlab=labX,ylab=labY,asp=1,axes=showAxes)
  
  xy_lim <- par("usr")
  
  if (!is.null(backgroundColor)){ # don't show background when zooming
    polygon(Xc,Yc,col=backgroundColor,lty=0)    
  }
  
  maxA <- max(catch$RN$A[which(X >= xy_lim[1]& X <= xy_lim[2] & Y >= xy_lim[3] & Y <= xy_lim[4])])
  
  for (i in AvailableNodes){
    
    reach <- OCN$RN$toAGReach[i]
    if (OCN$RN$A[i]>=OCN$thrA & X[i] >= xy_lim[1]-1*catch$cellsize & X[i] <= xy_lim[2]+1*catch$cellsize &
        Y[i] >= xy_lim[3]-1*catch$cellsize & Y[i] <= xy_lim[4]+1*catch$cellsize ) {
      if (  ((is.nan(theme[reach])==TRUE | is.na(theme[reach])==TRUE)) ||
            ( cutoff==TRUE && (theme[reach] < min(Breakpoints) || theme[reach] > max(Breakpoints))) )  {
        hexcolor <- nanColor
      } else {
        val <- theme[reach]
        colvalue <- which(Breakpoints > val)[1] - 1  
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        
        hexcolor <- colPalette[colvalue]
      }
      if (drawNodes==TRUE){
        lines(c(X[i],X[OCN$RN$downNode[i]]),c(Y[i],Y[OCN$RN$downNode[i]]),
              lwd=0.1+2.9*(OCN$RN$A[i]/maxA)^0.5,col=riverColor)
      } else {
        lines(c(X[i],X[OCN$RN$downNode[i]]),c(Y[i],Y[OCN$RN$downNode[i]]),
              lwd=0.1+2.9*(OCN$RN$A[i]/maxA)^0.5,col=hexcolor)}
    }
  }
  
  if (drawNodes==TRUE){
    for (i in 1:OCN$AG$nNodes){
      if (is.nan(theme[i]) || (cutoff==TRUE && (theme[i] < min(Breakpoints) || theme[i] > max(Breakpoints)) )){
        hexcolor <- nanColor
      } else {
        colvalue <- which(Breakpoints > theme[i])[1] - 1
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        hexcolor <- colPalette[colvalue]
      }
      if (nodeType=="upstream") {
        points(OCN$AG$X[i],OCN$AG$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      } else if (nodeType=="downstream") {
        points(OCN$AG$XReach[i],OCN$AG$YReach[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      }
    }
  }
  
  if (addLegend) {
    if (discreteLevels==FALSE){
      tmp <- par()$plt[2]
      pr <- 1 - tmp
      image.plot(col=colPalette,legend.only=TRUE,zlim=c(minval,maxval),
                 smallplot=c(0.88, 0.9,par()$plt[3],par()$plt[4]))
    } else {
      if (is.null(colLevels)){
        str <- NULL
        for (level in 1:N_colLevels){
          str <- c(str, as.character(round(1000*Breakpoints[level])/1000) )}
      } else { 
        str <- vector(mode="character", N_colLevels)
        for (level in 1:(N_colLevels-1)){
          str[level] <- paste("[",as.character(round(1000*Breakpoints[level])/1000),"; ",
                              as.character(round(1000*Breakpoints[level+1])/1000),")",sep="")
        } 
        str[N_colLevels] <- paste("[",as.character(round(1000*Breakpoints[N_colLevels])/1000),"; ",
                                  as.character(round(1000*Breakpoints[N_colLevels+1])/1000),"]",sep="")
      }
      
      legend(x=1.01*max(Xc),y=max(Yc),
             str,fill=colPalette,ncol=ceiling(N_colLevels/20), xpd=TRUE, cex=0.8, bty="n")
    }
  }
  invisible(xy_lim)
  
}