binplot_area <- function(y, A, n_bins=10, col.shade="#909090AA",col.line="black",ylim=NULL, add=F, ylab="",
                    shade=T){
  if (is.null(ylim)){ylim=c(min(y, na.rm=T), max(y, na.rm=T))}
  A <- log10(A)
  bin_lims <- seq(min(A), max(A), length.out=n_bins+1)
  center_bins <- qlow <- qmedian <- qhigh <- numeric(n_bins)
  for (i in 1:n_bins){
    tmp <- which(A >= bin_lims[i] & A < bin_lims[i+1])
    qlow[i] <- quantile(y[tmp], 0.025, na.rm=T)
    qmedian[i] <- quantile(y[tmp], 0.5, na.rm=T)
    qhigh[i] <- quantile(y[tmp], 0.975, na.rm=T)
    center_bins[i] <- 0.5*(bin_lims[i] + bin_lims[i+1])
  }
  if (add==F){
    plot(10^center_bins,qmedian, ylim=ylim, xlim=c(3e5,5e8), type="n",log="x", 
         xlab="Drainage Area [m^2]", ylab=ylab,xaxt="n")
    axis(1, at=c(5e5, 5e6, 5e7, 5e8))
  } else {
    points(10^center_bins,qmedian,  type="n")
  }
  if (shade==T){
    polygon(c(10^center_bins, rev(10^center_bins)), c(qlow, rev(qhigh)), col=col.shade, border=NA)
  } else {
    lines(10^center_bins, qlow, lty=3, col=col.shade)
    lines(10^center_bins, qhigh, lty=3, col=col.shade)
  }
  lines(10^center_bins,qmedian,col=col.line)
}

fit_lm_error <- function(vec, A){
  
  vecA <- seq(min(log10(A)), max(log10(A)), length.out=100)
  vecA <- 10^vecA  
  bin_lims <- c(min(A), quantile(A,seq(0.1,0.9,0.1)), max(A)+1e-6)
  quantA_ID <- numeric(length(A))
  for (i in 1:10){
    quantA_ID[which(A >= bin_lims[i] & A < bin_lims[i+1])] <- i
  }
  
  coef_mat <- matrix(0,100,2)
  for (k in 1:100){
    set.seed(k)
    subs <- numeric(0)
    for (i in 1:10){subs <- c(subs,sample(which(quantA_ID==i), 50, replace=T))}
    lmod <- lm(vec[subs] ~ A[subs])
    ss <- summary(lmod)
    coef_mat[k,] <- ss$coefficients[,1]
  }
  mm <-  matrix(coef_mat[,2],100,1) %*% matrix(vecA,1,100) + matrix(rep(coef_mat[,1],100),100,100) 
  QQ <- apply(mm,2,quantile,probs=c(0.025,0.975))
  
  if (sum(coef_mat[,2]>0)>=95){
    flag <- "+"
  } else if (sum(coef_mat[,2]>0)>=50){
    flag <- "(+)"
  } else if (sum(coef_mat[,2]>0)>5){
    flag <- "(-)"
  } else {
    flag <- "-"
  }
  
  ll <- list(QQ=QQ, flag=flag)
}