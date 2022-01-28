rm(list=ls())

# This also assumes that MPICH2 is properly installed on your machine and that TauDEM command line executables exist
# MPICH2.  Obtain from http://www.mcs.anl.gov/research/projects/mpich2/
# Install following instructions at http://hydrology.usu.edu/taudem/taudem5.0/downloads.html.  
# It is important that you install this from THE ADMINISTRATOR ACCOUNT.

# TauDEM command line executables.  
# If on a PC download from http://hydrology.usu.edu/taudem/taudem5.0/downloads.html
# The install package will install to c:\program files\taudem or c:\program files (x86)\taudem set a 
# system path.  If you want to do this manually you can download the command line executables and place where you wish.
# If on a different system, download the source code and compile for your system.

library(raster)
library(shapefiles)
library(spam)
library(fields)
library(rgdal)
library(BayesianTools)

pathname <- dirname(parent.frame(2)$ofile)
setwd(pathname)

source("support/neigh.R")
source("support/draw_thematic_catch.R")
source("support/prepare_catchment.R")
source("support/eval_hydrology.R")
source("support/eval_covariates.R")
source("support/locate_site.R")
source("support/eDITH_model.R")

if (!file.exists("support/ThurData.rda")){
  thrA <- 2.5e5
  x_outlet <- 735162
  y_outlet <- 261667
  maxReachLength <- 1000
  DEM <- raster("support/DTM.asc")
  
  catch <- prepare_catchment(DEM, thrA, x_outlet, y_outlet, maxReachLength)
  
  sites <- read.csv("support/sites.csv")
  
  #par(ask=TRUE)
  samplingNodesAG <- samplingNodesRN <-  numeric(length(sites$X))
  for (i in 1:dim(sites)[1]){
    #pdf(paste(sites$Site_id[i],'.pdf',sep=""),width=12,height=9,paper="special")
    str <- sprintf("%s  -  X %.8g  - Y %.8g",sites$Site_id[i],sites$X[i],sites$Y[i])
    tmp <- locate_site(sites$X[i],sites$Y[i],catch,showPlot = F, title=str)
    samplingNodesAG[i] <- tmp$AGnode
    samplingNodesRN[i] <- tmp$RNnode
    #dev.off()
  }
  sites[["nodeAG"]] <- samplingNodesAG
  sites[["nodeRN"]] <- samplingNodesRN

  
  # hydrology ####
  catch <- eval_hydrology(catch,sites,"support/hydrology_data.csv")
  
  # choose covariates ####
  covariates_morphology <- c("DrainageArea","StreamOrder","MeanElevLoc","SlopeLoc","MeanSlopeUps")
  covariates_geol <- c("Alluvial","Alpine","Loess","Molasses","Moraines","Peat","Scree","Water")
  covariates_landcover <- c("Forest","Lake","Orchard","Rock","Swamp","UrbanArea")
  
  geology_shp <- "support/shapefile geology/PY_Surface_Base"
  landcover_shp <- "support/shapefile landcover/merge_landcover"
  
  tmp <- eval_covariates(catch, geology_shp, landcover_shp, covariates_morphology, covariates_geol, covariates_landcover)
  Zcovariates <- tmp$Zcovariates
  
  save(catch,sites,Zcovariates,file="support/ThurData.rda")
  
} else {
  load("support/ThurData.rda")
}


# edit eDNA data to adapt format ####
datafile <- "16S"
tmp <- read.csv(paste(pathname,"/calibration/",datafile,"/",datafile,".csv",sep=""))
names(tmp)[1] <- "Genus"
if (datafile=="16S"){
  genusNames <- tmp$Genus  
} else {
  genusNames <- substr(tmp$Genus,4,100)}
tmp <- as.data.frame(t(tmp[,-1]))
colnames(tmp) <- genusNames

time_id <- substr(row.names(tmp),2,2)
time_id[time_id=="1"] <- "spring"
time_id[time_id=="2"] <- "summer"
time_id[time_id=="3"] <- "autumn"

site_id <- substr(row.names(tmp),3,6)

tmp$time_id <- time_id
tmp$site_id <- site_id

value <- numeric(0)
time_id <- numeric(0)
site_id <- numeric(0)
genus <- character(0)
for (g in 1:length(genusNames)){
  genus <- c(genus, rep(as.character(genusNames[g]), length(tmp$time_id)))
  value <- c(value, tmp[[genusNames[g]]])
  time_id <- c(time_id, tmp$time_id)
  site_id <- c(site_id, tmp$site_id)
}
rm(tmp)
eDNA_data <- data.frame(matrix(nrow=length(value),ncol=0))
eDNA_data[["Genus"]] <- genus
eDNA_data[["site_id"]] <- site_id
eDNA_data[["time_id"]] <- time_id
eDNA_data[["value"]] <- value


# run eDITH for a genus and a season
for (genus in unique(eDNA_data$Genus)){
  cat(sprintf("%s \n",genus))
  for (tt in 1:length(unique(eDNA_data$time_id))){
    
    time_id <- unique(eDNA_data$time_id)[tt] 
    cat(sprintf("    %s \n",time_id))
    fname <- paste(genus,"_",time_id,".rda",sep="")
    tmp <- which(eDNA_data$Genus == genus & eDNA_data$time_id == time_id)
    eDNA_values <- eDNA_data$value[tmp]
    
    eDNA_siteAG <- numeric(length(tmp))
    for (i in 1:length(tmp)){
      eDNA_siteAG[i] <- sites$nodeAG[which(sites$Site_id==eDNA_data$site_id[tmp[i]])]
    }
    data <- vector("list")
    data[["eDNA_values"]] <- eDNA_values
    data[["eDNA_siteAG"]] <- eDNA_siteAG
    out <- numeric(length(Zcovariates) + 2)
    
    if(!file.exists(paste(pathname,"/calibration/",datafile,"/",fname,sep=""))){
      save(out,data,file=paste(pathname,"/calibration/",datafile,"/",fname,sep=""))
      if (sum(eDNA_values)==0){ # 
        cat("     no eDNA found \n")
        gD <- NULL
        qq <- NULL
      } else {
        
        # real implementation
        density = function(par){
          d1 = dlnorm(par[1], log(5), sqrt(log(5)-log(4)), log =TRUE) # median of 5 h, mode of 4 h
          # for bacteria: try with median of 1 h, mode of 0.5 h
          d2 = dunif(par[2], 0, 1, log =TRUE)
          d_beta <- numeric(length(Zcovariates))
          for (i in 1:length(Zcovariates)){
            d_beta[i] <- dnorm(par[i+2], 0, 3, log=TRUE)
          }
          return(d1 + d2 + sum(d_beta))
        }
        
        sampler = function(n=1){
          d1 = rlnorm(n, log(5), sqrt(log(5)-log(4)))
          d2 = runif(n, 0, 1)
          d_beta <- numeric(length(Zcovariates))
          for (i in 1:length(Zcovariates)){
            d_beta[i] <- rnorm(n, 0, 3)
          }
          return(c(d1,d2,d_beta))
        }
        
        prior <- createPrior(density = density, sampler = sampler,
                             lower = c(0,0,rep(-10,length(Zcovariates))),
                             upper = c(24,1,rep(10,length(Zcovariates))) )
        
        likelihood <- function(param){
          tau <- param[1]*3600
          p0 <- param[2]
          beta <- param[-c(1,2)]
          p <- p0 *exp(as.numeric(as.matrix(Zcovariates) %*% as.matrix(beta)))
          Conc <- eDITH_model(p=p, tau=tau, width=catch$AG$width[[time_id]], leng=catch$AG$leng, velocity=catch$AG$velocity[[time_id]],
                              Q=catch$AG$Q[[time_id]], nNodes=catch$AG$nNodes, perm=catch$AG$perm, downNode=catch$AG$downNode)
          Conc <- Conc[data$eDNA_siteAG]
          loglik <- sum( log(1/(Conc + 1)) + data$eDNA_values*log(Conc/(Conc + 1)) )
          return(loglik)
        }
        
        set.seed(1)
        options(warn=-1) # ignore warnings
        setUp <- createBayesianSetup(likelihood, prior=prior)
        
        settings <- list(iterations = 3.5e6,  message = FALSE, burnin=5e5)
        t0 <- Sys.time()
        suppressMessages( out <- runMCMC(bayesianSetup = setUp, sampler = "DREAMzs", settings = settings) )
        
        cat(sprintf("     100%% done  -  Date: %16s  -  Elapsed time: %.2f s \n",
                    Sys.time(),difftime(Sys.time(),t0,units="secs")))
        
        qq <- getCredibleIntervals(rbind(out$chain[[1]],out$chain[[2]],out$chain[[3]]))
        gD <- gelmanDiagnostics(out)
        out <- MAP(out)
    
      } 
      save(out,gD,qq,data,file=paste(pathname,"/calibration/",datafile,"/",fname,sep=""))
    } else {
      cat("     Already processed \n")
    }
  }
  
  cat("\n")
}




