rm(list=ls())

pathname <- dirname(parent.frame(2)$ofile)
setwd(pathname)

library(raster)
library(shapefiles)
library(spam)
library(fields)
library(rgdal)
library(BayesianTools)
library(data.table)
library(betapart)
library(vegan)
library(car) # allows Anova with type III sum of squares

load("support/ThurData.rda")

source("support/eDITH_model.R")
source("support/draw_thematic_catch.R")
source("support/extract_unconnected_pairs.R")
source("support/fit_lm_error.R")

# read all eDNA datasets ####
datafile <- c("12S","COI","16S")
for (fn in datafile){
  tmp <- read.csv(paste(pathname,"/calibration/",fn,"/",fn,".csv",sep=""))
  
  names(tmp)[1] <- "Genus"
  if (fn=="16S"){
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
    genus <- c(genus, rep(genusNames[g], length(tmp$time_id)))
    value <- c(value, tmp[[genusNames[g]]])
    time_id <- c(time_id, tmp$time_id)
    site_id <- c(site_id, tmp$site_id)
  }
  rm(tmp)
  
  eval(parse(text=paste("data_",fn,"=data.frame(matrix(nrow=length(value),ncol=0))",sep="")))
  eval(parse(text=paste("data_",fn,"[['Genus']] <- genus",sep="")))
  eval(parse(text=paste("data_",fn,"[['time_id']] <- time_id",sep="")))
  eval(parse(text=paste("data_",fn,"[['site_id']] <- site_id",sep="")))
  eval(parse(text=paste("data_",fn,"[['value']] <- value",sep="")))
}

# read functions
funct.df <- read.csv(paste(pathname,"/calibration/Function_genus_20210120.csv",sep=""))
funct.df$description[funct.df$description=='Macro'] <- 'Invert'
funct.df$description[funct.df$description=='Micro'] <- 'Invert'

# remove Barbus, Gobio, Phoxinus from COI dataset
genus_vec <- intersect(data_12S$Genus, data_COI$Genus)
data_COI <- data_COI[is.na(match(data_COI$Genus, genus_vec)), ]  

# read eDITH model results and calculate PD, PA ####
threshold <- 0.5
p_matrix <- PD_matrix <- PA_matrix <- tau_matrix <- covariate_significance <-  vector("list",3)

time_id_vec <- c("spring","summer","autumn")
names(PD_matrix) <- names(PA_matrix) <- names(tau_matrix) <- names(covariate_significance) <- time_id_vec
allGenusNames <- c(unique(data_12S$Genus), unique(data_COI$Genus), unique(data_16S$Genus))

bacteria_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Bacteria"]
invert_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Invert"]
fish_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Fish"]

for (time_id in time_id_vec){
  p_matrix[[time_id]] <- vector("list",length(allGenusNames))
  names(p_matrix[[time_id]]) <- allGenusNames
  PD_matrix[[time_id]] <- vector("list",length(allGenusNames))
  names(PD_matrix[[time_id]]) <- allGenusNames
  PA_matrix[[time_id]] <- data.frame(matrix(0,catch$AG$nNodes,length(allGenusNames)))
  names(PA_matrix[[time_id]]) <- allGenusNames
  tau_matrix[[time_id]] <- data.frame(matrix(0,1,length(allGenusNames)))
  names(tau_matrix[[time_id]]) <- allGenusNames
  covariate_significance[[time_id]] <- data.frame(matrix(0,length(Zcovariates),length(allGenusNames)),row.names=names(Zcovariates))
  names(covariate_significance[[time_id]]) <- allGenusNames
  for (datafile in c("12S","COI","16S")){
    eval(parse(text=paste("eDNA_data=data_",datafile,sep="")))
    genus_vec <- unique(eDNA_data$Genus)
    for (genus in genus_vec){
      fname <- paste(pathname,"/calibration/",datafile,"/",genus,"_",time_id,".rda",sep="")
      if (file.exists(fname)){
        load(fname)
        cat(sprintf("%s  -   %s  -  %s \n",genus,datafile,time_id))
        if (is.list(out)){
          map <- out$parametersMAP
        } else {
          map <- out
        }
        tau <- as.numeric(map[1])*3600
        p0 <-  as.numeric(map[2])
        beta <- as.numeric(map[-c(1,2)])
        if (is.null(qq)){
          cI <- matrix(0,2,length(beta))
        } else {cI <- qq[,3:(length(beta)+2)]}
        
        p <- p0 *exp(as.numeric(as.matrix(Zcovariates) %*% as.matrix(beta)))
        
        conc <- eDITH_model(p=p, tau=tau, width=catch$AG$width[[time_id]], leng=catch$AG$leng, velocity=catch$AG$velocity[[time_id]],
                            Q=catch$AG$Q[[time_id]], nNodes=catch$AG$nNodes, perm=catch$AG$perm, downNode=catch$AG$downNode)
        
        Conc <- conc[data$eDNA_siteAG]
        loglik <- sum( log(1/(Conc + 1)) + data$eDNA_values*log(Conc/(Conc + 1)) )
        
        q <- numeric(catch$AG$nNodes)
        for (i in 1:catch$AG$nNodes){
          tmp <- which(catch$AG$downNode==i)
          sumQ <- 0
          if (length(tmp)>0){
            for (j in 1:length(tmp)){
              sumQ <- sumQ + catch$AG$Q[[time_id]][tmp[j]]
            }
          }
          q[i] <- catch$AG$Q[[time_id]][i] - sumQ
        }
        
        NU <- p*catch$AG$width[[time_id]]*catch$AG$leng/q*exp(-catch$AG$leng/catch$AG$velocity[[time_id]]/tau)
        PD <- NU/(1+NU)
        
        #cat("\n")
      } else {
        PD <- numeric(catch$AG$nNodes)
      }
      PA <- as.numeric(PD >= threshold)
      if (tau==0){tau <- NA}
      
      #PA[[time_id]] <- PA[[time_id]] + PD
      p_matrix[[time_id]][[genus]] <- p
      PD_matrix[[time_id]][[genus]] <- PD #as.numeric(PD > threshold)
      PA_matrix[[time_id]][[genus]] <- PA
      covariate_significance[[time_id]][[genus]] <- -as.numeric(cI[2,]<0)+as.numeric(cI[1,]>0)
      tau_matrix[[time_id]][[genus]] <- tau/3600
    }
  }
}

# alpha diversity ####
cat('\n'); cat('alpha diversity... \n')
alpha.diversity <- vector("list",3)
names(alpha.diversity) <- time_id_vec
functions_vec <- unique(funct.df$description)[1:3]
for (time_id in time_id_vec){
  alpha.diversity[[time_id]] <- data.frame(matrix(0,catch$AG$nNodes,length(functions_vec) + 1))
  names(alpha.diversity[[time_id]]) <- c(functions_vec,"All")
  for (func in functions_vec){
    genera <- allGenusNames[which(funct.df$description[match(allGenusNames, funct.df$Taxa)]==func)]
    alpha <- numeric(catch$AG$nNodes)
    for (genus in genera){
      alpha <- alpha + PA_matrix[[time_id]][[genus]]
    }
    alpha.diversity[[time_id]][[func]] <- alpha
    alpha.diversity[[time_id]]$All <- alpha.diversity[[time_id]]$All + alpha
  }
}

# evaluate beta diversity ####
cat('\n'); cat('beta diversity... \n')
beta.div.bact <- vector("list",3)
names(beta.div.bact) <- time_id_vec
for (time_id in time_id_vec){
  tmp <- PA_matrix[[time_id]][,match(funct.df$Taxa[which(funct.df$description=="Bacteria")],allGenusNames)]
  bb <- betapart.core(tmp)
  beta.div.bact[[time_id]] <- beta.pair(bb,index.family = "jac")
}

beta.div.invert <- vector("list",3)
names(beta.div.invert) <- time_id_vec
for (time_id in time_id_vec){
  tmp <- PA_matrix[[time_id]][,match(funct.df$Taxa[which(funct.df$description=="Invert")],allGenusNames)]
  bb <- betapart.core(tmp)
  beta.div.invert[[time_id]] <- beta.pair(bb,index.family = "jac")
  # NaN occure if alpha div of one of the two sites of the pair is 0;
  # if only one of the two sites has alpha=0, then beta.jac=1, but beta.jne and beta.jtu are undetermined
  # if both sites have alpha=0, then all beta are undetermined
  # for the moment, let's consider beta.jac=1, beta.jtu=1, beta.jne=0 in these cases
  beta.div.invert[[time_id]]$beta.jtu[is.nan(beta.div.invert[[time_id]]$beta.jtu)] <- 1
  beta.div.invert[[time_id]]$beta.jac[is.nan(beta.div.invert[[time_id]]$beta.jac)] <- 1
  beta.div.invert[[time_id]]$beta.jne[is.nan(beta.div.invert[[time_id]]$beta.jne)] <- 0
}

beta.div.fish <- vector("list",3)
names(beta.div.fish) <- time_id_vec
for (time_id in time_id_vec){
  tmp <- PA_matrix[[time_id]][,match(funct.df$Taxa[which(funct.df$description=="Fish")],allGenusNames)]
  bb <- betapart.core(tmp)
  beta.div.fish[[time_id]] <- beta.pair(bb,index.family = "jac")
  # NaN occure if alpha div of one of the two sites of the pair is 0;
  # if only one of the two sites has alpha=0, then beta.jac=1, but beta.jne and beta.jtu are undetermined
  # if both sites have alpha=0, then all beta are undetermined
  # for the moment, let's consider beta.jac=1, beta.jtu=1, beta.jne=0 in these cases
  beta.div.fish[[time_id]]$beta.jtu[is.nan(beta.div.fish[[time_id]]$beta.jtu)] <- 1
  beta.div.fish[[time_id]]$beta.jac[is.nan(beta.div.fish[[time_id]]$beta.jac)] <- 1
  beta.div.fish[[time_id]]$beta.jne[is.nan(beta.div.fish[[time_id]]$beta.jne)] <- 0
}

# pairwise beta diversity in headwaters vs large reaches ####
cat('\n'); cat('beta diversity upstream/downstream... \n')

#  pick only pairs of flow-unconnected sites such that each node appears only in one pair
if (!file.exists("support/beta_HeadMain.rda")){
  beta_HeadMain <- vector("list",4); names(beta_HeadMain) <- c("bact","invert","fish")
  # Bacteria
  beta_HeadMain$bact <- vector("list",3); names(beta_HeadMain$bact) <- time_id_vec
  for (time_id in time_id_vec){
    beta_HeadMain$bact[[time_id]]$p_values <- numeric(100)
    beta_HeadMain$bact[[time_id]]$mean_eff <- data.frame(matrix(0,100,2))
    names(beta_HeadMain$bact[[time_id]]$mean_eff) <- c("Up","Down")
  }
  set.seed(3)
  for (i in 1:100){
    cat(sprintf("%d \n",i))
    ls <- extract_unconnected_pairs(catch)
    par(mfrow=c(1,3))
    for (time_id in time_id_vec){
      betaMat <- as.matrix(beta.div.bact[[time_id]]$beta.jac)
      DF <- data.frame(beta=c(betaMat[ls$headIndices], betaMat[ls$downIndices]),
                       streamDist=c(catch$AG$downstreamLengthUnconnected[ls$headIndices], catch$AG$downstreamLengthUnconnected[ls$downIndices]),
                       location=c(rep("Up",length(ls$headIndices)), rep("Down",length(ls$downIndices))) )
      aa <- aov(beta ~ streamDist + location, data=DF)
      AA <- Anova(aa, type="III")
      boxplot(beta ~ location, data=DF, notch=T); title(sprintf("%s  -  p: %.4f",time_id,AA$`Pr(>F)`[3]))
      beta_HeadMain$bact[[time_id]]$p_values[i] <- AA$`Pr(>F)`[3]
      beta_HeadMain$bact[[time_id]]$mean_eff$Up[i] <- mean(DF$beta[DF$location=="Up"]) # this doesn't require ancova!
      beta_HeadMain$bact[[time_id]]$mean_eff$Down[i] <- mean(DF$beta[DF$location=="Down"])
    }
  }
  
  # Invert
  beta_HeadMain$invert <- vector("list",3); names(beta_HeadMain$invert) <- time_id_vec
  for (time_id in time_id_vec){
    beta_HeadMain$invert[[time_id]]$p_values <- numeric(100)
    beta_HeadMain$invert[[time_id]]$mean_eff <- data.frame(matrix(0,100,2))
    names(beta_HeadMain$invert[[time_id]]$mean_eff) <- c("Up","Down")
  }
  set.seed(4)
  for (i in 1:100){
    cat(sprintf("%d \n",i))
    par(mfrow=c(1,3))
    for (time_id in time_id_vec){
      ls <- extract_unconnected_pairs(catch, alpha.diversity[[time_id]]$Invert) # here pair extraction must be done for each season (some pairs have alpha=0)
      betaMat <- as.matrix(beta.div.invert[[time_id]]$beta.jac)
      DF <- data.frame(beta=c(betaMat[ls$headIndices], betaMat[ls$downIndices]),
                       streamDist=c(catch$AG$downstreamLengthUnconnected[ls$headIndices], catch$AG$downstreamLengthUnconnected[ls$downIndices]),
                       location=c(rep("Up",length(ls$headIndices)), rep("Down",length(ls$downIndices))) )
      aa <- aov(beta ~ streamDist + location, data=DF)
      AA <- Anova(aa, type="III")
      boxplot(beta ~ location, data=DF, notch=T); title(sprintf("%s  -  p: %.4f",time_id,AA$`Pr(>F)`[3]))
      beta_HeadMain$invert[[time_id]]$p_values[i] <- AA$`Pr(>F)`[3]
      beta_HeadMain$invert[[time_id]]$mean_eff$Up[i] <- mean(DF$beta[DF$location=="Up"])
      beta_HeadMain$invert[[time_id]]$mean_eff$Down[i] <- mean(DF$beta[DF$location=="Down"])
    }
  }
  
  # Fish
  beta_HeadMain$fish <- vector("list",3); names(beta_HeadMain$fish) <- time_id_vec
  for (time_id in time_id_vec){
    beta_HeadMain$fish[[time_id]]$p_values <- numeric(100)
    beta_HeadMain$fish[[time_id]]$mean_eff <- data.frame(matrix(0,100,2))
    names(beta_HeadMain$fish[[time_id]]$mean_eff) <- c("Up","Down")
  }
  set.seed(5)
  for (i in 1:100){
    cat(sprintf("%d \n",i))
    par(mfrow=c(1,3))
    for (time_id in time_id_vec){
      ls <- extract_unconnected_pairs(catch, alpha.diversity[[time_id]]$Fish) # here pair extraction must be done for each season (some pairs have alpha=0)
      betaMat <- as.matrix(beta.div.fish[[time_id]]$beta.jac)
      DF <- data.frame(beta=c(betaMat[ls$headIndices], betaMat[ls$downIndices]),
                       streamDist=c(catch$AG$downstreamLengthUnconnected[ls$headIndices], catch$AG$downstreamLengthUnconnected[ls$downIndices]),
                       location=c(rep("Up",length(ls$headIndices)), rep("Down",length(ls$downIndices))) )
      aa <- aov(beta ~ streamDist + location, data=DF)
      AA <- Anova(aa, type="III")
      boxplot(beta ~ location, data=DF, notch=T); title(sprintf("%s  -  p: %.4f",time_id,AA$`Pr(>F)`[3]))
      beta_HeadMain$fish[[time_id]]$p_values[i] <- AA$`Pr(>F)`[3]
      beta_HeadMain$fish[[time_id]]$mean_eff$Up[i] <- mean(DF$beta[DF$location=="Up"])
      beta_HeadMain$fish[[time_id]]$mean_eff$Down[i] <- mean(DF$beta[DF$location=="Down"])
    }
  }
  save(beta_HeadMain,file="support/beta_HeadMain.rda")
} else {load(file="support/beta_HeadMain.rda")}



# temporal variation in beta diversity ####
cat('\n'); cat('temporal beta diversity... \n')
beta_temp <- vector("list",4)
names(beta_temp) <- c("bact","invert","fish")
# Bacteria
beta_temp$bact <- vector("list",3)
names(beta_temp$bact) <- c("jtu","jne","jac")
beta_temp$bact$jac <- beta_temp$bact$jne <- beta_temp$bact$jtu <- data.frame(matrix(0,catch$AG$nNodes,3))
names(beta_temp$bact$jac) <- names(beta_temp$bact$jtu) <- names(beta_temp$bact$jne) <- c("spring_summer","spring_autumn","summer_autumn")
for (i in 1:catch$AG$nNodes){
  tmp_df <- data.frame(matrix(0,3,sum(funct.df$description=="Bacteria",na.rm=T) ))
  k <- 1 # season counter
  for (time_id in time_id_vec){
    tmp_df[k, ] <- PA_matrix[[time_id]][i,match(funct.df$Taxa[which(funct.df$description=="Bacteria")],allGenusNames)]
    k <- k + 1
  }
  bb <- betapart.core(tmp_df)
  bp <- beta.pair(bb,index.family = "jac")
  beta_temp$bact$jtu[i,] <- as.vector(bp$beta.jtu)
  beta_temp$bact$jne[i,] <- as.vector(bp$beta.jne)
  beta_temp$bact$jac[i,] <- as.vector(bp$beta.jac)
}

# Invert
beta_temp$invert <- vector("list",3)
names(beta_temp$invert) <- c("jtu","jne","jac")
beta_temp$invert$jac <- beta_temp$invert$jne <- beta_temp$invert$jtu <- data.frame(matrix(0,catch$AG$nNodes,3))
names(beta_temp$invert$jac) <- names(beta_temp$invert$jtu) <- names(beta_temp$invert$jne) <- c("spring_summer","spring_autumn","summer_autumn")
for (i in 1:catch$AG$nNodes){
  tmp_df <- data.frame(matrix(0,3,sum(funct.df$description=="Invert",na.rm=T) ))
  k <- 1 # season counter
  for (time_id in time_id_vec){
    tmp_df[k, ] <- PA_matrix[[time_id]][i,match(funct.df$Taxa[which(funct.df$description=="Invert")],allGenusNames)]
    k <- k + 1
  }
  bb <- betapart.core(tmp_df)
  bp <- beta.pair(bb,index.family = "jac")
  beta_temp$invert$jtu[i,] <- as.vector(bp$beta.jtu)
  beta_temp$invert$jne[i,] <- as.vector(bp$beta.jne)
  beta_temp$invert$jac[i,] <- as.vector(bp$beta.jac)
}

# Fish
beta_temp$fish <- vector("list",3)
names(beta_temp$fish) <- c("jtu","jne","jac")
beta_temp$fish$jac <- beta_temp$fish$jne <- beta_temp$fish$jtu <- data.frame(matrix(0,catch$AG$nNodes,3))
names(beta_temp$fish$jac) <- names(beta_temp$fish$jtu) <- names(beta_temp$fish$jne) <- c("spring_summer","spring_autumn","summer_autumn")
for (i in 1:catch$AG$nNodes){
  tmp_df <- data.frame(matrix(0,3,sum(funct.df$description=="Fish",na.rm=T) ))
  k <- 1 # season counter
  for (time_id in time_id_vec){
    tmp_df[k, ] <- PA_matrix[[time_id]][i,match(funct.df$Taxa[which(funct.df$description=="Fish")],allGenusNames)]
    k <- k + 1
  }
  bb <- betapart.core(tmp_df)
  bp <- beta.pair(bb,index.family = "jac")
  beta_temp$fish$jtu[i,] <- as.vector(bp$beta.jtu)
  beta_temp$fish$jne[i,] <- as.vector(bp$beta.jne)
  beta_temp$fish$jac[i,] <- as.vector(bp$beta.jac)
}

# evaluate alpha and temporal beta for war eDNA data ####
# alpha div for eDNA data
A_data <- catch$AG$A[sites$nodeAG[1:73]]
alpha.div.data <- vector("list",3)
names(alpha.div.data) <- time_id_vec
for (time_id in time_id_vec){
  alpha.div.data[[time_id]] <- vector("list",4)
  names(alpha.div.data[[time_id]]) <- c("Bacteria","Invert","Fish","All")
}

# Bacteria
for (time_id in time_id_vec){
  alpha.div.data[[time_id]]$Bacteria <- numeric(length(A_data))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    alpha.div.data[[time_id]]$Bacteria[i] <- sum(data_16S$value[data_16S$time_id==time_id & data_16S$site_id==site]>0)
  }
}

# Invert
vv <- match(data_COI$Genus,invert_names)>0
vv[is.na(vv)] <- FALSE
data_invert <- data_COI[vv,]
for (time_id in time_id_vec){
  alpha.div.data[[time_id]]$Invert <- numeric(length(A_data))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    alpha.div.data[[time_id]]$Invert[i] <- sum(data_invert$value[data_invert$time_id==time_id & data_invert$site_id==site]>0)
  }
}

# Fish
vv <- match(data_COI$Genus, setdiff(intersect(unique(data_COI$Genus),fish_names), intersect(unique(data_COI$Genus),unique(data_12S$Genus))))
data_fish2 <- data_COI[which(!is.na(vv)), ]

data_fish <- data.frame(Genus=c(data_12S$Genus, data_fish2$Genus), time_id=c(data_12S$time_id,data_fish2$time_id), 
                        site_id=c(data_12S$site_id,data_fish2$site_id), value=c(data_12S$value,data_fish2$value) )
for (time_id in time_id_vec){
  alpha.div.data[[time_id]]$Fish <- numeric(length(A_data))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    alpha.div.data[[time_id]]$Fish[i] <- sum(data_fish$value[data_fish$time_id==time_id & data_fish$site_id==site]>0)
  }
}

# temporal beta for eDNA data
beta_temp.data <- vector("list",3)
names(beta_temp.data) <- c("bact","invert","fish")
# fish
beta_temp.data$fish <- vector("list",3)
names(beta_temp.data$fish) <- c("spring_summer","spring_autumn","summer_autumn")

PA.data <- vector("list",3)
names(PA.data) <- c("fish","invert","bact")
PA.data$fish <- vector("list",3)
names(PA.data$fish) <- c("spring","summer","autumn")
for (time_id in time_id_vec){
  PA.data$fish[[time_id]] <- matrix(0,length(A_data),length(fish_names))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    for (j in 1:length(fish_names)){
      genus <- fish_names[j]
      tmp <- data_fish$value[data_fish$time_id==time_id & data_fish$site==site & data_fish$Genus==genus]
      if (length(tmp)==0){tmp <- 0}
      if (tmp>0){
        PA.data$fish[[time_id]][i,j] <- 1
      }
    }
  }
}
bb <- beta.temp(PA.data$fish$spring,PA.data$fish$summer, index.family="jaccard")
beta_temp.data$fish$spring_summer <- bb$beta.jac
bb <- beta.temp(PA.data$fish$spring,PA.data$fish$autumn, index.family="jaccard")
beta_temp.data$fish$spring_autumn <- bb$beta.jac
bb <- beta.temp(PA.data$fish$summer,PA.data$fish$autumn, index.family="jaccard")
beta_temp.data$fish$summer_autumn <- bb$beta.jac

# invert
PA.data$invert <- vector("list",3)
names(PA.data$invert) <- c("spring","summer","autumn")
for (time_id in time_id_vec){
  PA.data$invert[[time_id]] <- matrix(0,length(A_data),length(invert_names))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    for (j in 1:length(invert_names)){
      genus <- invert_names[j]
      tmp <- data_invert$value[data_invert$time_id==time_id & data_invert$site==site & data_invert$Genus==genus]
      if (length(tmp)==0){tmp <- 0}
      if (tmp>0){
        PA.data$invert[[time_id]][i,j] <- 1
      }
    }
  }
}
bb <- beta.temp(PA.data$invert$spring,PA.data$invert$summer, index.family="jaccard")
beta_temp.data$invert$spring_summer <- bb$beta.jac
bb <- beta.temp(PA.data$invert$spring,PA.data$invert$autumn, index.family="jaccard")
beta_temp.data$invert$spring_autumn <- bb$beta.jac
bb <- beta.temp(PA.data$invert$summer,PA.data$invert$autumn, index.family="jaccard")
beta_temp.data$invert$summer_autumn <- bb$beta.jac

# bact
PA.data$bact <- vector("list",3)
names(PA.data$bact) <- c("spring","summer","autumn")
for (time_id in time_id_vec){
  PA.data$bact[[time_id]] <- matrix(0,length(A_data),length(bacteria_names))
  for (i in 1:length(A_data)){
    site <- sites$Site_id[i]
    for (j in 1:length(bacteria_names)){
      genus <- bacteria_names[j]
      tmp <- data_16S$value[data_16S$time_id==time_id & data_16S$site==site & data_16S$Genus==genus]
      if (length(tmp)==0){tmp <- 0}
      if (tmp>0){
        PA.data$bact[[time_id]][i,j] <- 1
      }
    }
  }
}
bb <- beta.temp(PA.data$bact$spring,PA.data$bact$summer, index.family="jaccard")
beta_temp.data$bact$spring_summer <- bb$beta.jac
bb <- beta.temp(PA.data$bact$spring,PA.data$bact$autumn, index.family="jaccard")
beta_temp.data$bact$spring_autumn <- bb$beta.jac
bb <- beta.temp(PA.data$bact$summer,PA.data$bact$autumn, index.family="jaccard")
beta_temp.data$bact$summer_autumn <- bb$beta.jac


# tornado plots for covariates ####
cat('\n'); cat('covariate plots... \n')
tornado.plot <- function(data_raw, Zcovariates, xlim=NULL, n.comp=15){
  
  data <- matrix(0,2,length(Zcovariates))
  colnames(data) <- names(Zcovariates)
  data[1,] <- rowSums(as.matrix(data_raw)==1)
  data[2,] <- -rowSums(as.matrix(data_raw)==-1)
  #old.par <- par()
  #suppressWarnings(par(c(old.par,mar=c(5.1,7,4.1,2.1))))
  ll <- sort(colSums(abs(data)),index.return=T)
  if(is.null(xlim)){
    x <- pretty(c(min(data),max(data)))
    xlim <- c(min(x), max(x))
  } 
  x <- pretty(xlim)
  barplot(data[1,ll$ix[(length(ll$ix)-n.comp+1):length(ll$ix)]], horiz = T, las=1, xlim = xlim, xaxt='n', ylab = '',
          beside=T, col=c('springgreen'))
  barplot(data[2,ll$ix[(length(ll$ix)-n.comp+1):length(ll$ix)]], horiz = T, las=1, xlim = xlim, xaxt='n', ylab = '',
          beside=T, col=c('indianred2'), add = TRUE)
  axis(1, at=x, las=TRUE)
  #suppressWarnings(par(old.par))
}
seas.cols <- c('#B010D0','#D0B010','#D35000')

# Fig. 1 ####
Xgrid <- seq(min(catch$FD$X)-catch$cellsize, max(catch$FD$X)+catch$cellsize, catch$cellsize)
Ygrid <- seq(min(catch$FD$Y)-catch$cellsize, max(catch$FD$Y)+catch$cellsize, catch$cellsize)
Zmat <- matrix(NA,nrow=length(Ygrid),ncol=length(Xgrid))
for (i in 1:length(catch$FD$X)){
  ind_row <- which(Ygrid==catch$FD$Y[i])
  ind_col <- which(Xgrid==catch$FD$X[i])
  if (!is.na(catch$FD$toSC[i])){Zmat[ind_row,ind_col] <- catch$FD$Z[i]}
}

png(file="Fig1_relief.png",width=8/2.54, height=10/2.54,units="in",res=600)
par(mar=c(0,0,0,0))
image(Xgrid,Ygrid,t(Zmat),col=colorRampPalette(terrain.colors(20),bias=2)(2100),
      breaks=400:2500,xlab = " ", ylab = " ", axes = FALSE, asp = 1)
dev.off()

pdf(file="Fig1_contour.pdf",width=8/2.54, height=10/2.54)
par(mar=c(0,0,0,0))
plot(c(min(catch$FD$X),max(catch$FD$X)), c(min(catch$FD$Y),max(catch$FD$Y)),
     type="n",xlab = " ", ylab = " ", axes = FALSE, asp = 1)
for (i in 1:length(catch$RN$X)){
  lines(c(catch$RN$X[i],catch$RN$X[catch$RN$downNode[i]]),c(catch$RN$Y[i],catch$RN$Y[catch$RN$downNode[i]]),
        lwd=0.1+2.9*(catch$RN$A[i]/max(catch$RN$A))^0.5,col="#0070DD")
}
polygon(catch$CM$X,catch$CM$Y)
points(sites$X[1:73], sites$Y[1:73], pch=21, bg="#DD0070", cex=0.5)
points(sites$X[74:77], sites$Y[74:77], pch=22, bg="#808080", cex=0.5)
lines(c(711000,716000),c(221000,221000))
lines(c(711000,716000),c(221600,221600))
lines(c(711000,711000),c(221000,221600))
lines(c(712000,712000),c(221000,221600))
lines(c(713000,713000),c(221000,221600))
lines(c(714000,714000),c(221000,221600))
lines(c(715000,715000),c(221000,221600))
lines(c(716000,716000),c(221000,221600))
dev.off()

pdf(file="Fig1_colormap.pdf",width=8/2.54, height=10/2.54)
plot(c(min(catch$FD$X),max(catch$FD$X)), c(min(catch$FD$Y),max(catch$FD$Y)),
     type="n",xlab = " ", ylab = " ", axes = FALSE, asp = 1)
image.plot(Xgrid,Ygrid,t(Zmat),col=colorRampPalette(terrain.colors(20),bias=2)(2100),
           breaks=400:2500,legend.only=TRUE)
dev.off()

# Fig. 2 ####
group.names <- c("bact","invert","fish")
pdf(file="Fig2_insets.pdf",width=20/2.54,height=15/2.54)
par(mfrow=c(3,3))
p1 <- hist(alpha.diversity$spring$Fish,seq(-0.5,6.5,1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$spring$Fish,seq(-0.5,6.5,1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Blues 3")[2],  ylim=c(0,0.5), ylab="Frequency", xlab="",main="",axes=F); title('Fish')
plot(p2,col=rgb(0,0,0,0.25),  add=T, main="")
axis(1, pos=0); axis(2,pos=-0.5,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$spring$Invert,seq(0,30,3), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$spring$Invert,seq(0,30,3), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Greens 3")[2], ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); title('Invertebrates')
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$spring$Bacteria,seq(0,120,12), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$spring$Bacteria,seq(0,120,12), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Reds 2")[2],  ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); title('Bacteria')
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))

p1 <- hist(alpha.diversity$summer$Fish,seq(-0.5,6.5,1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$summer$Fish,seq(-0.5,6.5,1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Blues 3")[2],  ylim=c(0,0.5), ylab="Frequency", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25),  add=T, main="")
axis(1, pos=0); axis(2,pos=-0.5,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$summer$Invert,seq(0,30,3), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$summer$Invert,seq(0,30,3), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Greens 3")[2], ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$summer$Bacteria,seq(0,120,12), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$summer$Bacteria,seq(0,120,12), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Reds 2")[2],  ylim=c(0,0.5), ylab="", xlab="",main="",axes=F);
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))

p1 <- hist(alpha.diversity$autumn$Fish,seq(-0.5,6.5,1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$autumn$Fish,seq(-0.5,6.5,1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Blues 3")[2],  ylim=c(0,0.5), ylab="Frequency", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25),  add=T, main="")
axis(1, pos=0); axis(2,pos=-0.5,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$autumn$Invert,seq(0,30,3), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$autumn$Invert,seq(0,30,3), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Greens 3")[2], ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(alpha.diversity$autumn$Bacteria,seq(0,120,12), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(alpha.div.data$autumn$Bacteria,seq(0,120,12), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Reds 2")[2],  ylim=c(0,0.5), ylab="", xlab="",main="",axes=F);
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
dev.off()

vecA <- seq(min(log10(catch$AG$A)), max(log10(catch$AG$A)), length.out=100)
vecA <- 10^vecA
summary.alpha <- data.frame(matrix("",3,3),row.names=group.names)
names(summary.alpha) <- time_id_vec
pdf(file="Fig2.pdf",width=20/2.54,height=18/2.54)
par(mfrow=c(3,3))
plot(catch$AG$A, alpha.diversity$spring$Fish, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,10),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Blues 3")[2])
points(A_data,alpha.div.data$spring$Fish, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Fish")
ll <- fit_lm_error(alpha.diversity$spring$Fish, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$spring$Fish ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,9,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["fish",]$spring <- ll$flag

plot(catch$AG$A, alpha.diversity$spring$Invert, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,30),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Greens 3")[2])
points(A_data,alpha.div.data$spring$Invert, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Invertebrates")
ll <- fit_lm_error(alpha.diversity$spring$Invert, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$spring$Invert ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,27,sprintf('%s;     R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["invert",]$spring <- ll$flag

plot(catch$AG$A, alpha.diversity$spring$Bacteria, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,120),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Reds 2")[2])
points(A_data,alpha.div.data$spring$Bacteria, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Bacteria")
ll <- fit_lm_error(alpha.diversity$spring$Bacteria, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$spring$Bacteria ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,110,sprintf('%s;    R^2: %.1e',ll$flag,ss1$r.squared))
summary.alpha["bact",]$spring <- ll$flag

plot(catch$AG$A, alpha.diversity$summer$Fish, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,10),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Blues 3")[2])
points(A_data,alpha.div.data$summer$Fish, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5);
ll <- fit_lm_error(alpha.diversity$summer$Fish, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$summer$Fish ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,9,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["fish",]$summer <- ll$flag

plot(catch$AG$A, alpha.diversity$summer$Invert, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,30),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Greens 3")[2])
points(A_data,alpha.div.data$summer$Invert, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(alpha.diversity$summer$Invert, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$summer$Invert ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,27,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["invert",]$summer <- ll$flag

plot(catch$AG$A, alpha.diversity$summer$Bacteria, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,120),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Reds 2")[2])
points(A_data,alpha.div.data$summer$Bacteria, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(alpha.diversity$summer$Bacteria, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$summer$Bacteria ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,110,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["bact",]$summer <- ll$flag

plot(catch$AG$A, alpha.diversity$autumn$Fish, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9), 
     ylim=c(0,10),bty="n",xaxt="n",yaxt="n",xlab="Drainage area [km2]",ylab="",col=hcl.colors(3,"Blues 3")[2])
points(A_data,alpha.div.data$autumn$Fish, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(alpha.diversity$autumn$Fish, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$autumn$Fish ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,9,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["fish",]$autumn <- ll$flag

plot(catch$AG$A, alpha.diversity$autumn$Invert, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,30),bty="n",xaxt="n",yaxt="n",xlab="Drainage area [km2]",ylab="",col=hcl.colors(3,"Greens 3")[2])
points(A_data,alpha.div.data$autumn$Invert, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(alpha.diversity$autumn$Invert, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$autumn$Invert ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,27,sprintf('%s;    R^2: %.3f',ll$flag,ss1$r.squared))
summary.alpha["invert",]$autumn <- ll$flag

plot(catch$AG$A, alpha.diversity$autumn$Bacteria, pch=19, cex=0.5, log="x",  xlim=c(1e5,1e9),
     ylim=c(0,120),bty="n",xaxt="n",yaxt="n",xlab="Drainage area [km2]",ylab="",col=hcl.colors(3,"Reds 2")[2])
points(A_data,alpha.div.data$autumn$Bacteria, cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(alpha.diversity$autumn$Bacteria, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(alpha.diversity$autumn$Bacteria ~ catch$AG$A); ss1 <- summary(lmod)
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=1)
text(5e6,110,sprintf('%s;    R^2: %.1e',ll$flag,ss1$r.squared))
summary.alpha["bact",]$autumn <- ll$flag
dev.off()


# Fig. 3 ####
pdf(file="Fig3.pdf",width=15/2.54, height=14/2.54)
par(mfrow=c(3,3),mar=c(0,0,0,0))
draw_thematic_catch(catch, alpha.diversity$spring$Fish, colLevels = c(0,6,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Blues 3",rev=T), backgroundColor = "#303030"); #title("spring")
draw_thematic_catch(catch, alpha.diversity$spring$Invert, colLevels = c(0,25,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Greens 3",rev=T), backgroundColor = "#303030"); #title("spring")
draw_thematic_catch(catch, alpha.diversity$spring$Bacteria, colLevels = c(0,100,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Reds 2",rev=T), backgroundColor = "#303030"); #title("spring")

draw_thematic_catch(catch, alpha.diversity$summer$Fish, colLevels = c(0,6,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Blues 3",rev=T), backgroundColor = "#303030"); #title("summer")
draw_thematic_catch(catch, alpha.diversity$summer$Invert, colLevels = c(0,25,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Greens 3",rev=T), backgroundColor = "#303030"); #title("summer")
draw_thematic_catch(catch, alpha.diversity$summer$Bacteria, colLevels = c(0,100,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Reds 2",rev=T), backgroundColor = "#303030"); #title("summer") 

draw_thematic_catch(catch, alpha.diversity$autumn$Fish, colLevels = c(0,6,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Blues 3",rev=T), backgroundColor = "#303030"); #title("autumn")
draw_thematic_catch(catch, alpha.diversity$autumn$Invert, colLevels = c(0,25,1000), addLegend = F,
                    colPalette = hcl.colors(1000,"Greens 3",rev=T), backgroundColor = "#303030"); #title("autumn")
draw_thematic_catch(catch, alpha.diversity$autumn$Bacteria, colLevels = c(0,100,1000), addLegend = F, 
                    colPalette = hcl.colors(1000,"Reds 2",rev=T), backgroundColor = "#303030"); #title("autumn")

dev.off()


# Fig. 4 ####

summary.beta.spat <- data.frame(matrix("",3,3),row.names=group.names)
names(summary.beta.spat) <- time_id_vec
pdf(file="Fig4.pdf",width=20/2.54,height=8/2.54)
par(mfrow=c(1,3))
for (time in time_id_vec){
  mm <- cbind(beta_HeadMain$fish[[time]]$mean_eff, beta_HeadMain$invert[[time]]$mean_eff, beta_HeadMain$bact[[time]]$mean_eff)
  for (i in 1:6){ qq <- quantile(mm[,i],c(0.025,0.975));  mm[mm[,i]<qq[1],i] <- qq[1]; mm[mm[,i]>qq[2],i] <- qq[2]}
  for (j in 1:3){
    if (min(mm[,(j-1)*2+1])>max(mm[,2*j])){
      summary.beta.spat[4-j,][[time]] <- "-"
    } else if (min(mm[,2*j])>max(mm[,(j-1)*2+1])){
      summary.beta.spat[4-j,][[time]] <- "+"
    } else {summary.beta.spat[4-j,][[time]] <- "="}
  }
  boxplot(mm, range=0,ylim=c(0.4, 0.8), 
          col=c(hcl.colors(4,"Blues 3")[2:3], hcl.colors(4,"Greens 3")[2:3], hcl.colors(4,"Reds 2")[2:3]),
          ylab="Jaccard index"); title(time); abline(v=c(2.5,4.5), col="gray")
}
dev.off()

# Fig. 5 ####
seas.names <- c("spring_summer","summer_autumn","spring_autumn")
group.names <- c("bact","invert","fish")
summary.beta.temp <- data.frame(matrix("",3,3),row.names=group.names)
names(summary.beta.temp) <- seas.names
pdf(file="Fig5.pdf",width=20/2.54,height=14/2.54)
vecA <- seq(min(log10(catch$AG$A)), max(log10(catch$AG$A)), length.out=100)
vecA <- 10^vecA
par(mfrow=c(2,3))
plot(catch$AG$A, beta_temp$fish$jac$spring_summer, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Blues 3")[2])
points(A_data,beta_temp.data$fish$spring_summer,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Fish")
ll <- fit_lm_error(beta_temp$fish$jac$spring_summer, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$fish$jac$spring_summer ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.3f',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["fish",]$spring_summer <- ll$flag

plot(catch$AG$A, beta_temp$invert$jac$spring_summer, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Greens 3")[2])
points(A_data,beta_temp.data$invert$spring_summer,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Invertebrates")
ll <- fit_lm_error(beta_temp$invert$jac$spring_summer, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$invert$jac$spring_summer ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.3f',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["invert",]$spring_summer <- ll$flag

plot(catch$AG$A, beta_temp$bact$jac$spring_summer, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Reds 2")[2])
points(A_data,beta_temp.data$bact$spring_summer,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); title("Bacteria")
ll <- fit_lm_error(beta_temp$bact$jac$spring_summer, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$bact$jac$spring_summer ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.3f',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["bact",]$spring_summer <- ll$flag

plot(catch$AG$A, beta_temp$fish$jac$spring_autumn, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Blues 3")[2])
points(A_data,beta_temp.data$fish$spring_autumn,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5); 
ll <- fit_lm_error(beta_temp$fish$jac$spring_autumn, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$fish$jac$spring_autumn ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.3f',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["fish",]$spring_autumn <- ll$flag

plot(catch$AG$A, beta_temp$invert$jac$spring_autumn, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Greens 3")[2])
points(A_data,beta_temp.data$invert$spring_autumn,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5);
ll <- fit_lm_error(beta_temp$invert$jac$spring_autumn, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$invert$jac$spring_autumn ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.3f',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["invert",]$spring_autumn <- ll$flag

plot(catch$AG$A, beta_temp$bact$jac$spring_autumn, pch=19, cex=0.5, log="x", xlim=c(1e5,1e9),
     ylim=c(0,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="",col=hcl.colors(3,"Reds 2")[2])
points(A_data,beta_temp.data$bact$spring_autumn,cex=1.25)
axis(1,pos=0); axis(2,pos=1e5);
ll <- fit_lm_error(beta_temp$bact$jac$spring_autumn, catch$AG$A)
polygon(c(vecA, rev(vecA)), c(ll$QQ[1,], rev(ll$QQ[2,])), col="#BBBBBB", border=NA)
lmod <- lm(beta_temp$bact$jac$spring_autumn ~ catch$AG$A); ss1 <- summary(lmod)
text(5e6,0.1,sprintf('p: %.1e; - R^2: %.1e',ss1$coefficients[2,4],ss1$r.squared))
lines(vecA, ss1$coefficients[1,1]+ss1$coefficients[2,1]*vecA, lwd=2)
summary.beta.temp["bact",]$spring_autumn <- ll$flag
dev.off()

pdf(file="Fig5_inset.pdf",width=20/2.54,height=10/2.54)
par(mfrow=c(2,3))
p1 <- hist(beta_temp$fish$jac$spring_summer,seq(0,1,0.1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$fish$spring_summer,seq(0,1,0.1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Blues 3")[2],  ylim=c(0,0.5), ylab="Frequency", xlab="",main="",axes=F); title('Fish')
plot(p2,col=rgb(0,0,0,0.25),  add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(beta_temp$invert$jac$spring_summer,seq(0,1,0.1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$invert$spring_summer,seq(0,1,0.1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Greens 3")[2], ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); title('Invertebrates')
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(beta_temp$bact$jac$spring_summer,seq(0,1,0.1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$bact$spring_summer,seq(0,1,0.1),  plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Reds 2")[2],  ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); title('Bacteria')
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))

p1 <- hist(beta_temp$fish$jac$spring_autumn,seq(0,1,0.1),  plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$fish$spring_autumn,seq(0,1,0.1),  plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Blues 3")[2],  ylim=c(0,0.6), ylab="Frequency", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25),  add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.3,0.6))
p1 <- hist(beta_temp$invert$jac$spring_autumn,seq(0,1,0.1), plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$invert$spring_autumn,seq(0,1,0.1), plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Greens 3")[2], ylim=c(0,0.5), ylab="", xlab="",main="",axes=F); 
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))
p1 <- hist(beta_temp$bact$jac$spring_autumn,seq(0,1,0.1),  plot=F); p1$counts <- p1$counts/sum(p1$counts)
p2 <- hist(beta_temp.data$bact$spring_autumn,seq(0,1,0.1),  plot=F); p2$counts <- p2$counts/sum(p2$counts)
plot(p1, col=hcl.colors(3,"Reds 2")[2],  ylim=c(0,0.5), ylab="", xlab="",main="",axes=F);
plot(p2,col=rgb(0,0,0,0.25), add=T, main="")
axis(1, pos=0); axis(2,pos=0,at=c(0,0.25,0.5))

dev.off()

# Fig. 6 ####
pdf(file="Fig6.pdf",width=15/2.54,height=14/2.54)
par(mfrow=c(3,3),mar=c(0,0,0,0))
draw_thematic_catch(catch, beta_temp$fish$jac$spring_summer, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Blues 3",rev=T),nanColor="#999999",addLegend=F); #title("Fish")
draw_thematic_catch(catch, beta_temp$invert$jac$spring_summer, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Greens 3",rev=T),nanColor="#999999",addLegend=F); #title("Invertebrates")
draw_thematic_catch(catch, beta_temp$bact$jac$spring_summer, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Reds 2",rev=T),nanColor="#999999",addLegend=F); #title("Bacteria")
draw_thematic_catch(catch, beta_temp$fish$jac$spring_autumn, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Blues 3",rev=T), nanColor="#999999",addLegend=F); # 
draw_thematic_catch(catch, beta_temp$invert$jac$spring_autumn, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Greens 3",rev=T),nanColor="#999999",addLegend=F); 
draw_thematic_catch(catch, beta_temp$bact$jac$spring_autumn, colLevels = c(0,1,1000),
                    colPalette=hcl.colors(1000,"Reds 2",rev=T),nanColor="#999999",addLegend=F)
dev.off()


# Fig. S1 ####
# re-read data_COI (as Barbus, Gobio, Phoxinus had been removed)
tmp <- read.csv(paste0(pathname,"/calibration/COI/COI.csv"))

names(tmp)[1] <- "Genus"
genusNames <- substr(tmp$Genus,4,100)
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
  genus <- c(genus, rep(genusNames[g], length(tmp$time_id)))
  value <- c(value, tmp[[genusNames[g]]])
  time_id <- c(time_id, tmp$time_id)
  site_id <- c(site_id, tmp$site_id)
}
rm(tmp)

data_COI=data.frame(matrix(nrow=length(value),ncol=0))
data_COI[['Genus']] <- genus
data_COI[['time_id']] <- time_id
data_COI[['site_id']] <- site_id
data_COI[['value']] <- value

sumReads.bact <- c(sum(data_16S$value[which(funct.df$description[match(data_16S$Genus,funct.df$Taxa)]=="Bacteria" & data_16S$time_id=="spring")]),
         sum(data_16S$value[which(funct.df$description[match(data_16S$Genus,funct.df$Taxa)]=="Bacteria" & data_16S$time_id=="summer")]),
         sum(data_16S$value[which(funct.df$description[match(data_16S$Genus,funct.df$Taxa)]=="Bacteria" & data_16S$time_id=="autumn")]))
sumReads.invert <- c(sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Invert" & data_COI$time_id=="spring")]),
         sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Invert" & data_COI$time_id=="summer")]),
         sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Invert" & data_COI$time_id=="autumn")]))
sumReads.fish_COI <- c(sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Fish" & data_COI$time_id=="spring")]),
                     sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Fish" & data_COI$time_id=="summer")]),
                     sum(data_COI$value[which(funct.df$description[match(data_COI$Genus,funct.df$Taxa)]=="Fish" & data_COI$time_id=="autumn")]))
sumReads.fish_12S <- c(sum(data_12S$value[which(funct.df$description[match(data_12S$Genus,funct.df$Taxa)]=="Fish" & data_12S$time_id=="spring")]),
                     sum(data_12S$value[which(funct.df$description[match(data_12S$Genus,funct.df$Taxa)]=="Fish" & data_12S$time_id=="summer")]),
                     sum(data_12S$value[which(funct.df$description[match(data_12S$Genus,funct.df$Taxa)]=="Fish" & data_12S$time_id=="autumn")]))

bacteria_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Bacteria"]
sumPresence.bact <- c(0,0,0)
for (t in 1:3){
  for (b in 1:length(bacteria_names)){
    sumPresence.bact[t] <- sumPresence.bact[t] + (max(data_16S$value[which(data_16S$Genus==bacteria_names[b] & data_16S$time_id==time_id_vec[t])])>0) 
  }
}
invert_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Invert"]
sumPresence.invert <- c(0,0,0)
for (t in 1:3){
  for (b in 1:length(invert_names)){
    sumPresence.invert[t] <- sumPresence.invert[t] + (max(data_COI$value[which(data_COI$Genus==invert_names[b] & data_COI$time_id==time_id_vec[t])])>0) 
  }
}
fish_names <- allGenusNames[funct.df$description[match(allGenusNames,funct.df$Taxa)]=="Fish"]
sumPresence.fish_COI <- c(0,0,0)
for (t in 1:3){
  for (b in 1:length(fish_names)){
    sumPresence.fish_COI[t] <- suppressWarnings(sumPresence.fish_COI[t] + (max(data_COI$value[which(data_COI$Genus==fish_names[b] & data_COI$time_id==time_id_vec[t])])>0))
  }
}
sumPresence.fish_12S <- c(0,0,0)
for (t in 1:3){
  for (b in 1:length(fish_names)){
sumPresence.fish_12S[t] <- suppressWarnings(sumPresence.fish_12S[t] + (max(data_12S$value[which(data_12S$Genus==fish_names[b] & data_12S$time_id==time_id_vec[t])])>0)) 
  }
}
 
pdf(file="FigS1.pdf",width=18/2.54, height=12/2.54) 
par(mfrow=c(2,4))
barplot(sumReads.fish_12S,ylim=c(0,3e5), col=seas.cols, ylab="No. reads"); title('Fish - 12S')
barplot(sumReads.fish_COI, ylim=c(0,6e5), col=seas.cols); title('Fish - COI')
barplot(sumReads.invert, ylim=c(0,4e6), col=seas.cols); title('Invertebrates')
barplot(sumReads.bact, ylim=c(0,6e5), col=seas.cols); title('Bacteria')
barplot(sumPresence.fish_12S, ylim=c(0,10), col=seas.cols, ylab="No. detected genera", names.arg=time_id_vec);
barplot(sumPresence.fish_COI, ylim=c(0,10), col=seas.cols, names.arg=time_id_vec);
barplot(sumPresence.invert, ylim=c(0,80), col=seas.cols, names.arg=time_id_vec);
barplot(sumPresence.bact, ylim=c(0,250), col=seas.cols, names.arg=time_id_vec); 
dev.off()


# Fig. S2 ####
pdf(file="FigS2.pdf",height=22/2.54, width=16/2.54)
names(Zcovariates) <- c("M-DA","M-SO","M-LS","M-US","M-LE",
                        "G-AL","G-AP","G-LO","G-MO","G-PE","G-SC","G-WA",
                        "L-FO","L-LA","L-OR","L-RO","L-SW","L-UR",
                        "LUT","GON","WIS","NE3","GL2","GL1","DIE",
                        "TH8","TH7","TH6","TH5","NE2","TH4","NE1","TH3","TH2","TH1")
par(mfrow=c(3,3)); 
data <- covariate_significance$spring[,match(funct.df$Taxa[which(funct.df$description=="Fish")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-6,6)); title('Fish - Spring')
data <- covariate_significance$spring[,match(funct.df$Taxa[which(funct.df$description=="Invert")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-20,20)); title('Invert - Spring')
data <- covariate_significance$spring[,match(funct.df$Taxa[which(funct.df$description=="Bacteria")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-40,40)); title('Bacteria - Spring')

data <- covariate_significance$summer[,match(funct.df$Taxa[which(funct.df$description=="Fish")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-6,6)); title('Fish - Summer')
data <- covariate_significance$summer[,match(funct.df$Taxa[which(funct.df$description=="Invert")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-20,20)); title('Invert - Summer')
data <- covariate_significance$summer[,match(funct.df$Taxa[which(funct.df$description=="Bacteria")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-40,40)); title('Bacteria - Summer')

data <- covariate_significance$autumn[,match(funct.df$Taxa[which(funct.df$description=="Fish")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-6,6)); title('Fish - Autumn')
data <- covariate_significance$autumn[,match(funct.df$Taxa[which(funct.df$description=="Invert")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-20,20)); title('Invert - Autumn') 
data <- covariate_significance$autumn[,match(funct.df$Taxa[which(funct.df$description=="Bacteria")],allGenusNames)]
tornado.plot(data, Zcovariates,c(-40,40)); title('Bacteria - Autumn')
dev.off()

# Fig. S3 ####
pdf(file="FigS3.pdf",width=20/2.54,height=14/2.54)
par(mfrow=c(2,3))
plot(alpha.diversity$spring$Bacteria,alpha.diversity$spring$Invert,
     xlim=c(0,120),ylim=c(0,25), pch=19, col=seas.cols[1], cex=0.5, bty="n",xaxt="n",yaxt="n", 
     xlab="Bacteria richness",ylab="Invertebrate richness")
axis(1,pos=0); axis(2,pos=0); title("Spring")
lmod <- lm(alpha.diversity$spring$Invert ~ alpha.diversity$spring$Bacteria)
lines(c(0,120), lmod$coefficients[2]*c(0,120) + lmod$coefficients[1], col="black", lwd=2)
text(40,25,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(30,23.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(30,22,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$summer$Bacteria,alpha.diversity$summer$Invert, 
     xlim=c(0,120),ylim=c(0,25), pch=19, col=seas.cols[2], cex=0.5, bty="n",xaxt="n",yaxt="n",xlab="Bacteria richness",ylab="")
axis(1,pos=0); axis(2,pos=0); title("Summer")
lmod <- lm(alpha.diversity$summer$Invert ~ alpha.diversity$summer$Bacteria)
lines(c(0,120), lmod$coefficients[2]*c(0,120) + lmod$coefficients[1], col="black", lwd=2)
text(40,25,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(30,23.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(30,22,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$autumn$Bacteria,alpha.diversity$autumn$Invert, 
     xlim=c(0,120),ylim=c(0,25), pch=19, col=seas.cols[3], cex=0.5, bty="n",xaxt="n",yaxt="n",xlab="Bacteria richness",ylab="")
axis(1,pos=0); axis(2,pos=0); title("Autumn")
lmod <- lm(alpha.diversity$autumn$Invert ~ alpha.diversity$autumn$Bacteria)
lines(c(0,120), lmod$coefficients[2]*c(0,120) + lmod$coefficients[1], col="black", lwd=2)
text(40,25,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(30,23.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(30,22,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$spring$Invert,alpha.diversity$spring$Fish,
     xlim=c(0,25),ylim=c(0,10), pch=19, col=seas.cols[1], bty="n",xaxt="n",yaxt="n", cex=0.5,
     xlab="Invert richness", ylab="Fish richness")
axis(1,pos=0); axis(2,pos=0); 
lmod <- lm(alpha.diversity$spring$Fish ~ alpha.diversity$spring$Invert)
lines(c(0,25), lmod$coefficients[2]*c(0,25) + lmod$coefficients[1], col="black", lwd=2)
text(10,10,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(5,9.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(5,9,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$summer$Invert,alpha.diversity$summer$Fish,
     xlim=c(0,25),ylim=c(0,10), pch=19, col=seas.cols[2], bty="n",xaxt="n",yaxt="n", cex=0.5,
     xlab="Invert richness",ylab="")
axis(1,pos=0); axis(2,pos=0);
lmod <- lm(alpha.diversity$summer$Fish ~ alpha.diversity$summer$Invert)
lines(c(0,25), lmod$coefficients[2]*c(0,25) + lmod$coefficients[1], col="black", lwd=2)
text(10,10,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(5,9.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(5,9,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$autumn$Invert,alpha.diversity$autumn$Fish,
     xlim=c(0,25),ylim=c(0,10), pch=19, col=seas.cols[3], bty="n",xaxt="n",yaxt="n", cex=0.5,
     xlab="Invert richness",ylab="")
axis(1,pos=0); axis(2,pos=0)
lmod <- lm(alpha.diversity$autumn$Fish ~ alpha.diversity$autumn$Invert)
lines(c(0,25), lmod$coefficients[2]*c(0,25) + lmod$coefficients[1], col="black", lwd=2)
text(10,10,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(5,9.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(5,9,sprintf("R^2 = %.3f",summary(lmod)$r.squared))
dev.off()

# Fig. S4 ####
pdf(file="FigS4.pdf",width=20/2.54,height=14/2.54)
par(mfrow=c(2,3))

plot(alpha.diversity$spring$Fish,alpha.diversity$summer$Fish, cex=0.5,
     xlim=c(0,10),ylim=c(0,10), pch=19, col=hcl.colors(3,"Blues 3")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in summer")
axis(1,pos=0); axis(2,pos=0); title("Fish")
lmod <- lm(alpha.diversity$summer$Fish ~ alpha.diversity$spring$Fish)
lines(c(0,10), lmod$coefficients[2]*c(0,10) + lmod$coefficients[1], col="black", lwd=2)
text(3,10,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(2,9,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(2,8,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$spring$Invert,alpha.diversity$summer$Invert, cex=0.5,
     xlim=c(0,25),ylim=c(0,25), pch=19, col=hcl.colors(3,"Greens 3")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in summer")
axis(1,pos=0); axis(2,pos=0); title("Invert")
lmod <- lm(alpha.diversity$summer$Invert ~ alpha.diversity$spring$Invert)
lines(c(0,25), lmod$coefficients[2]*c(0,25) + lmod$coefficients[1], col="black", lwd=2)
text(8,25,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(7,23.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(6,22,sprintf("R^2 = %.3f",summary(lmod)$r.squared))


plot(alpha.diversity$spring$Bacteria,alpha.diversity$summer$Bacteria, cex=0.5,
     xlim=c(0,120),ylim=c(0,120), pch=19, col=hcl.colors(3,"Reds 2")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in summer")
axis(1,pos=0); axis(2,pos=0); title("Bacteria")
lmod <- lm(alpha.diversity$summer$Bacteria ~ alpha.diversity$spring$Bacteria)
lines(c(0,120), lmod$coefficients[2]*c(0,120) + lmod$coefficients[1], col="black", lwd=2)
text(40,120,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(30,110,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(30,100,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$spring$Fish,alpha.diversity$autumn$Fish, cex=0.5,
     xlim=c(0,10),ylim=c(0,10), pch=19, col=hcl.colors(3,"Blues 3")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in autumn")
axis(1,pos=0); axis(2,pos=0); 
lmod <- lm(alpha.diversity$autumn$Fish ~ alpha.diversity$spring$Fish)
lines(c(0,10), lmod$coefficients[2]*c(0,10) + lmod$coefficients[1], col="black", lwd=2)
text(3,10,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(2,9,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(2,8,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$spring$Invert,alpha.diversity$autumn$Invert, cex=0.5,
     xlim=c(0,25),ylim=c(0,25), pch=19, col=hcl.colors(3,"Greens 3")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in autumn")
axis(1,pos=0); axis(2,pos=0); 
lmod <- lm(alpha.diversity$autumn$Invert ~ alpha.diversity$spring$Invert)
lines(c(0,25), lmod$coefficients[2]*c(0,25) + lmod$coefficients[1], col="black", lwd=2)
text(8,25,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(7,23.5,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(6,22,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

plot(alpha.diversity$spring$Bacteria,alpha.diversity$autumn$Bacteria, cex=0.5,
     xlim=c(0,120),ylim=c(0,120), pch=19, col=hcl.colors(3,"Reds 2")[2], bty="n",xaxt="n",yaxt="n", 
     xlab="Richness in spring",ylab="Richness in autumn")
axis(1,pos=0); axis(2,pos=0); 
lmod <- lm(alpha.diversity$autumn$Bacteria ~ alpha.diversity$spring$Bacteria)
lines(c(0,120), lmod$coefficients[2]*c(0,120) + lmod$coefficients[1], col="black", lwd=2)
text(40,120,sprintf("y = %.2f x + %.2f",lmod$coefficients[2],lmod$coefficients[1]))
text(30,110,sprintf("p = %.1e",summary(lmod)$coefficients[2,4]))
text(30,100,sprintf("R^2 = %.3f",summary(lmod)$r.squared))

dev.off()

# Fig. S5 ####
pdf(file="FigS5.pdf",width=18/2.54,height=10/2.54)
par(mfrow=c(1,2))
plot((1:length(fish_names))/length(fish_names),sort(colSums(PA_matrix$spring[,fish_names]),decreasing=T)/catch$AG$nNodes, 
     type="l", ylim=c(0,1), xlim=c(0,1),xaxt="n",yaxt="n",bty="n",
     col=hcl.colors(3,"Blues 3")[2],xlab="Fraction of taxa",ylab="Fraction of reaches occupied")
axis(1,pos=0); axis(2,pos=0)
lines((1:length(invert_names))/length(invert_names),sort(colSums(PA_matrix$spring[,invert_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Greens 3")[2])
lines((1:length(bacteria_names))/length(bacteria_names),sort(colSums(PA_matrix$spring[,bacteria_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Reds 2")[2]); #title('Spring')
lines((1:length(fish_names))/length(fish_names),sort(colSums(PA_matrix$summer[,fish_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Blues 3")[2],xlim=c(0,1),lty=2)
lines((1:length(invert_names))/length(invert_names),sort(colSums(PA_matrix$summer[,invert_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Greens 3")[2],xlim=c(0,1),lty=2)
lines((1:length(bacteria_names))/length(bacteria_names),sort(colSums(PA_matrix$summer[,bacteria_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Reds 2")[2],lty=2); #title('Summer')
lines((1:length(fish_names))/length(fish_names),sort(colSums(PA_matrix$autumn[,fish_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Blues 3")[2],xlim=c(0,1),lty=4)
lines((1:length(invert_names))/length(invert_names),sort(colSums(PA_matrix$autumn[,invert_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Greens 3")[2], xlim=c(0,1),lty=4)
lines((1:length(bacteria_names))/length(bacteria_names),sort(colSums(PA_matrix$autumn[,bacteria_names]),decreasing=T)/catch$AG$nNodes, 
      col=hcl.colors(3,"Reds 2")[2],lty=4); #title('Autumn')
legend("topright", legend=c("Fish spring","Fish summer","Fish autumn","Invert spring","Invert summer","Invert autumn",
                            "Bacteria spring","Bacteria summer","Bacteria autumn"), 
       lty=c(1,2,4,1,2,4,1,2,4), col=c(rep(hcl.colors(3,"Blues 3")[2],3),rep(hcl.colors(3,"Greens 3")[2],3),
                                       rep(hcl.colors(3,"Reds 2")[2],3)))

plot((1:length(fish_names))/length(fish_names),
     sort(colSums((PA_matrix$spring[,fish_names]+PA_matrix$summer[,fish_names]+
                     PA_matrix$autumn[,fish_names])==3),decreasing = T)/catch$AG$nNodes,
     type="l", ylim=c(0,1), xlim=c(0,1),xaxt="n",yaxt="n",bty="n",
     col=hcl.colors(3,"Blues 3")[2], xlab="Fraction of taxa",ylab="Fraction of reaches occupied")
axis(1,pos=0); axis(2,pos=0)
lines((1:length(invert_names))/length(invert_names),
      sort(colSums((PA_matrix$spring[,invert_names]+PA_matrix$summer[,invert_names]+
                      PA_matrix$autumn[,invert_names])==3),decreasing = T)/catch$AG$nNodes,
      col=hcl.colors(3,"Greens 3")[2])
lines((1:length(bacteria_names))/length(bacteria_names),
      sort(colSums((PA_matrix$spring[,bacteria_names]+PA_matrix$summer[,bacteria_names]+
                      PA_matrix$autumn[,bacteria_names])==3),decreasing = T)/catch$AG$nNodes,
      col=hcl.colors(3,"Reds 2")[2])
dev.off()

# Fig. S6 ####
pdf(file="FigS6.pdf",width=18/2.54, height=21/2.54) 
par(mfrow=c(3,3), mai=c(0,0,0,0))
draw_thematic_catch(catch, PA_matrix$spring$Gobio, discreteLevels = T, colPalette = c('white','red')); title('Spring')
draw_thematic_catch(catch, PA_matrix$summer$Gobio, discreteLevels = T, colPalette = c('white','red')); title('Summer')
draw_thematic_catch(catch, PA_matrix$autumn$Gobio, discreteLevels = T, colPalette = c('white','red')); title('Autumn')
draw_thematic_catch(catch, PA_matrix$spring$Phoxinus, discreteLevels = T, colPalette = c('white','red')); 
draw_thematic_catch(catch, PA_matrix$summer$Phoxinus, discreteLevels = T, colPalette = c('white','red')); 
draw_thematic_catch(catch, PA_matrix$autumn$Phoxinus, discreteLevels = T, colPalette = c('white','red')); 
draw_thematic_catch(catch, PA_matrix$spring$Salmo, discreteLevels = T, colPalette = c('white','red')); 
draw_thematic_catch(catch, PA_matrix$summer$Salmo, discreteLevels = T, colPalette = c('white','red')); 
draw_thematic_catch(catch, PA_matrix$autumn$Salmo, discreteLevels = T, colPalette = c('white','red')); 
dev.off()

betaspat_temp <- data.frame(matrix(0,3,3), row.names=c('Fish','Invert','Bacteria'))
names(betaspat_temp) <- time_id_vec
for (time in time_id_vec){
  betaspat_temp['Fish', time] <- mean(beta.div.fish[[time]]$beta.jac,na.rm=T)
  betaspat_temp['Invert', time] <- mean(beta.div.invert[[time]]$beta.jac,na.rm=T)
  betaspat_temp['Bacteria', time] <- mean(beta.div.bact[[time]]$beta.jac,na.rm=T)
}

betaspat_partitioning <- data.frame(matrix(0,3,3), row.names=c('Fish','Invert','Bacteria'))
names(betaspat_partitioning) <- time_id_vec
for (time in time_id_vec){
betaspat_partitioning['Fish', time] <- mean(beta.div.fish[[time]]$beta.jtu/beta.div.fish[[time]]$beta.jac,na.rm=T)
betaspat_partitioning['Invert', time] <- mean(beta.div.invert[[time]]$beta.jtu/beta.div.invert[[time]]$beta.jac,na.rm=T)
betaspat_partitioning['Bacteria', time] <- mean(beta.div.bact[[time]]$beta.jtu/beta.div.bact[[time]]$beta.jac,na.rm=T)
}

betatemp_mean <- data.frame(matrix(0,3,2), row.names=c('Fish','Invert','Bacteria'))
names(betatemp_mean) <- c('spring_summer','spring_autumn')
for (time in names(betatemp_mean)){
  betatemp_mean['Fish',time] <- mean(beta_temp$fish$jac[[time]],na.rm=T)
  betatemp_mean['Invert',time] <- mean(beta_temp$invert$jac[[time]],na.rm=T)
  betatemp_mean['Bacteria',time] <- mean(beta_temp$bact$jac[[time]],na.rm=T)
}

betatemp_partitioning <- data.frame(matrix(0,3,2), row.names=c('Fish','Invert','Bacteria'))
names(betatemp_partitioning) <- c('spring_summer','spring_autumn')
for (time in names(betatemp_mean)){
  betatemp_partitioning['Fish',time] <- mean(beta_temp$fish$jtu[[time]]/beta_temp$fish$jac[[time]],na.rm=T)
  betatemp_partitioning['Invert',time] <- mean(beta_temp$invert$jtu[[time]]/beta_temp$invert$jac[[time]],na.rm=T)
  betatemp_partitioning['Bacteria',time] <- mean(beta_temp$bact$jtu[[time]]/beta_temp$bact$jac[[time]],na.rm=T)
}


# Fig. S7 ####
pdf(file="FigS7.pdf",width=18/2.54,height=7/2.54)
par(mfrow=c(1,2))
hist(catch$AG$streamOrder[catch$AG$A<=median(catch$AG$A)],seq(0.5,5.5,1),ylim=c(0,1000),
     xlab="Stream Order",xlim=c(0.5,5.5),xaxt="n",yaxt="n", main='Upstream' )
axis(1,pos=0); axis(2,pos=0.5)
hist(catch$AG$streamOrder[catch$AG$A>median(catch$AG$A)],seq(0.5,5.5,1),ylim=c(0,1000),
     xlab="Stream Order",xlim=c(0.5,5.5),xaxt="n",yaxt="n", main='Downstream',ylab="" )
axis(1,pos=0); axis(2,pos=0.5)
dev.off()



