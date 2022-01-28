eval_hydrology <- function(catch,sites,filename){
  hydrodata <- read.csv(filename)
  
  hydrosites_id <- unique(hydrodata$site_id)
  time_id <- unique(hydrodata$time_id)
  
  Area <- numeric(length(hydrosites_id))
  for (i in 1:length(hydrosites_id)){
    Area[i] <- catch$RN$A[sites$nodeRN[which(sites$Site_id==hydrosites_id[i])]]
  }
  
  width <- Q <- depth <- velocity <- matrix(0,catch$AG$nNodes,length(time_id))
  
  for (t in 1:length(time_id)){
    width_meas <- extract_sites(hydrodata,"W",time_id[t])
    lmod <- lm(log(width_meas) ~ log(Area))
    width[,t] <- exp(lmod$coefficients[1])*((catch$AG$A+catch$AG$AReach)/2)^lmod$coefficients[2]
    
    Q_meas <- extract_sites(hydrodata,"Q",time_id[t])
    lmod <- lm(log(Q_meas) ~ log(Area))
    Q[,t] <- exp(lmod$coefficients[1])*((catch$AG$A+catch$AG$AReach)/2)^lmod$coefficients[2]
    
    D_meas <- extract_sites(hydrodata,"D",time_id[t],"H2414")
    lmod <- lm(log(D_meas) ~ log(Area))
    depth[,t] <- exp(lmod$coefficients[1])*((catch$AG$A+catch$AG$AReach)/2)^lmod$coefficients[2]
  }
  velocity <- Q/width/depth 
  
  width_df <- as.data.frame(width); colnames(width_df) <- time_id
  Q_df <- as.data.frame(Q); colnames(Q_df) <- time_id
  depth_df <- as.data.frame(depth); colnames(depth_df) <- time_id
  velocity_df <- as.data.frame(velocity); colnames(velocity_df) <- time_id
  
  # calculate path velocities
  catch$AG$pathvelocities <- vector("list")
  for (t in 1:length(time_id)){
    sprintf("Calculate path velocities for time_id %s...",time_id[t])
    catch$AG$pathvelocities[[time_id[t]]] <- spam(0,catch$AG$nNodes,catch$AG$nNodes)
    set_nodes <- matrix(0,catch$AG$nNodes^2,2)
    set_values <- numeric(catch$AG$nNodes^2)
    k <- 1
    for (i in 1:catch$AG$nNodes){
      for (j in 1:catch$AG$nNodes){
        path <- catch$AG$downstreamPath[[i]][[j]]
        if (!is.null(path) && !(i == which(catch$AG$downNode==0) && j == which(catch$AG$downNode==0))){
          set_values[k] <- catch$AG$downstreamPathLength[i,j] / (sum(catch$AG$leng[path] / velocity_df[path,t]))
          set_nodes[k,] <- c(i,j) 
          k <- k + 1
        } else if (i == which(catch$AG$downNode==0) && j == which(catch$AG$downNode==0)) {
          set_values[k] <- velocity_df[i,t]
          set_nodes[k,] <- c(i,j) 
          k <- k + 1
        }
      }
    }
    set_values <- set_values[-(k:catch$AG$nNodes^2)]
    set_nodes <- set_nodes[-(k:catch$AG$nNodes^2),]
    catch$AG$pathvelocities[[time_id[t]]][set_nodes] <- set_values
    for (i in 1:catch$AG$nNodes){
      catch$AG$pathvelocities[[time_id[t]]][i,i] <- velocity_df[i,t] # patch to correct when length of path is null
    }
    
  }
  
  catch$AG[["width"]] <- width_df
  catch$AG[["Q"]] <- Q_df
  catch$AG[["depth"]] <- depth_df
  catch$AG[["velocity"]] <- velocity_df
  
  invisible(catch)
}


extract_sites <- function(hydrodata,type,
                          time_id=hydrodata$time_id,
                          discard_site=""){
  
  hydrosites_id <- unique(hydrodata$site_id)
  
  discarded <- !logical(length(hydrodata$site_id))
  for (i in 1:length(discard_site)){
    discarded <- discarded & hydrodata$site_id!=discard_site[i] 
  }
  
  indices <- which(hydrodata$type==type & hydrodata$time_id==time_id & discarded) 
  values <- hydrodata$value[indices] 
  
  kkk <- NA*numeric(length(hydrosites_id))
  for (i in 1:length(hydrosites_id)){
    tmp <- which(hydrodata$site_id[indices]==hydrosites_id[i])
    if (length(tmp)>0){
      kkk[i] <- values[tmp]
    }
  }
  invisible(kkk)
}