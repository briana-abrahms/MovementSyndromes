#' Calculation of Animal Movement Statistics #######
#' From: Abrahms et al. 2017, Movement Ecology #####
#' DOI: 10.1186/s40462-017-0104-2 ##################

#' Calculate 5 metrics:
#' 1. Turn angle correlation (TAC)
#' 2. Mean Residence Time (RT)
#' 3. Mean Time to Return (T2R)
#' 4. Mean Volume of Intersection (VI)
#' 5. Max Net Squared Displacement (MNSD)

################################ SETUP ###############################

# load libraries
library(adehabitatLT)
library(adehabitatHR)
library(lubridate)

# Convert dataframe to ltraj object
ltraj <- as.ltraj(xy=Data[,c("lon", "lat")], 
                  date=as.POSIXct(Data$Timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
                  id=Data$id)


############################## CALCULATE METRICS ############################

#################################
# 1. Turn angle correlation (TAC) 
#################################

# for-loop for calculating TAC for multiple individuals
# NB: a) all trajectories in the ltraj object must be 'regular'

TAC<-matrix(ncol=1, nrow=length(ltraj)) # create empty data frame to populate with for-loop

for (i in 1:length(ltraj)){
  SA <- adehabitatLT::acfang.ltraj(ltraj[i], which = "relative") 
  TAC[i,] <- 1/(SA[[1]][1,])
}

View(TAC)


######################################################
# 2-3. Residence Times (RT) and Times to Return (T2R)
######################################################

#' Residence time = the number of hours the animal spends inside a circle of a given radius 
#' centered on each location without leaving the radius for more than a specified cut-off time
#' Time-to-return = the number of hours the animal spends beyond a specified cut-off time 
#' before its return to a circle of a given radius centered on each location
#' *Adapted from van Moorter et al. 2015, Journal of Animal Ecology

# define Residence Times and Times to Return functions
RTandT2R <- function(x, radius, maxt, units="hour", addinfo = F){
  fR <- function(x, dframe, radius, maxt, units=units){
    tmp <- dframe[c(x:nrow(dframe)),]
    dists <- sqrt((tmp$x - tmp$x[1])^2 + (tmp$y - tmp$y[1])^2)
    dists <- as.numeric(dists<=radius)
    ext <- which(dists[-length(dists)] > dists[-1])+1
    entr <-  which(dists[-length(dists)] < dists[-1])+1
    bts <- difftime(tmp$date[entr], tmp$date[ext[c(1:length(entr))]], units=units)    
    tmp1 <- as.numeric(difftime(tmp$date[ext[(as.numeric(bts)>maxt)][1]], tmp$date[1], units=units)) #first exit
    if (is.na(tmp1) & length(ext)>0) tmp1 <- as.numeric(difftime(tmp$date[ext[length(ext)]], tmp$date[1], units=units))  
    tmp2 <- as.numeric(difftime(tmp$date[entr[(as.numeric(bts)>maxt)][1]], tmp$date[1], units=units)) #first re-entry
    return(c(tmp1, tmp2))
  } 
  res <- data.frame(do.call(rbind, lapply(c(1:nrow(x)), fR, dframe=x, radius=radius, maxt=maxt, units=units)))
  names(res) <- c(paste("RT", radius, maxt, sep="_"), paste("T2R", radius, maxt, sep="_"))
  
  if (!addinfo) return(res)
  if (addinfo) {
    attributes(x)$infolocs <- cbind(attributes(x)$infolocs, res)
    return(x) 
  }
}

# create for-loop for calculating RT and T2R for multiple individuals
lres <- list()
for (j in 1:length(ltraj)){
  res <- ltraj[[j]][,c("x","y","date")]
  meanDist<- mean(ltraj[[j]][1:nrow(ltraj[[j]])-1,"dist"], na.rm=T)
  rads <- c(meanDist) 
  maxts <- c(12) 
  params <- expand.grid(rads=rads, maxts=maxts)
  for (i in 1:nrow(params)){
    nams <- names(res)
    tmp <- RTandT2R(ltraj[[j]], radius = params$rads[i], maxt=params$maxts[i], units="hour", addinfo = F)
    res <- cbind(res, tmp)
    names(res) <- c(nams, paste("RT", params$rads[i], params$maxts[i], sep="_"), paste("T2R", params$rads[i], params$maxts[i], sep="_"))
  }
  lres[[j]] <- res
}

# Produce a mean statistic for each individual
meanRTs <- sapply(lapply(lres, "[[", 4), function(x) mean(x, na.rm=T)); View(meanRTs)
meanT2Rs <- sapply(lapply(lres, "[[", 5), function(x) mean(x, na.rm=T)); View(meanT2Rs)


#####################################
# 4. Mean Volume of Intersection (VI)
#####################################

# For-loop to calculate overlap within an individual's 95% home range by month for multiple individuals

VI_monthly<-matrix(ncol=1, nrow=length(ltraj)) # create an empty matrix to populate with a for-loop

for (i in 1:length(ltraj)){
  ltraj[[i]] <- ltraj[[i]][complete.cases(ltraj[[i]][,c("x","y")]),] # remove NAs from coordinates
  ltraj[[i]]$month <- month(ltraj[[i]]$date)
  kudoverlap_monthly <- adehabitatHR::kerneloverlap(SpatialPointsDataFrame(ltraj[[i]][,c("x","y")], 
                                                    data=data.frame(id=ltraj[[i]]$month)), 
                                                    grid=200, method="VI")
  mw <- matrix(nr = nrow(kudoverlap_monthly), nc = nrow(kudoverlap_monthly))
  mw<-(row(mw) == col(mw) - 1) + 0 
  monthval<-kudoverlap_monthly * mw
  avg_kudoverlap_monthly<-sum(monthval)/(nrow(monthval)-1) # average month-to-month volume of intersection 
  VI_monthly[i,]<-c(avg_kudoverlap_monthly) 
}

VI_monthly<-data.frame(VI_monthly); View(VI_monthly)


############################################
# 5. Maximum Net Squared Displacement (MNSD)
############################################

# use for-loop to calculate MNSD for multiple individuals
df_ltraj<-ld(ltraj) # convert ltraj object back to data frame
id<-levels(df_ltraj$id)
maxNSD<-matrix(ncol=1, nrow=length(id)) # create empty data frame to populate with for-loop

for (i in 1:length(id)){
  NSD <- data.frame(maxNSD = df_ltraj$R2n[which(df_ltraj$id==id[i])]) # net squared displacements
  maxNSD[i,] <- sapply(NSD, function(x) max(NSD$maxNSD, na.rm=TRUE)) # calculate max NSD
}

View(maxNSD)
