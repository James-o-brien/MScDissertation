#################
### New DispX ###
#################

library(EnvStats)
library(devtools)
library(tidyverse)
library(magrittr)
library(parallel)
library(doParallel)
library(foreach)
install_github('nathanvan/parallelsugar')
library(parallelsugar)

ExtractCentering<-function(dir, haz="TC",saver=T, input_folder='IIDIPUS_Input/'){
  
  if(saver & file.exists(paste0(dir, input_folder, "centerings"))) 
    return(readRDS(paste0(dir, input_folder, "centerings")))
  
  path<-paste0(dir, input_folder, "ODDobjects/")
  ufiles<-list.files(path=path,pattern=haz,recursive = T,ignore.case = T)
  ufiles<-ufiles[grepl(ufiles,pattern = haz)]
  GDP<-nGDP<-0
  for(fff in ufiles){
    ODDy<-readRDS(paste0(path,fff))
    
    GDP<-GDP+sum(log(ODDy@data$GDP[ODDy@data$GDP>0]),na.rm=T)
    nGDP<-nGDP+length(ODDy@data$GDP[ODDy@data$GDP>0 & !is.na(ODDy@data$GDP)])
    
  }
  
  center<-list(Pdens=log(301),dollar=GDP/nGDP)
  # center<-list(Gov=98.7,Vuln=51,CC=44,MPI=53.7,Pinf=103,Pexp=112,Sinc=0.2152956,Ik=0.4,A=3.6,H=1.65)
  print(unlist(center))
  
  saveRDS(center,paste0(dir,input_folder, "centerings"))
  
  return(center)
}

source('RCode/Model (1).R')
source('RCode/Functions.R')

setClass("ODD", 
         slots = c(dir="character",
                   hazard="character",
                   cIndies="data.frame",
                   fIndies="list",
                   IDPs="data.frame", # includes measurement dates
                   gmax="list",
                   alerts="data.frame",
                   I0="numeric",
                   hazdates="Date",
                   eventid="numeric",
                   predictDisp="data.frame",
                   modifier="list"),
         contains = "SpatialPixelsDataFrame")

setClass("ODDpolys", 
         slots = c(pixelsIndices = "list",
                   sourceInfo = "data.frame",
                   valDF = "data.frame"
         ),
         contains = "SpatialPolygonsDataFrame")

AlgoParams<-list(Np = 20, # Number of Monte Carlo particles
                 cores = 1, # Number of parallelised threads per event
                 NestedCores = 1, # How many cores are to be used inside the ODD displacement calculations?
                 AllParallel = F, # Do you want fully-nested (data) parallelisation?
                 itermax = 2000, # How many iterations do we want?
                 ABC = -500, # Approximate Bayesian Computation rejection
                 cap = -300, # if log values are too low, then log(mean(exp(LL)))=-Inf
                 GreedyStart = 0, # How sure are we of the initial covariance matrix for accepted parameters? (Larger -> more confident)
                 Pstar = 0.234, # Adaptive metropolis acceptance rate
                 gamzy0 = 0.2, # How quickly do the rejected parameters start having an influence on the covariance? (like GreedyStart) 
                 epsilon = 50, # Do we still want values at larger numbers of iterations to have an influence on the covariance?
                 minVar = 1e-4, # Prevent certain parameters from being too sure of themselves
                 t_0 = 200,
                 eps = 0.000000001,
                 epsilon_min = c(0.15,0.03,0.1,0.1), #(log) standard deviation of ABC kernels for displacement, mortality, and building destruction
                 epsilon_max = c(0.45,0.09,0.3,0.3), #initial standard deviation of ABC kernels if the kernel is shrunk over time
                 kernel = 'lognormal', #options are lognormal or loglaplace
                 smc_steps = 100, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 120, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
)

ExtractCIndy<- function(ODD,iso = NULL,var=NULL){
  cIndies<-ODD@cIndies
  if(!is.null(iso)) cIndies%<>%filter(iso3%in%iso)
  if(!is.null(var)) cIndies%<>%filter(variable%in%var)
  cIndies
}

FormParams<-function(ODD,listy){
  
  return(c(listy,list(I0=ODD@I0,fIndies=ODD@fIndies)))
  
  listy%<>%c(list(I0=ODD@I0,fIndies=ODD@fIndies))
  Ivars<-unique(ODD@cIndies$variable)
  Params<-listy
  tParams<-list()
  for (iso3c in unique(ODD@cIndies$iso3)){
    # Extract the income distribution stochastic diffusion enhancement variable
    tParams$Ik<-ExtractCIndy(ODD,iso = iso3c,var = "Ik")$value
    # Extract all the constant country specific variables
    tParams$var<-ExtractCIndy(ODD,iso = iso3c,
                              var = Ivars[!(Ivars=="Ik" | endsWith(Ivars,"p100"))])%>%
      dplyr::select(-iso3)
    # tParams$var%<>%rbind(data.frame(value=NA,variable="dollar"))
    Params[[iso3c]]<-tParams
  }
  # names(Params)[(length(listy)+1):length(Params)]<-unique(ODD@cIndies$iso3)
  return(Params)
}

LL_IDP_new<-function(Y, epsilon, kernel, cap){
  LL <- 0
  k <- 10
  impacts = c('displacement', 'mortality', 'b_destroyed', 'b_damaged') # c('gmax', 'mortality', 'buildDestroyed') #move this outside
  predictions = c('disp_predictor', 'mort_predictor', 'nBDes_predictor', 'nBDam_predictor') #c('disp_predictor', 'mort_predictor', 'nBD_predictor')
  impacts_observed <- intersect(which(impacts %in% colnames(Y)), which(predictions %in% colnames(Y)))
  
  if (kernel == 'loglaplace'){ #use a laplace kernel 
    for (i in impacts_observed){
      iso_observed <- which(!is.na(Y[,impacts[i]]) & !is.na(Y[,predictions[i]]))
      LL_impact = log(dloglap(Y[iso_observed,impacts[i]]+k, location.ald = log(Y[iso_observed,predictions[i]]+k), scale.ald = epsilon[i], tau = 0.5, log = FALSE)/
                        (1-ploglap(k, location.ald = log(Y[iso_observed,predictions[i]]+k), scale.ald = epsilon[i], tau = 0.5, log = FALSE)))
      if(all(is.finite(LL_impact))){
        LL = LL + LL_impact
      } else {
        LL = LL + cap
      }
    }
  }
  else if (kernel == 'lognormal'){ #use a lognormal kernel 
    for (i in impacts_observed){
      iso_observed <- which(!is.na(Y[,impacts[i]]) & !is.na(Y[,predictions[i]]))
      LL_impact = log(dlnormTrunc(Y[iso_observed,impacts[i]]+k, log(Y[iso_observed,predictions[i]]+k), sdlog=epsilon[i], min=k))
      if(all(is.finite(LL_impact))){
        LL = LL + LL_impact
      } else {
        LL = LL + cap
      }
    }
  }
  else {
    print(paste0("Failed to recognise kernel", AlgoParams$kernel))
    return(-Inf)
  }
  if (any(is.infinite(LL))) {
    inf_id <- which(is.infinite(LL))
    LL[inf_id] <- cap
  }
  return(LL)
}

ODDpolys <- readRDS(paste0(dir, "IIDIPUS_Input/ODDpolys"))
valDF <- ODDpolys@valDF

setGeneric("DispX_new", function(ODDpixels, valDF, Omega, center, BD_params, LL, Method, epsilon=c(0.15,0.03,0.1,0.1))
  standardGeneric("DispX_new") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX_new", "ODD", function(ODDpixels, valDF, Omega, center, BD_params, LL=F,
                                       Method=list(Np=20,cores=1,cap=-300), epsilon=c(0.15,0.03,0.1,0.1)
){

  # Extract 0D parameters & speed up loop
  Params<-FormParams(ODDpixels,list(Np=Method$Np,center=center))
  # Income distribution percentiles & extract income percentile  
  SincN<-seq(10,90,by = 10); Sinc<-ExtractCIndy(ODDpixels,var = paste0("p",SincN,"p100"))
  # Speed-up calculation (through accurate cpu-work distribution) to only values that are not NA
  if(!LL) {notnans<-which(!(is.na(ODDpixels@data$Population) | is.na(ODDpixels@data$ISO3C) | is.na(ODDpixels@data$GDP)))
  } else notnans<-which(!(is.na(ODDpixels@data$Population) | is.na(ODDpixels@data$ISO3C) | is.na(ODDpixels@data$GDP) |
                            !ODDpixels@data$ISO3C%in%ODDpixels@gmax$iso3))
  
  BD_data_present <- ifelse(all(is.na(ODDpixels@data$nBuildings)) , F, T)
  # Calculate non-local linear predictor values
  #LP<-GetLP(ODDpixels,Omega,Params,Sinc,notnans)
  LP<-rep(1,10) # Setting linear predictor term to 1
  # Speed things up a little
  hrange<-grep("hazMean",names(ODDpixels),value = T)
  # Function to predict displacement per gridpoint
  CalcDam<-function(ij){
    iso3c<-ODDpixels@data$ISO3C[ij]
    # Calculate local linear predictor (NOTE: is a vector due to income distribution)
    # locallinp<-LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] #LOOSEEND
    locallinp<-rep(1,10) # Setting local linear predictor term to 1
    
    # Sample population per income distribution (Assumes 9 percentiles):
    lPopS <- SplitSamplePop(Pop=ODDpixels@data$Population[ij],Method$Np) 
    tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
    tPop[3,]=colSums(lPopS)
    for(h in hrange){
      if(is.na(ODDpixels@data[ij,h])) next
      # Resample population based on who is remaining
      ind<-tPop[3,]>0
      if(h!=hrange[1]) {
        if(sum(ind)==0) break #if no remaining population, skip modelling
        if(length(lPopS[,ind])==0) break #if no remaining population, skip modelling
        lPopS[,ind]<-SplitSamplePop(Pop=tPop[3,ind])
      }
      I_ij<-ODDpixels@data[ij,h]
      
      # Separate into income distributions (as each have 10% of population, order doesn't matter)
      for (s in 1:length(SincN)){
        if(all(lPopS[s,]==0)) next
        # Predict damage at coordinate {i,j} (vector with MC particles)
        Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[s], error=function(e) NA)
        if(any(is.na(Damage))) print(ij)
        
        D_MortDisp <- D_MortDisp_calc(Damage, Omega) #First row of D_MortDisp is D_Mort, second row is D_Disp
        
        # Accumulate the number of people displaced/deceased, but don't accumulate the remaining population
        tPop[3,ind]<-0
        
        tPop[,ind]<-tPop[,ind] + Fbdam(lPopS[s,ind],D_MortDisp[2,ind], D_MortDisp[1,ind], (1-D_MortDisp[1,ind]-D_MortDisp[2,ind]))
      } 
    }
    #ensure the total displaced, deceased or remaining does not exceed total population
    tPop[tPop>ODDpixels@data$Population[ij]]<-floor(ODDpixels@data$Population[ij])
    
    #if no building damage/destruction data:
    if(!BD_data_present) return(rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np), rep(NA, Method$Np))) #return total displacement and mortality, set number of buildings damaged and destroyed to NA
    
    #otherwise, sample from the model for the number of buildings destroyed:
    #we take locallinp[5] which corresponds to locallinp for the median GDP
    nBuildings <- rep(ODDpixels@data$nBuildings[ij], Method$Np)
    nBDes <- rep(0, Method$Np)
    nBDam <- rep(0, Method$Np)
    Builds <- rbind(rep(0, Method$Np), rep(0, Method$Np))
    for (h in hrange){
      
      if(is.na(ODDpixels@data[ij,h])) next
      if(is.na(ODDpixels@data[ij,h])) next
      if(all(nBuildings==0)) break #if no remaining buildings, skip modelling LOOSEEND
        
      I_ij <- ODDpixels@data[ij,h]
      Damage <- tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[5], error=function(e) NA) #calculate unscaled damage (excluding GDP)
      D_extent <- plnorm(Damage, meanlog = Omega$Lambda3$nu, sdlog = Omega$Lambda3$omega)
        
      for(i in 1:length(nBuildings)){
        if(D_extent[i] > 0.5){
          nBDes[i] <- nBuildings[i]
          } else if(D_extent[i] >= 0.2 & D_extent[i] <= 0.5){
            nBDam[i] <- nBuildings[i]
          } else next 
        }
    
        moreBD = rbind(nBDes, nBDam)
        Builds = Builds + moreBD
      
      }
    
    return(rbind(tPop[1:2,,drop=FALSE], Builds))
  }
  
  Dam<-array(0,c(nrow(ODDpixels),Method$Np,4)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Destroyed, 
  # Dam[,,4] = Buildings Damaged
  
  if(Method$cores>1) { Dam[notnans,,]<-aperm(simplify2array(mclapply(X = notnans,FUN = CalcDam,mc.cores = Method$cores)), perm=c(3,2,1))
  } else  Dam[notnans,,]<- aperm(simplify2array(lapply(X = notnans,FUN = CalcDam)), perm=c(3,2,1))
  
  funcy <- function(i, LLout = T, epsilon = AlgoParams$epsilon_min, kernel = AlgoParams$kernel, cap = AlgoParams$cap) {
    tmp <- data.frame(polygon = ODDpixels@data$polygon, IDPs = Dam[,i,1], mort = Dam[,i,2], nBDes = Dam[,i,3]/Method$Np, nBDam = Dam[,i,4]/Method$Np) %>% # average over the MC particles
      group_by(polygon) %>% summarise(disp_predictor = floor(sum(IDPs, na.rm = T)), 
                                   mort_predictor = floor(sum(mort, na.rm = T)),
                                   nBDes_predictor = floor(sum(nBDes, na.rm = T)),
                                   nBDam_predictor = floor(sum(nBDam, na.rm = T)),
                                   .groups = 'drop_last')
    tmp <- tmp[!is.na(tmp$polygon), ]
    tmp %<>% merge(valDF, by = "polygon") %>% arrange(desc(disp_predictor))
    tmp$nBUnaf_predictor <- tmp$b_total - tmp$nBDam_predictor - tmp$nBDes_predictor
    
    if(LLout) {
      return(LL_IDP_new(tmp, epsilon,  kernel, cap))
    }
    return(tmp)
  }
  
  if (LL == F & Method$Np == 1){
    ODDpixels@predictDisp <- funcy(1,LLout = F) 
    
    return(ODDpixels)
  }
  outer <- sapply(1:Method$Np, funcy, epsilon = epsilon) 
  outer[outer < (-745)]<- -745
  
  # Find the best fit solution
  if(length(unique(ODDpixels@data$polygon)) > 1) {
    if(LL)  return(log(rowMeans(exp(outer), na.rm = T)))
    MLE <- which.max(log(colSums(exp(outer), na.rm = T)))
  }  else {
    return(outer) #SMC-CHANGE
    if(LL)  return(log(mean(exp(outer), na.rm = T)))
    MLE <- which.max(log(exp(outer)))
  }
  if(Method$Np == 1){
    MLE = 1
  }
  # Save into ODDpixels object
  ODDpixels@data$Disp <- Dam[,MLE,1]  
  ODDpixels@data$Mort <- Dam[,MLE,2]
  ODDpixels@data$nBDes <- Dam[,MLE,3]
  ODDpixels@data$nBDam <- Dam[,MLE,4]
  ODDpixels@predictDisp <- funcy(MLE, LLout = F) 
  
  return(ODDpixels)
  
})

ODDpixels <- readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects_v4/TC20200404VUT_2208_agg_4"))
Omega <- list(Lambda1 = list(nu=0.2,omega=0.7),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=1.5,omega=3),
              theta = list(e=0.4),
              eps = list(eps=0.95))
ODDpixels %<>% DispX_new(valDF = valDF, Omega = Omega, center = Model$center, LL=F, Method = AlgoParams)
ODDpixels@predictDisp
