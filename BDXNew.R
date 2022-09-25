###############
### BDX new ###
###############

library(nnet)
library(EnvStats)
library(purrr)
library(tidyverse)
library(magrittr)

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

source('RCode/Functions.R')
source('RCode/Model (1).R')

Omega <- list(Lambda1 = list(nu=0.2,omega=0.7),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=1.5,omega=3),
              theta = list(e=0.4),
              eps = list(eps=0.95))

AlgoParams<-list(Np=20, # Number of Monte Carlo particles
                 cores=1, # Number of parallelised threads per event
                 NestedCores=1, # How many cores are to be used inside the ODD displacement calculations?
                 AllParallel=F, # Do you want fully-nested (data) parallelisation?
                 itermax=2000, # How many iterations do we want?
                 ABC=-500, # Approximate Bayesian Computation rejection
                 cap=-300, # if log values are too low, then log(mean(exp(LL)))=-Inf
                 GreedyStart=0, # How sure are we of the initial covariance matrix for accepted parameters? (Larger -> more confident)
                 Pstar=0.234, # Adaptive metropolis acceptance rate
                 gamzy0=0.2, # How quickly do the rejected parameters start having an influence on the covariance? (like GreedyStart) 
                 epsilon=50, # Do we still want values at larger numbers of iterations to have an influence on the covariance?
                 minVar=1e-4, # Prevent certain parameters from being too sure of themselves
                 t_0 =200,
                 eps = 0.000000001,
                 epsilon_min=c(0.15,0.03,0.1,0.1), #(log) standard deviation of ABC kernels for displacement, mortality, and building destruction
                 epsilon_max=c(0.45,0.09,0.3,0.3), #initial standard deviation of ABC kernels if the kernel is shrunk over time
                 kernel='lognormal', #options are lognormal or loglaplace
                 smc_steps = 100, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 120, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
)

setClass("BD", 
         slots = c(hazard="character",
                   cIndies="data.frame",
                   fIndies="list",
                   I0="numeric",
                   hazdates="Date",
                   eventid="numeric",
                   coefs="numeric",
                   buildingsfile="character",
                   modifier="list"),
         contains = "SpatialPointsDataFrame")

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

setGeneric("BDX_new", function(BD,Omega,Model,Method,LL)
  standardGeneric("BDX_new") )
setMethod("BDX_new", "BD", function(BD,Omega,Model,Method=list(Np=20,cores=1),LL=F){
  # Only calculate buildings with all key parameters
  notnans <- which(!(is.na(BD@data$Population) | is.na(BD@data$grading)))
  if(nrow(BD) == 0){
    if(LL){return(0)}
    else return(BD)
  }
  
  BD<-BD[notnans,] ;notnans<-1:nrow(BD)
  Params<-FormParams(BD,list(Np=Method$Np,center=Model$center))
  
  hrange<-grep("hazMean",names(BD),value = T)
  UnscaledVals <- 0
  CalcBD<-function(ij){
    locallinp<-rep(1,10)
    for(h in hrange){
      if(length(BD@data[ij,h])==0) next
      if(is.na(BD@data[ij,h])) next
      I_ij<-BD@data[ij,h]
      UnscaledVals <- fDamUnscaled(I=I_ij, Params=Params, Omega=Omega)
    }
    return(UnscaledVals)
  }
  
  UnscaledVals <- unlist(lapply(lapply(X = notnans,FUN = CalcBD), mean))
  
  #
  ### Fitting the MLR model
  #
  multinom_model <- multinom(BD@data$grading ~ UnscaledVals + BD@data$Population, data = BD@data[notnans,])
  
  # Predicting the class
  BD@data$ClassPredicted <- predict(multinom_model, newdata = BD@data[notnans,], "class")
  
  # Building classification table
  tab <- table(BD@data$grading, BD@data$ClassPredicted)
  accuracy <- round((sum(diag(tab))/sum(tab))*100, 2)
  
  if(LL){return(-(100-accuracy))}
  
  return(BD)
  
})

BD <- readRDS(paste0(dir, "IIDIPUS_Input/BDobjects_v5/TC20200404VUT_BD"))
BDX_results <- BDX_new(BD, Omega, Model, Method = AlgoParams, LL = F)
tab <- table(BDX_results@data$grading, BDX_results@data$ClassPredicted)
