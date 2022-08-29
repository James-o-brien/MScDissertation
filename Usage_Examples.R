# Usage/Examples

# Extract Environment Variables
source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
# Extract model functions and priors
source('RCode/Model.R')
# Extract the model parameterisation algorithm, default = Adaptive MCMC
source('RCode/Method.R')
# Implement the new ODD obj
#source('RCode/Val_obj.R')

# cl <- makePSOCKcluster(4)
# setDefaultCluster(cl)
# clusterExport(NULL, c('SplitSamplePop', 'fDamUnscaled','Fbdam', 'fBD', 'stochastic', 'rgammaM'))
# clusterEvalQ(NULL, library(tidyverse))
# Dam[notnans,,]<-aperm(simplify2array(parLapply(NULL, X = notnans, fun = CalcDam)), perm=c(3,2,1))
# stopCluster(cl)

#library(devtools)
#install_github('nathanvan/parallelsugar')
#library(parallelsugar)

# This file holds all the required data - hazard, exposure, vulnerability, as well as event information and observed displacement estimates in the 'ODD' class format
ODDy<-readRDS(paste0(dir,"IIDIPUS_Input/ODDobjects/EQ20210814HTI_10919_example"))
#ODDy@data$coefs <- rep(NA, length(ODDy@data$Population))
#ODDy@data$buildingsfile <- rep(NA, length(ODDy@data$Population))
#ODDy@coefs <- rep(NA, length(ODDy@cIndies))
#ODDy@buildingsfile

# This is the model parameterisation, currently trained only on earthquakes on a global level
Omega<-readRDS(paste0(dir,"IIDIPUS_Results/Omega_v2_20210828.Rdata"))

Omega <- list(Lambda1 = list(nu=0,omega=0.5),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

AlgoParams<-list(Np=20, # Number of Monte Carlo particles
                 cores=3, # Number of parallelised threads per event
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
                 epsilon_min=c(0.15,0.03,0.1), #(log) standard deviation of ABC kernels for displacement, mortality, and building destruction
                 epsilon_max=c(0.45,0.09,0.3), #initial standard deviation of ABC kernels if the kernel is shrunk over time
                 kernel='lognormal', #options are lognormal or loglaplace
                 smc_steps = 200, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 100, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
)


# Test to see if the displacement prediction calculations are working
ODDy%<>%DispX(Omega = Omega,center = Model$center,LL=F,Method = AlgoParams)
ODDpixels_agg_for_DispX1
bbox2 <- c(-73, 18.3, -72.7, 18.6)
nbuildings_HTI_2 <- ExtractOSMbuild(ODDy@bbox)

saveRDS(ODDpixels_agg, paste0(dir, 'Data_misc/ODDpixels_agg_eg.rds'))
saveRDS(ODDpolys, paste0(dir, 'Data_misc/ODDpolys_19_08.rds'))

ODDpolys <- readRDS(paste0(dir, 'Data_misc/ODDpolys_19_08.rds'))
ODDpolys@polygons
ODDpixels_agg <- readRDS(paste0(dir, 'Data_misc/ODDpixels_agg_eg.rds'))

ODDpixels_agg %<>% DispX_new(ODDpolys = ODDpolys, Omega = Omega, center = Model$center, LL=F, Method = AlgoParams)
ODDpixels_agg@predictDisp


delmoral <- readRDS(paste0(dir,'RCode/delmoral_1000part_withbdcont'))
BDX_results <- readRDS(paste0(dir,'Data_misc/BDX_results_19_08.rds'))
sum(ODDpixels_agg$Disp, na.rm = T)
bbox3 <- c(-73, 18.4, -72.8, 18.6)
nbuildings_HTI_3 <- ExtractOSMbuild(bbox3)
xy_HTI <- nbuildings_HTI_3[,c(3,4)]
nbuildings_HTI_3_spdf <- SpatialPointsDataFrame(coords = xy_HTI, data = nbuildings_HTI_3,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

ODDpixels_HTI_agg <- ParAggnBuildings(ODDy, nbuildings_HTI_3_spdf)
saveRDS(ODDpixels_HTI_agg, paste0(dir, 'Data_misc/ODDpixels_HTI_agg.rds'))

ODDpixels_HTI_agg <- readRDS(paste0(dir, 'Data_misc/ODDpixels_HTI_agg.rds'))

ODDpixels_HTI_agg%<>%DispX(Omega = Omega,center = Model$center,LL=T,Method = AlgoParams)


ODDpixels_HTI_agg@predictDisp
sum(ODDpixels_HTI_agg$nBuildings)
ODDpixels_HTI_agg




# Aggregating per polygon 19/08/22

ODDpixels_agg <- readRDS(paste0(dir, 'Data_misc/ODDpixels_agg_eg.rds'))
ODDpixels_agg %<>% DispX(ODDpolys = ODDpolys, Omega = Omega, center = Model$center, LL = F, Method = AlgoParams)
ODDpixels_agg@predictDisp

ODDpixels_agg@data$Poly <- rep(NA, length(ODDpixels_agg@data$nBuildings))
for(i in 1:length(ODDpolys@polyIndices)){
  ODDpixels_agg@data$Poly[ODDpolys@polyIndices[[i]]] <- names(ODDpolys@polygons)[i]
}
ODDpixels_agg %<>% DispX_new(Omega = Omega, center = Model$center, LL=F, Method = AlgoParams)
ODDpixels_agg@predictDisp















table(ODDpixels_agg@data$Poly)
names(ODDpolys@polygons)



ODDpixels@data$Poly <- rep(NA, length(ODDpixels@data$nBuildings))
ODDpixels 

ODDpixels[ODDpolys@polyIndices[[1]],]

for(i in 1:length(ODDpolys@polyIndices)){
  ODDpixels_agg@data$Poly <- ODDpixels[ODDpolys@polyIndices[[1]],]
}

ODDpixels_agg@data$Poly <- names(which(ODDpixels[ODDpolys@polyIndices] %in% ODDpixels@coords))



nbuildings_HTI_4 <- ExtractOSMbuild(ODDy@bbox)






nbuildings_HTI_2 <- ExtractOSMbuild(bbox)


nbuildings_HTI <- ExtractOSMbuild(ODDy@bbox)
xy_HTI <- nbuildings_HTI[,c(4,5)]
nbuildings_HTI_spdf <- SpatialPointsDataFrame(coords = xy_HTI, data = nbuildings_HTI,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))











ODDpixels_agg_for_DispX8 <- ODDpixels_agg
ODDpixels_agg_for_DispX7%<>%DispX(Omega = Omega,center = Model$center,LL=F,Method = AlgoParams)
ODDpixels_agg_for_DispX7@predictDisp

ODDpixels_agg_for_DispX8

ODDy@predictDisp


# Plot the ODD object using base R functions:
plot(ODDy) # default plots the CIESIN population data
plot(ODDy["hazMean1"]) # plots the principle hazard intensity
# Or you can use some of the ODD class methods
p<-MakeODDPlots(ODDy) # plots hazard and population side-by-side
p<-plotODDy(ODDy) # plots hazard contour lines ontop of displaced population surface plot
# Then we can also tune the plot to our greatest desires
p<-plotODDy(ODDy,breakings = c(0,10,50,100,500,1000,5000,10000),bbox=c(-74.5,17.9,-72.5,19),zoomy = 9)
p
#######################################################################


# Download and install the necessary packages, and load source files & data:
source('RCode/GetODDPackages.R')
# Search for an earthquake by providing a date range and longitude-latitude bounding box
input<-list(
  sdate=as.Date("2010-12-16"), # "YYYY-MM-DD"
  fdate=as.Date("2019-12-17"), # "YYYY-MM-DD"
  iso3="PHL", # Country code in ISO-3C form
  datadir=dir, # Location of the main folder to access the data 
  plotdir=paste0(dir,"Plots/") # Location for plots as paste0(datadir,plotdir)
)
# Or extract the data purely based on the USGS id number
input<-list(USGSid="at00qxtxcn",
            datadir=dir, # Location of the main folder to access the data 
            plotdir=paste0(dir,"Plots/") # Location for plots as paste0(datadir,plotdir)
)
#%%%%%%%%%%%%% Variables and functions from ODDRIN files %%%%%%%%%%%%%%%#
input%<>%append(list(Model=Model, # Model parameters - equations, parameters, ... (Model.R)
                     Omega=Omega, # Parameterisation of the model (GetInitialValues.R)
                     Method=AlgoParams)) # Number of CPUs, particles in SMC, etc. (Method.R)

ODDy<-AutoQuake(input)



##################################################################################
##################################################################################
###################### Running and stepping through BDX ##########################
##################################################################################
##################################################################################

BDy <- readRDS(paste0(dir,"IIDIPUS_Input/BDobjects_v3/BD_TC20200404VUT_7325"))
BDX()
BDX_results <- BDX(BD = BDy, Omega = Omega, Model = Model, Method = AlgoParams, LL = F)
BD <- BDy
setGeneric("BDX", function(BD,Omega,Model,Method,LL)
  standardGeneric("BDX") )
setMethod("BDX", "BD", function(BD,Omega,Model,Method=list(Np=20,cores=4),LL=T){
  # Only calculate buildings with all key parameters
  if(!LL) {notnans<-which(!(is.na(BD@data$Population) | is.na(BD@data$ISO3C) | is.na(BD@data$GDP)))
  } else notnans<-which(!(is.na(BD@data$Population) | is.na(BD@data$ISO3C) | is.na(BD@data$GDP) | 
                            is.na(BD@data$grading)))
  if(nrow(BD) ==0){
    if(LL){return(0)}
    else return(BD)
  }
  BD<-BD[notnans,] ;notnans<-1:nrow(BD)
  # Get parameters for model
  Params<-FormParams(BD,list(Np=Method$Np,center=Model$center))
  # Income distribution percentiles & extract income percentile  
  SincN<-seq(20,90,by = 10); Sinc<-ExtractCIndy(BD,var = paste0("p",SincN,"p100"))
  # Load buildings file
  # buildings<-readRDS(BD@buildingsfile)
  # Sample income distribution by area*building height?
  # BD%<>%SampleBuildings(buildings,F)
  hrange<-grep("hazMean",names(BD),value = T)
  # Calculate non-local linear predictor values
  LP<-GetLP(BD,Omega,Params,Sinc,notnans)
  # for each building in list,
  CalcBD<-function(ij){
    
    
    #ij <- notnans[1]
    #which(!is.na(BD@data$ISO3C))
    
    
    iso3c<-BD@data$ISO3C[ij]
    # Calculate local linear predictor (NOTE: is a scalar - we randomly sample one value)
    locallinp<-tryCatch(sample(LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]],size=Method$Np, replace=TRUE)*
                          LP$Plinp[ij]*LP$linp[[iso3c]],         error=function(e) NA) #LOOSEEND: Assumes that a house is equally likely to be from each income bracket. 
    #locallinp<- LP$dGDP$linp[5]* LP$Plinp[ij]*LP$linp[[iso3c]]
    # if(is.na(locallinp)) stop(ij)
    # locallinp<-1.
    bDamage<-0
    for(h in hrange){
      if(length(BD@data[ij,h])==0) next
      if(is.na(BD@data[ij,h])) next
      # calculate the sampled hazard intensity I_ij
      # I_ij<-rnorm(n = Method$Np,
      #             mean = BD@data[ij,h],
      #             sd = BD@data[ij,paste0("hazSD",h)]/10)
      I_ij<-BD@data[ij,h]
      # calculate the scaled damage fDamageBuilding(I,Params,Omega)
      bDamage<-bDamage+fDamageBuilding(BD@data[ij,],I_ij,
                                       Params[c("I0","Np","center")],
                                       Omega,
                                       locallinp,
                                       Params[[iso3c]]$Ik)
    }
    
    bDamage[bDamage>=1]<-1-1e-7
    bDamage[bDamage<=1e-7]<-1e-7
    # Use l_beta functions to evaluate probability of classified in group
    if(LL)  return(LL_BD(b=bDamage,classified=BD@data$grading[ij],BD_params=Model$BD_params))
    return(predBD(b=bDamage,BD_params=Model$BD_params))
  }
  if(LL) {
    if(Method$cores>1) {return(colMeans(t(matrix(unlist(mclapply(X = notnans,FUN = CalcBD,mc.cores = Method$cores)),ncol=length(notnans)))))
    } else return(colMeans(t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))))
  }
  # classified<-t(matrix(unlist(mclapply(X = notnans,FUN = predBD,mc.cores = Method$cores)),ncol=length(notnans)))
  classified<-t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))
  # Save into the file
  BD$ClassPred<-names(Model$BD_params$functions)[round(classified)]
  return(BD)
  
})

