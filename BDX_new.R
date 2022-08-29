Omega <- list(Lambda1 = list(nu=0,omega=0.5),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

AlgoParams<-list(Np=1, # Number of Monte Carlo particles
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
                 epsilon_min=c(0.15,0.03,0.1), #(log) standard deviation of ABC kernels for displacement, mortality, and building destruction
                 epsilon_max=c(0.45,0.09,0.3), #initial standard deviation of ABC kernels if the kernel is shrunk over time
                 kernel='lognormal', #options are lognormal or loglaplace
                 smc_steps = 200, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 100, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
)

ODDpoints <- readRDS(paste0(dir,'Data_misc/TC20200404VUT_7325_points'))

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


ODDpoints <- readRDS(paste0(dir,'Data_misc/TC20200404VUT_7325_points'))
ODDpoints@data$Poly <- rep(NA, length(ODDpoints@data$grading))
for(i in 1:length(ODDpolys@pointsIndices)){
  ODDpoints@data$Poly[ODDpolys@pointsIndices[[i]]] <- names(ODDpolys@polygons)[i]
}

BDX_results_new <- BDX_new(BD = ODDpoints, Omega = Omega, Model = Model, Method = AlgoParams, LL = F)
table(BDX_results_new$ClassPred)


setGeneric("BDX_new", function(BD,Omega,Model,Method,LL)
  standardGeneric("BDX_new") )
setMethod("BDX_new", "BD", function(BD,Omega,Model,Method=list(Np=1,cores=1),LL=F){
  # Only calculate buildings with all key parameters
  if(!LL) {notnans<-which(!(is.na(BD@data$Population) | is.na(BD@data$ISO3C) | is.na(BD@data$GDP)))
  } else notnans<-which(!(is.na(BD@data$Population) | is.na(BD@data$ISO3C) | is.na(BD@data$GDP) | 
                            is.na(BD@data$grading)))
  if(nrow(BD) == 0){
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
    if(Method$cores>1) {
      LL_impact <- 0
      for(i in 1:length(unique(ODDpoints@data$Poly[notnans]))){
        polygon <- unique(ODDpoints@data$Poly[notnans])[i]
        inds <- which(ODDpoints[notnans,]@data$Poly == polygon)
        LL_impact[i] <- colMeans(t(matrix(unlist(mclapply(X = inds, FUN = CalcBD)), ncol=length(inds))))
      }
    } else 
      LL_impact <- 0
      for(i in 1:length(unique(ODDpoints@data$Poly[notnans]))){
        polygon <- unique(ODDpoints@data$Poly[notnans])[i]
        inds <- which(ODDpoints[notnans,]@data$Poly == polygon)
        LL_impact[i] <- colMeans(t(matrix(unlist(lapply(X = inds, FUN = CalcBD)), ncol=length(inds))))
      }
    return(LL_impact)
  }
  
  # classified<-t(matrix(unlist(mclapply(X = notnans,FUN = predBD,mc.cores = Method$cores)),ncol=length(notnans)))
  classified<-t(matrix(unlist(lapply(X = notnans,FUN = CalcBD)),ncol=length(notnans)))
  # Save into the file
  BD$ClassPred<-names(Model$BD_params$functions)[round(classified)]
  
  return(BD)
  
})

if (length(unique(ODDpoints@data$Poly[notnans])) != length(unique(ODDpixels@data$Poly[notnans]))){ # these are different notnans!!!!!!!!!!1!!!!
  # do something that fixes this!
  # might be adding a 0 in the position of the vector that is missing in ODDpoints. Will probably have to order the rows in both @data slots
  # to be alphabetically sorted according to the Polys column in both, this should make the code output the vector in the same order for both
}