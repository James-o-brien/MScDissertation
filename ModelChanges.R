#####################
### Model changes ###
#####################

dir<-directory<-"C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/"
setwd(directory)

source('RCode/Model.R')
source('RCode/Functions.R')
source('RCode/DispXNew.R')
source('RCode/BDXNew.R')

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

Omega <- list(Lambda1 = list(nu=0.2,omega=0.7),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=1.5,omega=3),
              theta = list(e=0.4),
              eps = list(eps=0.95))

#
### LL_Displacement_new
#

# Log-likelihood for displacement (ODD) objects
LL_Displacement_new <- function(LL,valDF,dir,Model,proposed,AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1,0.1)){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects_v4/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task parallelism: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many ODD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))]) #looseend
    
    tmpFn<-function(filer){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,filer))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
  #    ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX_new(ODDpixels = ODDy, valDF = valDF,Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams, epsilon=epsilon),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",filer));return(-Inf)}
      
      return(tLL) #SMC-CHANGE
      # Weight the likelihoods based on the number of events for that country
      cWeight<-1
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL))
      else return(cWeight*mean(tLL,na.rm=T))
    }
    return(colSums(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores)))) # SMC-CHANGE
    return(sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
    
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    for(i in 1:length(ufiles)){
      # Extract the ODD object
      ODDy<-readRDS(paste0(folderin,ufiles[1]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      #ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX_new(ODDpixels = ODDy, valDF = valDF, Omega = proposed,center = Model$center, LL = T,Method = AlgoParams, epsilon=epsilon),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate Disp LL of ",ufiles[i]));return(-Inf)}
      # Weight the likelihoods based on the number of events for that country
      cWeight<-1
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Add the likelihood to the list of all events.
      if(expLL) {LL<-LL+cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)
      } else LL<-LL+cWeight*mean(tLL,na.rm=T)
    }
    
    return(LL)
    
  }
  
}

#
### LL_Buildings_new
#

# Log-likelihood for building damage (BD) objects
LL_Buildings_new<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
  # Load BD files
  folderin<-paste0(dir,"IIDIPUS_Input/BDobjects_v5/")
  ufiles<-list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)
  
  # Parallelise appropriately
  if(AlgoParams$AllParallel){
    # Task Parallelisation: this parallelisation calculates each event side-by-side, which is ideal if we have many CPU threads available and many BD objects
    cores<-AlgoParams$cores
    AlgoParams$cores<-AlgoParams$NestedCores
    # When using task parallelisation, put the heaviest files first for optimisation reasons
    x <- file.info(paste0(folderin,ufiles))
    ufiles<-na.omit(ufiles[match(length(ufiles):1,rank(x$size))])
    
    tmpFn<-function(filer){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,filer))
      if(nrow(BDy@data)==0){return(0)}
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
                    error=function(e) NA)
      
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",filer));return(-Inf)}
      return(tLL) #SMC-CHANGE
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Return the average log-likelihood
      cWeight<-1 #looseend
      if(is.na(cWeight)){return(-Inf)}                               
      if(expLL) return(cWeight*(log(mean(exp(tLL-maxLL),na.rm=T))+maxLL)) 
      else return(cWeight*mean(tLL,na.rm=T))
      
    }
    return(LL + colSums(do.call(rbind, mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores)))) #SMC-CHANGE
    return(LL + sum(unlist(mclapply(X = ufiles,FUN = tmpFn,mc.cores = cores))))
    
  } else {
    
    # Data parallelism: this is nested parallelisation, ideal if we have low CPU threads and large but few ODD files
    for(i in 1:length(ufiles)){
      # Extract the BD object
      BDy<-readRDS(paste0(folderin,ufiles[1]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
 #     BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX_new(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
                    error=function(e) NA)
      # If all is good, add the LL to the total LL
      if(any(is.infinite(tLL)) | all(is.na(tLL))) {print(paste0("Failed to calculate BD LL of ",ufiles[i]));return(-Inf)}
      # We need the max to ensure that exp(Likelihood)!=0 as Likelihood can be very small
      maxLL<-max(tLL,na.rm = T)
      # Add the likelihood to the list of all events.
      if(expLL) LL<-LL+log(mean(exp(tLL-maxLL),na.rm=T))+maxLL
      else LL<-LL+mean(tLL,na.rm=T)
    }
    
    return(LL)
  }
}

#
### logTarget_new
#

# Bayesian Posterior distribution
logTarget_new <- function(dir,Model,proposed,AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1,0.1)){
  
  # Apply higher order priors
  #  if(!is.null(Model$HighLevelPriors)){
  #HP<-Model$HighLevelPriors(proposed,Model)
  #    print(paste0("Higher Level Priors = ",HP))
  #  }
  
  # Approximate Bayesian Computation rejection
  #  if(HP>AlgoParams$ABC) return(-5000000) # Cannot be infinite, but large (& negative) to not crash the Adaptive MCMC algorithm
  #  if (HP> AlgoParams$ABC) return(-Inf)
  
  # Add the log-likelihood values from the ODD (displacement) objects
  LL<-LL_Displacement_new(0,valDF,dir,Model,proposed,AlgoParams,expLL=T, epsilon)
  if (any(LL == 0)){ #SMC-CHANGE
    return(Inf)
  }
  # Add the log-likelihood values from the BD (building damage) objects
  LL%<>%LL_Buildings_new(dir,Model,proposed,AlgoParams,expLL=T)
  
  posterior<-LL 
  # Add Bayesian priors
  if(!is.null(Model$priors)){
    posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
  }

  return(posterior)
  
}

#
### Uncomment to run an example below: WARNING (this make take a couple of minutes) requires downloading the data and storing in the appropriately named folders
#

#logTarget_new(dir,Model,proposed=Omega,AlgoParams,expLL=T)
