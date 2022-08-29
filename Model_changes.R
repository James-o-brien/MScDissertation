# Model changes

AlgoParams<-list(Np = 20, # Number of Monte Carlo particles
                 cores = 2, # Number of parallelised threads per event
                 NestedCores = 2, # How many cores are to be used inside the ODD displacement calculations?
                 AllParallel = T, # Do you want fully-nested (data) parallelisation?
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
                 epsilon_min = c(0.15,0.03,0.1), #(log) standard deviation of ABC kernels for displacement, mortality, and building destruction
                 epsilon_max = c(0.45,0.09,0.3), #initial standard deviation of ABC kernels if the kernel is shrunk over time
                 kernel = 'lognormal', #options are lognormal or loglaplace
                 smc_steps = 2, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 2, #Number of particles in the ABC-SMC algorithm
                 smc_alpha = 0.9
)

Omega <- list(Lambda1 = list(kappa=5,nu=10,omega=10),
              Lambda2 = list(nu= 1.3, omega=0.9),
#              Lambda3 = list(nu=1.5,omega=3),
              zeta = list(k=1.5, lambda=1),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.9),
              eps = list(eps=0.95))

# Omega <- list(Lambda1 = list(nu=0,omega=0.5),
#               Lambda2 = list(nu= 1.3, omega=0.9),
# #              Lambda3 = list(nu=0.4,omega=0.6),
#               zeta = list(k=2.978697, lambda=1.405539),
#               Pdens = list(M=0.02988616, k = 6.473428),
#               dollar = list(M = -1.051271, k = 6.473428),
#               theta = list(e=0.2359788),
#               eps = list(eps=0.01304351))


#
### LL_Displacement
#
ODDy <- ODDpixels_2208_agg_4
ODDpixels_2208_agg_4 <- readRDS(paste0(dir, "IIDIPUS_Input/ODDobjects_v4/TC20200404VUT_2208_agg_4")) # this has nBuildings and Poly

ODDpixels_2208_agg_4 %<>% DispX_new(ODDpolys = ODDpolys, Omega = Omega,center = Model$center,LL=F,Method = AlgoParams)
ODDpixels_2208_agg_4
ODDpixels_2208_agg_4@predictDisp

Omega <- list(Lambda1 = list(kappa=5,nu=10,omega=10),
              Lambda2 = list(nu= 1.3, omega=0.9),
              #              Lambda3 = list(nu=1.5,omega=3),
              zeta = list(k=1.5, lambda=1),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.9),
              eps = list(eps=0.95))


# new structure
#LL_Displacement_new(0,ODDpolys=ODDpolys,dir=dir,Model=Model,proposed=Omega,AlgoParams=AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1))
ODDy <- ODDpixels
proposed <- Omega

saveRDS(ODDpolys, file=paste0(dir,"IIDIPUS_Input/ODDpolys"))
rm("ODDpolys")
gc()
LL_Displacement_new(0,dir=dir,Model=Model,valDF=valDF,proposed=Omega,AlgoParams=AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1,0.1,0.1))

# --------------------------------- New versions of functions in Model.R ------------------------------ # 

# Log-likelihood for displacement (ODD) objects
LL_Displacement_new <- function(LL,dir,Model,valDF,proposed,AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1,0.1,0.1)){
  
  # Load ODD files
  folderin<-paste0(dir,"IIDIPUS_Input/ODDobjects_v4/")
  ufiles<-na.omit(list.files(path=folderin,pattern=Model$haz,recursive = T,ignore.case = T)) #looseend
  #ODDpolys <- readRDS(paste0(dir,"IIDIPUS_Input/ODDpolys"))

    
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
#      ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX_new(ODDpixels = ODDy, valDF = valDF, Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams, epsilon=epsilon),
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
      ODDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
#      ODDy@fIndies<-Model$fIndies
      ODDy@gmax%<>%as.data.frame.list()
      # Apply DispX
      tLL<-tryCatch(DispX_new(ODDpixels = ODDy, valDF=valDF, Omega = proposed,center = Model$center, BD_params = Model$BD_params, LL = T,Method = AlgoParams, epsilon=epsilon),
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
### LL_Buildings
#

BDX_results <- LL_Buildings(0, dir, Model, proposed=Omega, AlgoParams, expLL=T)

# Log-likelihood for building damage (BD) objects
LL_Buildings<-function(LL,dir,Model,proposed,AlgoParams,expLL=T){
  # Load BD files
  folderin<-paste0(dir,"IIDIPUS_Input/BDobjects_v3/")
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
      #BDy@fIndies<-Model$fIndies
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
      BDy<-readRDS(paste0(folderin,ufiles[i]))
      # Backdated version control: old IIDIPUS depended on ODDy$fIndies values and gmax different format
      BDy@fIndies<-Model$fIndies
      # Apply BDX
      tLL<-tryCatch(BDX(BD = BDy,Omega = proposed,Model = Model,Method=AlgoParams, LL=T),
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
### logTarget
#


Omega <- list(Lambda1 = list(nu=10,omega=10),
              Lambda2 = list(nu= 1.3, omega=0.9),
              #              Lambda3 = list(nu=1.5,omega=3),
              zeta = list(k=1.5, lambda=1),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.9),
              eps = list(eps=0.95))

Omega <- list(Lambda1 = list(kappa = 10, nu=10,omega=0.5),
              Lambda2 = list(nu= 2.3, omega=2.9),
#              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=0.978697, lambda=0.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.0359788),
              eps = list(eps=0.01304351))


HLPrior_sample <- function(Model, AlgoParams){
  #Draws a sample from the prior that accepts the Higher Level Prior
  #Returns sample on the transformed space (not the physical space)
  n_x <- length(unlist(Model$links))
  HP <- AlgoParams$ABC + 1
  sample <- rep(0, n_x)
  while (HP > AlgoParams$ABC){
    sample <- runif(n_x, Model$par_lb, Model$par_ub) #generate proposal on the physical space
    HP <- Model$HighLevelPriors(relist(sample,skeleton=Model$skeleton),Model) #check higher level prior of the proposal
  }
  return(sample %>% relist(Model$skeleton)%>% Physical2Proposed(Model) %>% unlist())
}


Model$HighLevelPriors<-function(Omega,Model,modifier=NULL){
  
  if(!is.null(modifier)) lp<-exp(as.numeric(unlist(modifier))) else lp<-1.
  lp <- c(0, 1, 10) # lp_range - 0.361022 corresponds to lp for minimum GDP and PDens scaling, 2.861055 corresponds to maximum
  if(Model$haz=="EQ"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 4.5,theta = Omega$theta) 
    Damfun<-function(I_ij, type) {
      D <- Dfun(I_ij)
      if (type=='Damage'){
        return(D)
      }
      if (type=='Buildings Destroyed'){
        return(plnorm(D,Omega$Lambda3$nu, Omega$Lambda3$omega))
      }
      DispMort <- plnorm(D,Omega$Lambda1$nu, Omega$Lambda1$omega)
      Mort <- plnorm(D,Omega$Lambda2$nu, Omega$Lambda2$omega)
      Disp <- ifelse(DispMort - Mort < 0, 0, DispMort - Mort)
      if (type=='Displacement'){
        return(Disp)
      }
      else if (type=='Mortality'){
        return(Mort)
      }
      else if (type == 'DispMort'){
        return(DispMort)
      }
    }
    adder<-0
    # Add lower bound priors:
    adder<-sum(pweibull(Damfun(4.6, type='Displacement')*lp[1],3,0.01), pweibull(Damfun(4.6, type='Mortality')*lp[1],3,0.01),
               pweibull(Damfun(4.6, type='Buildings Destroyed')*lp[1],3,0.01), pweibull(Damfun(4.6, type='Damage')*lp[1],3,0.01))
    # Middle range priors:
    adder<-adder+sum(pweibull(Damfun(6, type='Displacement')*lp[2],15,0.8), 1- pweibull(Damfun(6, type='Displacement')*lp[2],3,0.005),
                     pweibull(Damfun(6, type='Mortality')*lp[2],15,0.1),# 1- pweibull(Damfun(6, type='Mortality')*lp[2],3,0.005),
                     pweibull(Damfun(6, type='Buildings Destroyed')*lp[2],15,0.8), 1- pweibull(Damfun(6, type='Buildings Destroyed')*lp[2],3,0.005),
                     pweibull(Damfun(6, type='Damage')*lp[2],15,0.8), 1- pweibull(Damfun(6, type='Damage')*lp[2],3,0.005),
                     pweibull(Damfun(6, type='DispMort')*lp[2], 15, 0.8)
    )
    # Upper bound priors:
    adder <- adder + pweibull(Damfun(7, type='Mortality')*lp[2],15,0.2)
    adder<-adder+sum(1-pweibull(Damfun(9, type='Displacement')*lp[3],5,0.3), 1-pweibull(Damfun(9, type='Mortality')*lp[3],5,0.2),
                     1-pweibull(Damfun(9, type='Buildings Destroyed')*lp[3],5,0.5), 1-pweibull(Damfun(9, type='Damage')*lp[3],30,0.85))
    
    # adder<-rep(0,length(lp))
    # for(i in 1:length(adder)){
    #   # Add lower bound priors:
    #   adder[i]<-sum(pweibull(Damfun(4.6, type='Displacement')*lp[i],3,0.01), pweibull(Damfun(4.6, type='Mortality')*lp[i],3,0.01),
    #                 pweibull(Damfun(4.6, type='Buildings Destroyed')*lp[i],3,0.01), pweibull(Damfun(4.6, type='Damage')*lp[i],3,0.01))
    #   # Middle range priors:
    #   adder[i]<-adder[i]+sum(pweibull(Damfun(6, type='Displacement')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Displacement')*lp[i],3,0.005),
    #                          pweibull(Damfun(6, type='Mortality')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Mortality')*lp[i],3,0.005),
    #                          pweibull(Damfun(6, type='Buildings Destroyed')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Buildings Destroyed')*lp[i],3,0.005),
    #                          pweibull(Damfun(6, type='Damage')*lp[i],15,0.8), 1- pweibull(Damfun(6, type='Damage')*lp[i],3,0.005)
    #                          )
    #   # Upper bound priors:
    #   adder[i]<-adder[i]+sum(1-pweibull(Damfun(9, type='Displacement')*lp[i],5,0.3), 1-pweibull(Damfun(9, type='Mortality')*lp[i],5,0.2),
    #                          1-pweibull(Damfun(9, type='Buildings Destroyed')*lp[i],5,0.5), 1-pweibull(Damfun(9, type='Damage')*lp[i],30,0.85))
    # }
    return(sum(adder)) #looseend: need to address when including modifiers
    #return(adder)
    
  } else if(Model$haz=="TC"){
    
    Dfun<-function(I_ij) h_0(I = I_ij,I0 = 3,theta = Omega$theta) 
    Dispfun<-function(I_ij) c(BinR(Dfun(I_ij)*Dfun(I_ij)*Omega$Lambda1$kappa+Omega$Lambda1$nu*Dfun(I_ij) + Omega$Lambda1$omega,Omega$zeta)%o%lp)
    Damfun<-function(I_ij) c(BinR(Dfun(I_ij),Omega$zeta)%o%lp)
    
    # Add lower bound priors:
    adder<-sum(50*pweibull(c(Dispfun(3.05),Damfun(3.05)),3,0.001))
    # Middle range priors:
    # adder<-adder+sum(15*(1- pweibull(c(Dispfun(3.53),Damfun(3.53)),3,0.005)) +15*pweibull(c(Dispfun(3.53),Damfun(3.53)),15,0.8) )
    # Upper bound priors:
    adder<-adder+sum(50*(1-pweibull(c(Dispfun(45),Damfun(45)),30,0.85)))
    
    return(-adder)
    
  }
  # I<-seq(from=4.05,to=9.5,length.out = 200)
  # Intensity<-data.frame(I_ij=rep(I,2),value=c(vapply(I,Dispfun,numeric(1)),vapply(I,BDfun,numeric(1))),
  #                       term=c(rep("Displacement",200),rep("Building Damage",200)))
  
}

logTarget_new(dir,Model,proposed,AlgoParams,expLL=T,epsilon=c(0.15,0.03,0.1,0.1,0.1))

# Bayesian Posterior distribution: This is for a group of ODD objects with observed data
logTarget_new <- function(dir,Model,proposed,AlgoParams,expLL=T, epsilon=c(0.15,0.03,0.1,0.1,0.1)){
  
#  Apply higher order priors
#  if(!is.null(Model$HighLevelPriors)){
##    HP<-Model$HighLevelPriors(proposed,Model)
#    print(paste0("Higher Level Priors = ",HP))
#  }
  
  # Approximate Bayesian Computation rejection
#  if(HP>AlgoParams$ABC) return(-5000000) # Cannot be infinite, but large (& negative) to not crash the Adaptive MCMC algorithm
#  if (HP> AlgoParams$ABC) return(-Inf)
  
  # Add the log-likelihood values from the ODD (displacement) objects
  LL<-LL_Displacement_new(0,dir,Model,valDF,proposed,AlgoParams,expLL=T, epsilon)
  #print(paste0("LL Displacements = ",LL)) ; sLL<-LL
  if (any(LL == 0)){ #SMC-CHANGE
    return(Inf)
  }
  # Add the log-likelihood values from the BD (building damage) objects
  LL%<>%LL_Buildings(dir,Model,proposed,AlgoParams,expLL=T)
  #print(paste0("LL Building Damages = ",LL-sLL))
  
  posterior<-LL #+HP
  # Add Bayesian priors
  if(!is.null(Model$priors)){
    posterior<-posterior+sum(Priors(proposed,Model$priors),na.rm=T)
  }
  #print(paste0("Posterior = ",posterior))
  
  return(posterior)
  
}


# #Omega <- list(Lambda1 = list(nu=10,omega=10),
# #              Lambda2 = list(nu= 1.3, omega=0.9),
# #              Lambda3 = list(nu=1.5,omega=3),
#               zeta = list(k=5, lambda=5),
#               Pdens = list(M=0.02988616, k = 6.473428),
#               dollar = list(M = -1.051271, k = 6.473428),
#               theta = list(e=0.9),
#               eps = list(eps=0.95))



Proposed2Physical<-function(proposed,Model,index=NULL){
  
  Model$links%<>%unlist()
  
  if(is.null(index)) index<-1:length(Model$links)
  # Link functions to convert values into useable/physical values
  for (i in index)  {
    proposed[i] <- match.fun(Model$links[[names(proposed)[i]]])(proposed[i], Model$par_lb[i],Model$par_ub[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
}

Physical2Proposed<-function(proposed,Model,index=NULL){
  
  proposed%<>%unlist()
  Model$unlinks%<>%unlist()
  
  if(is.null(index)) index<-1:length(Model$unlinks)
  # Link functions to convert values into useable/physical values
  for (i in index)  {
    proposed[i] <- match.fun(Model$unlinks[[names(proposed)[i]]])(proposed[i], Model$par_lb[i],Model$par_ub[i])
  }
  # Reshape into desired structure
  proposed%>%relist(skeleton=Model$skeleton)
  
}

abcSmc_delmoral_new(AlgoParams,Model)
abcSmc_delmoral_new <- function(AlgoParams, Model, unfinished=F, oldtag=''){
  #Input: 
  # - AlgoParams: Parameters describing the ABC-SMC Algorithm (e.g. the ABC rejection threshold for higher level priors)
  # - Model: Describes the data simulation and calculation of the distance measure 
  # - Unfinished: If TRUE, then include oldtag - the tag (end of the filename) of an unfinished ABC-SMC run to be completed.
  # Output:
  # - A list containing the parameter values for each particle at each step of the ABC-SMC algorithm, as well as the
  #   corresponding weights, distances and final perturbation kernel. 
  # Details:
  # - Uses the algorithm described in Del Moral et al., 2011: https://link.springer.com/content/pdf/10.1007/s11222-011-9271-y.pdf
  # - Uses the multivariate perturbation kernel described in Filippi et al., 2012: https://arxiv.org/pdf/1106.6280.pdf
  
  steps = AlgoParams$smc_steps
  Npart = AlgoParams$smc_Npart
  M = AlgoParams$Np
  
  start_time <- Sys.time()
  n_x <- length(Model$par_lb) #n_x = number of parameters
  Omega_sample <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the transformed space
  Omega_sample_phys <- array(NA, dim=c(Npart, n_x, steps)) #store sampled parameters on the untransformed space
  W <- array(NA, dim=c(Npart, steps)) #Weights
  d <- array(Inf, dim=c(Npart, M, steps)) #Distances
  npa=rep(0,Npart) #Number of alive particles (have at least one distance less than the tolerance)
  tolerance = 10000000 #Ensure all initial distances are less than the tolerance
  tolerancetarget=1 #LOOSEEND: need to add in an adaptive stopping rule
  tolerancestore=c(tolerance)
  essstore=c(Npart)
  N_T=Npart/2
  alpha=AlgoParams$smc_alpha
  
  tpa<-function(tolerance,x){ #calculate the total proportion of alive particles (have at least one distance less than the tolerance)
    d<-abs(x)
    for (i in 1:Npart){
      count<-0
      for (j in 1:M){
        if (d[i,j]<tolerance){
          count<-count+1
        }
      }
      npa[i]<-count
    }
    return(sum(npa>0)/Npart)
  }
  
  tag<-gsub(gsub(Sys.time(),pattern = " ", replacement = "_"),pattern = ":",replacement = "")
  if(unfinished==F){ #Initialize and perform sampling for s=1
    W[,1] <- 1/Npart # Give each particle an equal weight
    for (n in 1:Npart){
      print(paste('Step: 1, Particle:', n))
      #draw initial particle values from the prior and ensure that they satisfy the higher level prior
      Omega_sample[n,,1] <- HLPrior_sample(Model, AlgoParams)
      Omega_sample_phys[n,,1] <-  Omega_sample[n,,1] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
      #calculate distance
      d[n,,1] <- -logTarget_new(dir = dir,Model = Model,
                            proposed = Omega_sample_phys[n,,1] %>% relist(skeleton=Model$skeleton), 
                            AlgoParams = AlgoParams)
    }
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=array(0, dim=c(n_x,n_x)),
           W=W, 
           tolerance = tolerance),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
    s_start = 2 # continue the algorithm from s = 2
    n_start = 1
  } else { #Collect relevant information from the unfinished sample
    output_unfinished <- readRDS(paste0(dir,"IIDIPUS_Results/abcsmc_",oldtag))
    s_finish = max(which(colSums(!is.na(output_unfinished$omegastore[,1,]))>0)) #find last completed step
    n_finish = max(which(!is.na(output_unfinished$omegastore[,1,s_finish])))
    W[,1:s_finish] <- output_unfinished$W[,1:s_finish]
    d[,,1:s_finish] <- output_unfinished$d[,,1:s_finish]
    propCOV <- output_unfinished$propCOV
    Omega_sample_phys[,,1:s_finish] <- output_unfinished$omegastore[,,1:s_finish]
    for (n in 1:Npart){
      for (s in 1:s_finish){
        #convert particles from the physical to the transformed space
        if(!is.na(Omega_sample_phys[n,1,s]))
          Omega_sample[n,,s] <- Omega_sample_phys[n,,s] %>% relist(skeleton=Model$skeleton) %>% Physical2Proposed(Model) %>% unlist()
      }
    }
    s_start = ifelse(n_finish==Npart, s_finish+1, s_finish) #identify the appropriate step from which to continue the algorithm
    n_start = ifelse(n_finish==Npart, 1, n_finish + 1) #identify the appropriate particle from which to continue the algorithm
    tolerance = output_unfinished$tolerance
  }
  
  for (s in s_start:steps){
    
    #record and print time for each step
    end_time <- Sys.time()
    print(paste('Time:', end_time-start_time))
    start_time <- Sys.time()
    
    #Find the new tolerance such that alpha proportion of the current alive particles stay alive
    toleranceold <- tolerance
    reflevel <- alpha * tpa(toleranceold, d[,,s-1])
    tolerance<-uniroot(function(tolerance) tpa(tolerance,x=d[,,s-1])-reflevel,c(0,toleranceold))$root
    print(paste('      Step:', s, ', New tolerance is:', tolerance))
    if (tolerance<tolerancetarget){
      tolerance<-tolerancetarget
    }
    tolerancestore<-c(tolerancestore, tolerance)
    essstore<-c(essstore, 1/sum(W[,s-1]^2)) #effective sample size
    
    #compute the associated weights
    npa_old<-rowSums(d[,,s-1]<toleranceold)
    npa<-rowSums(d[,,s-1]<tolerance)
    a<-which(npa_old>0)
    b<-which(npa_old==0)
    W[a,s]<-W[a,s-1]*npa[a]/npa_old[a]
    W[b,s]<-rep(0,length(b))
    W[,s]<-W[,s]/sum(W[,s])
    
    # resample if necessary
    if ((sum(W[,s]^2)*N_T)>1){
      choice<-sample(1:Npart,Npart,replace= TRUE, prob = W[,s])
      Omega_sample[,,s]<-Omega_sample[choice,,s-1] 
      Omega_sample_phys[,,s]<-Omega_sample_phys[choice,,s-1] 
      d[,,s]<-d[choice,,s-1] 
      W[,s]<-rep(1/Npart,Npart) 
    } else { #these will be perturbed later via a MCMC step. 
      Omega_sample[,,s]<-Omega_sample[,,s-1] 
      Omega_sample_phys[,,s]<-Omega_sample_phys[,,s-1] 
      d[,,s]<-d[,,s-1] 
    }
    
    #calculate perturbation covariance based on Filippi et al., 2012
    tilda_i <- which(rowSums(d[,,s-1]<tolerance)>0) #identify old particles that fall within the new tolerance
    Omega_tilda <- Omega_sample[tilda_i,,s-1] 
    W_tilda <- W[tilda_i,s-1]
    W_tilda <- W_tilda/sum(W_tilda) #normalise weights
    
    propCOV <- matrix(0, n_x, n_x)
    for(i in 1:Npart){
      for(k in 1:length(tilda_i)){
        propCOV <- propCOV + W[i,s-1] * W_tilda[k] * ((Omega_tilda[k,]-Omega_sample[i,,s-1]) %*% t(Omega_tilda[k,]-Omega_sample[i,,s-1]))
      }
    }
    
    for(n in n_start:Npart){
      print(paste(' Step:', s, ', Particle:', n))
      if(W[n,s]>0){
        Omega_prop <- multvarNormProp(xt=Omega_sample[n,,s], propPars=propCOV) #perturb the proposal
        d_prop <-  -logTarget_new(dir = dir,Model = Model,
                              proposed = Omega_prop %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model), 
                              AlgoParams = AlgoParams) #calculate distance
        #calculate the acceptance probability:
        acc <- sum(d_prop<tolerance)/sum(d[n,,s]<tolerance) * modifyAcc(Omega_prop, Omega_sample[n,,s], Model)
        u <- runif(1)
        if(u < acc){
          Omega_sample[n,,s] <- Omega_prop
          Omega_sample_phys[n,,s] <- Omega_sample[n,,s] %>% relist(skeleton=Model$skeleton) %>% unlist()%>% Proposed2Physical(Model) %>% unlist()
          d[n,,s] <- d_prop
        }
        saveRDS(
          list(d=d, 
               omegastore=Omega_sample_phys,
               propCOV=propCOV,
               W=W,
               tolerance=toleranceold),
          paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
        )
      }
    }
    
    n_start = 1
    
    #plot the output so far:
    i_zeroweight <- which(W==0,  arr.ind=TRUE)
    plot(rep(1, Npart), apply(d[,,1], c(1), min), xlim=c(1, s), ylim=c(min(d), max(d[which(d<Inf)])))
    for (i in 2:s){
      i_s_zeroweight <- i_zeroweight[which(i_zeroweight[,2]==i),1]
      if (length(i_s_zeroweight)==0){
        points(rep(i, Npart), apply(d[,,i], c(1), min))
      }
      else{
        points(rep(i, Npart-length(i_s_zeroweight)), apply(d[-i_s_zeroweight,,i], c(1), min))
      }
    }
    
    saveRDS(
      list(d=d, 
           omegastore=Omega_sample_phys,
           propCOV=propCOV,
           W=W, 
           tolerance=tolerance),
      paste0(dir,"IIDIPUS_Results/abcsmc_",tag)
    )  
  }
}

abcsmc_2022_08_23_142322 <- readRDS(paste0(dir, "IIDIPUS_Results/abcsmc_2022-08-23_142322"))
d <- abcsmc_2022_08_23_142322$d
Npart <- AlgoParams$smc_Npart
s <- 2

i_zeroweight <- which(abcsmc_2022_08_23_142322$W==0,  arr.ind=TRUE)

plot(rep(1, Npart), apply(abcsmc_2022_08_23_142322$d[,,1], c(1), min), xlim=c(1, s), ylim=c(min(abcsmc_2022_08_23_142322$d), max(abcsmc_2022_08_23_142322$d[which(abcsmc_2022_08_23_142322$d<Inf)])))

for (i in 2:s){
  i_s_zeroweight <- i_zeroweight[which(i_zeroweight[,2]==i),1]
  if (length(i_s_zeroweight)==0){
    points(rep(i, Npart), apply(d[,,i], c(1), min))
  }
  else{
    points(rep(i, Npart-length(i_s_zeroweight)), apply(d[-i_s_zeroweight,,i], c(1), min))
  }
}
