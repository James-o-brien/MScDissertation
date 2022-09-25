#################################
### Run the ABC-SMC algorithm ###
#################################

dir<-directory<-"C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/"
setwd(directory)

source('RCode/Functions.R')
source('RCode/Model.R')
source('RCode/DispXNew.R')
source('RCode/BDXNew.R')
source('RCode/ModelChanges.R')

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

Sample_Omega <- function(Model, AlgoParams){
  #Returns a sample of Omega on the transformed space (not the physical space)
  n_x <- length(unlist(Model$links))
  sample <- rep(0, n_x)
  sample <- runif(n_x, Model$par_lb, Model$par_ub) #generate proposal on the physical space
  return(sample %>% relist(Model$skeleton)%>% Physical2Proposed(Model) %>% unlist())
}

if(is.null(AlgoParams$AllParallel)){
  if(AlgoParams$cores>4) { AlgoParams$AllParallel<-T
  } else AlgoParams$AllParallel<-F
}

multvarNormProp <- function(xt, propPars){
  # purpose : A multivariate Gaussian random walk proposal for Met-Hastings
  #           MCMC
  # inputs  : xt       - The value of the chain at the previous time step 
  #           propPars - The correlation structure of the proposal
  return(array(mvtnorm::rmvnorm(1, mean=xt, sigma=propPars),dimnames = list(names(xt))))
}

modifyAcc <- function(xNew, xPrev, Model){
  xNew %<>% relist(Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
  xPrev %<>% relist(Model$skeleton) %>% unlist() %>% Proposed2Physical(Model) %>% unlist()
  Model$acceptTrans%<>%unlist()
  index<-1:length(unlist(Model$unlinks))
  product <- 1
  for (i in index){
    product <- product * match.fun(Model$acceptTrans[[names(xNew)[i]]])(xNew[i], xPrev[i], Model$par_lb[i],Model$par_ub[i])
  }
  return(product)
}

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
      Omega_sample[n,,1] <- Sample_Omega(Model, AlgoParams)
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
        acc <- sum(d_prop<tolerance)/sum(d[n,,s]<tolerance) * Acc(Omega_prop, Omega_sample[n,,s], Model)
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

#
### Uncomment to run an example below: WARNING (this will take several days) and requires downloading the data and storing in the appropriately named folders
#

# memory.limit(56000)
# gc()
# abcSmc_delmoral_new(AlgoParams, Model)
