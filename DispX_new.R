setClass("ODDpixels", 
         slots = c(dir = "character",
                   hazard = "character",
                   cIndies = "data.frame",
                   fIndies = "list",
                   IDPs = "data.frame", # includes measurement dates
                   gmax = "list",
                   alerts = "data.frame",
                   I0 = "numeric",
                   hazdates = "Date",
                   eventid = "numeric",
                   predictDisp = "data.frame",
                   modifier = "list"),
         contains = c("SpatialPixelsDataFrame"))

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


Omega <- list(Lambda1 = list(nu=0,omega=0.5),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=0.4,omega=0.6),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.2359788),
              eps = list(eps=0.01304351))

AlgoParams<-list(Np = 20, # Number of Monte Carlo particles
                 cores = 4, # Number of parallelised threads per event
                 NestedCores = 1, # How many cores are to be used inside the ODD displacement calculations?
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
                 kernel = 'loglaplace', #options are lognormal or loglaplace
                 smc_steps = 20, #Number of steps in the ABC-SMC algorithm
                 smc_Npart = 100, #Number of particles in the ABC-SMC algorithm
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


fnBUnaf <- function(D_extent){
  NB <- rep(0, length(D_extent))
  for (i in 1:length(D_extent)){
    if(D_extent[i] < 0.1){
      NB[i] <- 1  
    }
  }
  return(NB)
}

fnBDam <- function(D_extent){
  NB <- rep(0, length(D_extent))
  for (i in 1:length(D_extent)){
    if(D_extent[i] >= 0.1 & D_extent[i] <= 0.5){
      NB[i] <- 1  
    }
  }
  return(NB)
}

fnBDes <- function(D_extent){
  NB <- rep(0, length(D_extent))
  for (i in 1:length(D_extent)){
    if(D_extent[i] > 0.5){
      NB[i] <- 1  
    }
  }
  return(NB)
}

ODDpolys <- readRDS(paste0(dir,"IIDIPUS_Input/ODDpolys"))
valDF <- ODDpolys@valDF
ODDpixels %<>% DispX_new(valDF, Omega,center=Model$center, LL=F, Method = AlgoParams)
ODDpixels@predictDisp
ODDpixels@data$polygon

setGeneric("DispX_new", function(ODDpixels, valDF, Omega, center, BD_params, LL, Method, epsilon=c(0.15,0.03,0.1,0.1,0.1))
  standardGeneric("DispX_new") )
# Code that calculates/predicts the total human displacement 
setMethod("DispX_new", "ODD", function(ODDpixels, valDF, Omega, center, BD_params, LL=F,
                                       Method=list(Np=20,cores=3,cap=-300), epsilon=c(0.15,0.03,0.1,0.1,0.1)
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
  LP<-rep(1,10)
  #GetLP(ODDpixels,Omega,Params,Sinc,notnans)
  # Speed things up a little
  hrange<-grep("hazMean",names(ODDpixels),value = T)
  # Function to predict displacement per gridpoint
  CalcDam<-function(ij){
    iso3c<-ODDpixels@data$ISO3C[ij]
    # Calculate local linear predictor (NOTE: is a vector due to income distribution)
    # locallinp<-LP$dGDP$linp[LP$dGDP$ind==LP$iGDP[ij]]*LP$Plinp[ij]*LP$linp[[iso3c]] #LOOSEEND
    locallinp<-rep(1,10) #reducing parameter space while I'm figuring out the MCMC
    
    # Sample population per income distribution (Assumes 9 percentiles):
    lPopS <- SplitSamplePop(Pop=ODDpixels@data$Population[ij],Method$Np) 
    tPop <-array(0,c(3, Method$Np)) #row 1 = tDisp, #row 2 = tMort, #row 3 = tRem
    tPop[3,]=colSums(lPopS)
    for(h in hrange){
      # for(h in c(1)){
      if(is.na(ODDpixels@data[ij,h])) next
      # Resample population based on who is remaining
      ind<-tPop[3,]>0
      if(h!=hrange[1]) {
        if(sum(ind)==0) break #if no remaining population, skip modelling
        if(length(lPopS[,ind])==0) break #if no remaining population, skip modelling
        #if(sum(ind)>1) sumz<-colSums(lPopS[,ind])
        #else sumz<-sum(lPopS[,ind])
        #lPopS[,!ind]<-0
        lPopS[,ind]<-SplitSamplePop(Pop=tPop[3,ind])
      }
      # Sample hazard Intensity 
      # the uncertainty is too high... so I scale it to get some interpretable results (I know, I'm not really a statistician, I don't even have a degree, I was actually just having a look around the department when they confused me for the interviewee. I didn't have the heart to say anything. You don't hate me as much as I do)
      # I_ij<-rnorm(n = Method$Np,
      #             mean = ODDpixels@data[ij,paste0("hazMean",h)],
      #             sd = ODDpixels@data[ij,paste0("hazSD",h)]/10)
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
    
    #if no building destruction data:
    if(!BD_data_present) return(rbind(tPop[1:2,, drop=FALSE], rep(NA, Method$Np))) #return total displacement and mortality, set number of buildings destroyed to NA
    
    #otherwise, sample from the model for the number of buildings destroyed:
    #we take locallinp[5] which corresponds to locallinp for the median GDP
    nBuildings = rep(ODDpixels@data$nBuildings[ij], Method$Np)
    nBDes = rep(0, Method$Np)
    nBDam = rep(0, Method$Np)
    nBUnaf = rep(0, Method$Np)
    for (h in hrange){
      if(is.na(ODDpixels@data[ij,h])) next
        if(all(nBuildings==0)) break #if no remaining buildings, skip modelling LOOSEEND
      
      
      I_ij<-ODDpixels@data[ij,h]
      Damage <-tryCatch(fDamUnscaled(I_ij,list(I0=Params$I0, Np=Params$Np),Omega)*locallinp[5], error=function(e) NA) #calculate unscaled damage (excluding GDP)
      D_extent <- BinR(Damage, Omega$zeta)

    if(length(D_extent[D_extent < 0.1]) == 0) {
      D_extent_Unaf <- 0
    } else {
        D_extent_Unaf <- D_extent[D_extent < 0.1]
        }
    
    if(length(D_extent[D_extent <= 0.5 & D_extent >= 0.1]) == 0) {
      D_extent_Dam <- 0
    } else {
        D_extent_Dam <- D_extent[D_extent <= 0.5 & D_extent >= 0.1]
        }
    
    if(length(D_extent[D_extent > 0.5]) == 0) {
      D_extent_Des <- 0
    } else {
        D_extent_Des <- D_extent[D_extent > 0.5]
        }
    
#    D_extent_Dam <- D_extent[D_extent <= 0.5 & D_extent >= 0.1]
#    D_extent_Des <- D_extent[D_extent > 0.5]
    
    moreBUnaf = fnBUnaf(D_extent_Unaf)
    nBUnaf = nBUnaf + moreBUnaf
    
    moreBDam = fnBDam(D_extent_Dam)
    nBDam = nBDam + moreBDam
    
    moreBDes = fnBDes(D_extent_Des)
    nBDes = nBDes + moreBDes
    
    nBuildings = nBuildings - moreBUnaf - moreBDam - moreBDes    
  
    }
    return(rbind(tPop[1:2,,drop=FALSE], nBDes, nBDam, nBUnaf))
  }
  
  Dam<-array(0,c(nrow(ODDpixels),Method$Np,5)) # Dam[,,1] = Displacement, Dam[,,2] = Mortality, Dam[,,3] = Buildings Destroyed, 
  # Dam[,,4] = Buildings Damaged, Dam[,,] = Buildings Unaffected
  
  #Method$cores is equal to AlgoParams$NestedCores (changed in Model file)
  if(Method$cores>1) { Dam[notnans,,]<-aperm(simplify2array(mclapply(X = notnans,FUN = CalcDam,mc.cores = Method$cores)), perm=c(3,2,1))
  } else  Dam[notnans,,]<- aperm(simplify2array(lapply(X = notnans,FUN = CalcDam)), perm=c(3,2,1))
  
  # return(Disp)
  
  # If the IDMC estimate is foreseen to be a lower or upper bound, or a generally poor estimate
  # for(c in ODDpixels@gmax$iso3){
  #   ind<-ODDpixels@data$ISO3C==c  & !is.na(ODDpixels@data$ISO3C)
  #   Disp[ind,]%<>%qualifierDisp(qualifier = ODDpixels@gmax$qualifier[ODDpixels@gmax$iso3==c],mu = Omega$mu)
  # }
  
  funcy <- function(i, LLout = T, epsilon = AlgoParams$epsilon_min, kernel = AlgoParams$kernel, cap = AlgoParams$cap) {
    tmp <- data.frame(polygon = ODDpixels@data$polygon, IDPs = Dam[,i,1], mort = Dam[,i,2], nBDes = Dam[,i,3], nBDam = Dam[,i,4], nBUnaf = Dam[,i,5]) %>% 
      group_by(polygon) %>% summarise(disp_predictor = floor(sum(IDPs, na.rm = T)), 
                                   mort_predictor = floor(sum(mort, na.rm = T)),
                                   nBDes_predictor = floor(sum(nBDes, na.rm = T)),
                                   nBDam_predictor = floor(sum(nBDam, na.rm = T)),
                                   nBUnaf_predictor = floor(sum(nBUnaf, na.rm = T)),
                                   .groups = 'drop_last')
    tmp <- tmp[!is.na(tmp$polygon), ]
    tmp %<>% merge(valDF, by="polygon") %>% arrange(desc(disp_predictor))
    
    #print(paste(tmp$nBDes_predictor, tmp$buildDestroyed))
    #print(tmp)
    if(LLout) {
      return(LL_IDP_new(tmp, epsilon,  kernel, cap))
    }
    return(tmp)
  }
  if (LL == F & Method$Np == 1){
    ODDpixels@predictDisp <- funcy(1,LLout = F) 
    
    return(ODDpixels)
  }
  outer <- sapply(1:Method$Np, funcy, epsilon = epsilon) #data.frame(length(unique(ODDpixels@gmax$iso3)))
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
  # ODDpixels@data$Disp<-Disp[,MLE]*sum(ODDpixels@gmax$gmax)/mean(sum(Disp[,MLE])) %>% round()
  ODDpixels@data$Disp <- Dam[,MLE,1]  #should this be named data$DispPred or something instead?
  ODDpixels@data$Mort <- Dam[,MLE,2]
  ODDpixels@data$nBDes <- Dam[,MLE,3]
  ODDpixels@data$nBDam <- Dam[,MLE,4]
  ODDpixels@data$nBUnaf <- Dam[,MLE,5]
  # I know this is kind of repeating code, but I want the other function as fast as possible
  ODDpixels@predictDisp <- funcy(MLE, LLout = F) 
  
  return(ODDpixels)
  
})


Omega <- list(Lambda1 = list(nu=0.2,omega=0.7),
              Lambda2 = list(nu= 1.3, omega=0.9),
              Lambda3 = list(nu=1.5,omega=3),
              zeta = list(k=2.978697, lambda=1.405539),
              Pdens = list(M=0.02988616, k = 6.473428),
              dollar = list(M = -1.051271, k = 6.473428),
              theta = list(e=0.9),
              eps = list(eps=0.95))




plot(ODDpixels) # default plots the CIESIN population data
plot(ODDpixels[c("hazMean100")]) # plots the principle hazard intensity
# Or you can use some of the ODD class methods
p <- MakeODDPlots(ODDpixels) # plots hazard and population side-by-side
p <- plotODDy(ODDpixels) # plots hazard contour lines ontop of displaced population surface plot
# Then we can also tune the plot to our greatest desires
p <- plotODDy(ODDpixels,breakings = c(0,10,50,100,500,1000,5000,10000),bbox=c(-74.5,17.9,-72.5,19),zoomy = 9)
p
