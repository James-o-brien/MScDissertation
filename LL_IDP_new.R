# Get the log-likelihood for the displacement data
LL_IDP_new<-function(Y, epsilon, kernel, cap){
  LL <- 0
  k <- 10
  impacts = c('displacement', 'mortality', 'b_destroyed', 'b_damaged', 'b_unaffected') 
  predictions = c('disp_predictor', 'mort_predictor', 'nBDes_predictor', 'nBDam_predictor', 'nBUnaf_predictor')
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
  # if(impacts_observed[1] & !is.na(Y$gmax)){
  #   #LL_disp <- log(dloglap(Y$gmax+k, location.ald = log(Y$disp_predictor+k), scale.ald = epsilon, tau = 0.5, log = FALSE)/
  #   #  (1-ploglap(k, location.ald = log(Y$disp_predictor+k), scale.ald = epsilon, tau = 0.5, log = FALSE)))
  #   LL_disp <- log(dlnormTrunc(Y$gmax+k, log(Y$disp_predictor+k), sdlog=epsilon[1], min=k))
  #   LL = LL + LL_disp
  # }
  # 
  # if("mortality" %in% colnames(Y) & !is.na(Y$mortality)){
  #   #LL_mort <- log(dloglap(Y$mortality+k, location.ald = log(Y$mort_predictor+k), scale.ald = epsilon, tau = 0.5, log = FALSE)/
  #   #  (1-ploglap(k, location.ald = log(Y$mort_predictor+k), scale.ald = 0.3 * epsilon, tau = 0.5, log = FALSE)))
  #   LL_mort <- log(dlnormTrunc(Y$mortality+k, log(Y$mort_predictor+k), sdlog=epsilon[2], min=k))
  #   LL = LL + LL_mort
  # }
  # 
  # if("buildDestroyed" %in% colnames(Y) & !is.na(Y$buildDestroyed)){
  #   #LL_BD <- log(dloglap(Y$buildDestroyed+k, location.ald = log(Y$nBD_predictor+k), scale.ald = epsilon, tau = 0.5, log = FALSE)/
  #   #  (1-ploglap(k, location.ald = log(Y$nBD_predictor+k), scale.ald = epsilon, tau = 0.5, log = FALSE)))
  #   LL_BD <- log(dlnormTrunc(Y$buildDestroyed+k, log(Y$nBD_predictor+k), sdlog=epsilon[3], min=k))
  #   LL = LL + LL_BD
  # }
  #print(LL)
  if (any(is.infinite(LL))) {
    inf_id <- which(is.infinite(LL))
    LL[inf_id] <- cap
  }
  return(LL)
}


#iso_observed <- which(!is.na(Y[,impacts[1]]) & !is.na(Y[,predictions[1]]))
#LL_impact = log(dloglap(Y[iso_observed,impacts[1]]+k, location.ald = log(Y[iso_observed,predictions[1]]+k), scale.ald = epsilon[1], tau = 0.5, log = FALSE)/
#                  (1-ploglap(k, location.ald = log(Y[iso_observed,predictions[1]]+k), scale.ald = epsilon[1], tau = 0.5, log = FALSE)))
#if(is.finite(LL_impact)){
 # LL = LL + LL_impact
##} else {
#  LL = LL + cap
#}



