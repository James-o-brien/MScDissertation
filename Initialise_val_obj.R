#
### Function to initialise new ODD_objects
#

# Below, we need to change the functions CheckODD, AddGDP and AddHazSDF (and maybe multiply_by) to make it relevant to the new ODD_objs

# Initialisation of the ODD object
setMethod(f="initialize_new", signature="ODD_objs",
          # definition=function(.Object,bbox,lhazSDF,dater=NULL,dir=directory,
          definition=function(.Object,
                              lhazSDF=NULL,
                              DamageData=NULL,
                              dir="./",
                              Model=list(
            INFORM_vars=c("CC.INS.GOV.GE", # Government Effectiveness
                          "VU.SEV.AD", # Economic Dependency Vulnerability
                          "CC.INS.DRR", # Disaster Risk Reduction
                          "VU.SEV.PD", # Multi-dimensional Poverty
                          "CC.INF.PHY" # Physical Infrastructure
            ), 
            fIndies=list(CC.INS.GOV.GE=returnX, # Government Effectiveness
                         VU.SEV.AD=returnX, # Economic Dependency Vulnerability
                         CC.INS.DRR=returnX, # Disaster Risk Reduction
                         VU.SEV.PD=returnX, # Multi-dimensional Poverty
                         CC.INF.PHY=returnX, # Physical Infrastructure
                         dollar=returnX, # IncomeDistribution*GDP
                         Pdens=returnX), # IncomeDistribution*GDP
            WID_perc=   c("p10p100", # top 90% share of Income Distribution
                          "p20p100", # top 80% share of Income Distribution
                          "p30p100", # top 70% share of Income Distribution
                          "p40p100", # top 60% share of Income Distribution
                          "p50p100", # top 50% share of Income Distribution
                          "p60p100", # top 40% share of Income Distribution
                          "p70p100", # top 30% share of Income Distribution
                          "p80p100", # top 20% share of Income Distribution
                          "p90p100" # top 10% share of Income Distribution
            ))) {
            
            if(is.null(lhazSDF)) return(.Object)
            if(!class(lhazSDF[[length(lhazSDF)]])[1]=="HAZARD") return(.Object)
            
            if(lhazSDF$hazard_info$hazard=="EQ") Model$INFORM_vars%<>%c("HA.NAT.EQ")
            else if(lhazSDF$hazard_info$hazard=="TC") Model$INFORM_vars%<>%c("HA.NAT.TC")
            else if(lhazSDF$hazard_info$hazard=="FL") Model$INFORM_vars%<>%c("HA.NAT.FL")
            else stop("Not currently prepared for hazards other than EQ, TC or FL")
            
            .Object@dir<-dir
            .Object@hazard<-lhazSDF$hazard_info$hazard
            
            if(length(unique(DamageData$eventid))==1) .Object@eventid<-unique(DamageData$eventid)
            if(.Object@hazard%in%c("EQ","TC")){
              .Object@gmax<-DamageData%>%group_by(iso3)%>%
                summarise(gmax=max(gmax),qualifier=qualifierDisp[which.max(gmax)], #LOOSEEND change to be displacement specific
                          mortality=max(mortality),qualifierMort=qualifierMort[which.max(mortality)],
                          buildDestroyed=max(buildDestroyed), qualifierBD = qualifierBD[which.max(buildDestroyed)])
              .Object@IDPs<-DamageData[,c("sdate","gmax","qualifier")]%>%
                transmute(date=sdate,IDPs=gmax,qualifier=qualifier)
            } else {
              # THIS IS READY FOR THE IPM APPROACH FOR MID/LONG DURATION HAZARDS
              .Object@IDPs<-DamageData%>%group_by(sdate)%>%summarise(IDPs=max(IDPs),.groups = 'drop_last')%>%
                transmute(date=sdate,IDPs=IDPs)
              # Note that I leave .Object@gmax intentionally blank
            }
            # This bounding box is taken as the minimum region that encompasses all hazard events in HAZARD object:
            bbox<-lhazSDF$hazard_info$bbox
            dater<-min(lhazSDF$hazard_info$sdate)
            .Object@hazdates<-lhazSDF$hazard_info$eventdates
            
            year<-AsYear(dater)
            
            print("Fetching population data")
            obj<-GetPopulationBbox(.Object@dir,bbox=bbox)
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            
            
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
            # I think the structure above is specific to SpatialPixelsDataFrames,
            # so I'll have to extend to incorporate SpatialPointsDataFrames, etc.
            
            
            
            
            
            .Object@data$nBuildings <- ExtractOSMnBuildings(bbox=bbox)
            
            print("Adding hazard events")
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-!is.na(.Object$Population)
            
            print("Fetching GDP-PPP data")
            .Object%<>%AddGDP(inds)
            
            print("Filter spatial data per country")
            .Object@data$ISO3C<-NA_character_
            .Object@data$ISO3C[inds]<-coords2country(.Object@coords[inds,])
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            
            print("Interpolate population & GDP values")
            # Note there are as many values returned as iso3c codes (returns as data.frame with columns 'iso3' and 'factor')
            Popfactors<-InterpPopWB(iso3c,dater)
            
            
            
            
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
            # Note that InterpPopWB and InterpGDPWB is in GetWorldBank.R, I don't think it should need changing, it calls another function in
            # that file as well, but I think these should be fine to leave as is, but just have to make sure it works with the Cyclone Harold data
            
            GDPfactors<-InterpGDPWB(iso3c,dater)
            for (iso in iso3c){
              indie<-.Object@data$ISO3C==iso & !is.na(.Object@data$ISO3C)
              .Object@data$Population[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
              
              
              
              #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
              # Not sure if multiply_by is a function that is defined within ODDRIN or not. 
              
              
              
              
              
              
              .Object@data$GDP[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
            }
            
            print("Extract country indicators - INFORM:")
            
            
            
            
            
            
            
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            # InterpINFORMdata is in GetINFORM.R, again I don't think it'll need changing.
            
            
            # INFORM (Joint Research Center - JRC) data:
            INFORM<-InterpINFORMdata(Model$INFORM_vars,max(dater,as.Date("2014-10-22")),iso=iso3c)
            # World Income Database (WID) data:
            if(year==AsYear(Sys.Date())) year<-AsYear(Sys.Date())-1
            print("Extract country indicators - WID:")
            
            
            
            
            
            
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
            # GetWID_perc is in GetSocioEconomic.R, again I don't think it'll need changing.
            
            
            
            
            
            
            WID<-GetWID_perc(Model$WID_perc,iso3c,year)
            # Bind it all together!
            .Object@cIndies<-rbind(INFORM,WID)
            .Object@fIndies<-Model$fIndies
            
            linp<-rep(list(1.),length(unique(.Object@cIndies$iso3)))
            names(linp)<-unique(.Object@cIndies$iso3)
            .Object@modifier<-linp
            
            print("Checking ODD values")
            checkODD(.Object)
            
            return(.Object)
          }
)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Below is the functions which are called within the initialise method above which will have to be adjusted  to allow for the new validation
# object ODD_objs

checkODD<-function(object) {
  
  if(!is.null(object$Population) & any(object$Population<0,na.rm = T)) return(F) 
  if(!is.null(object$GDP) & any(object$GDP<0,na.rm = T)) return(F) 
  if(any(is.na(object@cIndies$value))) 
    print(paste0("WARNING: missing country indicator elements for ",object@cIndies$iso3[is.na(object@cIndies$value)]))
  
  TRUE
}

# Add GDP data to the ODD object by interpolating onto the grid using cubic splines
setGeneric("AddHazSDF_new", function(ODD,lhazSDF) 
  standardGeneric("AddHazSDF_new") )
setMethod("AddHazSDF_new", "ODD", function(ODD,lhazSDF){
  
  ODD@I0<-lhazSDF$hazard_info$I0
  # interpolate data onto the grid
  coords<-Genx0y0(ODD)
  lenny<-length(lhazSDF) ; start<-2
  alertscores<-alertlevels<-c() ; dates<-rep(lhazSDF$hazard_info$sdate,lenny-start+1)
  
  polysave<-array(F,dim=c(nrow(ODD),(lenny-start+1)))
  
  for (i in start:lenny){
    print(i-start+1)
    hsdf<-lhazSDF[[i]]
    # Extract detail of specific hazard
    dates[i-start+1]<-hsdf@eventdate
    alertlevels%<>%c(hsdf@alertlevel)
    alertscores%<>%c(hsdf@alertscore)
    
    if(lhazSDF$hazard_info$hazard=="TC"){
      
      layer<-with(as.data.frame(hsdf),akima::interp(x=Longitude,y=Latitude,z=mean,
                                                    xo=coords$xo,yo=coords$yo,
                                                    linear=F,extrap = F))
      
      
      layer<-c(layer$z)
      layer[layer<ODD@I0]<-NA
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
    } else {
      
      # extract polycontour of I<I0
      pcontour<-ExtractI0poly(hsdf=hsdf,ODD=ODD)
      # Find all ODD coordinates not inside polycontour
      insidepoly<-rep(F,nrow(ODD))
      if (length(unique(pcontour$id)) > 0){
        for(p in 1:length(unique(pcontour$id))){
          tcont<-filter(pcontour,id==p)
          insidepoly<-insidepoly | sp::point.in.polygon(ODD@coords[,1],
                                                        ODD@coords[,2],
                                                        tcont$Longitude,
                                                        tcont$Latitude)>0
        }
      }
      rm(tcont)
      hsdf%<>%as.data.frame
      # Interpolate BOTH MEAN & SD onto the ODD grid
      print("mean")
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=mean,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      if(all(is.na(layer))) next
      
      ODD@data[paste0("hazMean",i-start+1)]<-layer
      
      print("sd")
      layer<-with(hsdf,akima::interp(x=Longitude,y=Latitude,z=sd,
                                     xo=coords$xo,yo=coords$yo,
                                     linear=F,extrap = F))
      layer<-c(layer$z)
      layer[!insidepoly]<-NA
      ODD@data[paste0("hazSD",i-start+1)]<-layer
      
      polysave[,i-start+1]<-insidepoly
    }
  }
  
  if(lhazSDF$hazard_info$hazard=="TC") { 
    ind<-unname(apply(ODD@data,2, function(x) sum(!is.na(x))))
    ODD@data<-ODD@data[,ind>0]
  }else ODD$Population[rowSums(polysave)==0]<-NA
  
  # ODD@hazdates<-dates
  ODD@alerts<-data.frame(alertscores=alertscores,alertlevels=alertlevels)
  
  return(ODD)
  
})



AddGDP_new<-function(ODDobj,inds=NULL,GDP=NULL){
  
  if(is.null(GDP)) GDP<-GetKummu(ODDobj@dir,c(ODDobj@bbox))
  # Minimise computation if only one GDP value is found
  if(length(unique(GDP@data$GDP))==1) { #LOOSEEND: should be GDP@data$X2015 ? 
    ODDobj$GDP<-rep(unique(GDP@data$GDP),length(ODDobj@data$Population))
    ODDobj$GDP[is.na(ODDobj$Population)]<-NA
    return(ODDobj)
  }
  ODDobj$GDP<-NA
  # interpolate data onto a regular grid
  if(!is.null(inds)) {
    ODDobj$GDP[inds]<-GDP%>%raster%>%raster::extract(ODDobj@coords[inds,])
  } else {
    ODDobj$GDP<-GDP%>%raster%>%raster::extract(ODDobj@coords)
  }
  
  # The Kummu dataset is not high resolution, 
  # Therefore, sometimes Pop data exists but GDP doesn't
  # We fix this by using the nearest-neighbour extrapolation
  GDP%<>%as.data.frame()
  for (i in which(!is.na(ODDobj$Population)&is.na(ODDobj$GDP))){
    # Find the index of the closest non-NA GDP value by longitude & latitude
    # NOTE: the end term is to ensure NA's are removed
    iminnie<-which.min((ODDobj@coords[i,1]-GDP[,2])^2*(ODDobj@coords[i,2]-GDP[,3])^2 + GDP[,1]/GDP[,1])
    ODDobj$GDP[i]<-GDP[iminnie,1]
  }
  
  return(ODDobj)
  
  
}
