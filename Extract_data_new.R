#
### New Extract Data function
#

val_data <- read.csv(paste0(dir, "Harold_data/news_sources/tc_harold_val_data_v2.csv"), header = TRUE)

#val_data$source_date %<>%as.Date()
#val_data$sdate %<>%as.Date()
#val_data$fdate %<>%as.Date()
EOSDIS <- data.frame("date" = c("02-04-2020","03-04-2020", "04-04-2020", "05-04-2020", "06-04-2020", "07-04-2020", "08-04-2020", "09-04-2020"),
             "wind_speed" = c(35, 60, 110, 125, 135, 110, 102.5, 95)) # in knots. First row is the correct values from EOSDIS
# ----------------------------------------- will have to change this LOOSEND

95 + (110 - 95)/2 # interpolating the second to last value for wind speed that is NA
EOSDIS$wind_speed <- as.numeric(EOSDIS$wind_speed)
#EOSDIS$date <- as.Date(EOSDIS$date, format = "%d/%m/%Y")

# The above is also in a spreadsheet. Could also use UNOSAT_Population_Exposure_HAROLD 
# in C:\Users\james\OneDrive\Documents\Oxford\Dissertation\Code\ODDRIN\Harold_data\misc\UNOSAT

CM <- read.csv(paste0(dir, "Harold_data/tc_harold_IDMC.csv"), header = TRUE, fileEncoding="UTF-8-BOM")
EMDAT <- read.csv(paste0(dir, "Harold_data/tc_harold_emdat.csv"), header = TRUE, fileEncoding="UTF-8-BOM")

CM %<>% transmute(
  iso3 = ISO3, 
  hazard = ifelse(Hazard.Type == "Storm", "TC", "UK"),
  sdate = as.Date(Start.Date),
  #fdate = as.Date(""),
  #eventid = event_id,
  qualifierDisp = 'total',
  displacement = ifelse(is.na(Internal.Displacements), 0, Internal.Displacements))
EMDAT %<>% transmute(
  iso3 = ISO,
  mortality = ifelse(is.na(Total.Deaths), 0, Total.Deaths), #assume blank cells correspond to 0 deaths ??? 
  qualifierMort = 'total',
  hazard = ifelse(Disaster.Type == "Earthquake", "EQ", ifelse(Disaster.Type == "Storm", "TC", "UK")), #Earthquake, Tropical Cyclone or Unknown
  sdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'),
  #eventid = paste(Year, Seq, sep='-'), 
  fdate = as.Date(paste(End.Day, End.Month, End.Year, sep='-'), format='%d-%m-%Y'), 
  inHelix = FALSE #tracks to avoid duplicates
)
CM <- CM[,-3] # remove sdate for CM. This doesn't agree with EMDAT or with EOSDIS which both say 02/04/2020, so I'm going with the majority vote.
DamageData <- dplyr::left_join(CM, EMDAT, by = c("iso3", "hazard")) # would have to do this by dates etc. if it was more than just Harold.
DamageData$eventid <- rep(1, length(DamageData$iso3))
# Extract GDACS database on the event for further validation & alertscore benchmarking
dfGDACS<-FilterGDACS(haz = "TC", syear = min(AsYear(DamageData$sdate)), fyear=max(AsYear(DamageData$sdate)), red=T)
# Extract all building damage points
Damage <- ExtractBDfiles(dir = dir, haz = "TC")
max(Damage$Longitude)
#lhazSDF <- simulateODDSim_TC(miniDam = DamageData, I0 = 3)
#lhazdat <- list(hazard_info=list(bbox=bbox, sdate=DamageData$sdate, fdate=DamageData$fdate, NumEvents=lenny, hazard="EQ", I0=4.5, eventdates=c()))
makeODDobject <- function(miniDam, I0 = 3, Damage){ #changed I0 from 4.5 to 3 
  
  #Damage_bbox <- c(min(Damage$Longitude), min(Damage$Latitude), max(Damage$Longitude), max(Damage$Latitude)) # bbox for Damage
  #Damage_bbox_1 <- c(min(Damage$Longitude), min(Damage$Latitude), -180, max(Damage$Latitude))
  #Damage_bbox_2 <- c(max(Damage$Longitude), min(Damage$Latitude), 180, max(Damage$Latitude))
  mins <- c(-180, -20.66485)
  maxes <- c(-178.74634, -14.93190)
  column.names <- c("min", "max")
  row.names <- c("Longitude", "Latitude")
  Damage_bbox_1 <- array(c(mins, maxes), dim = c(2, 2),
                         dimnames = list(row.names, column.names))
  
  r_1 <- raster(DamageData)
  r_1 <- raster(ncol=1, nrow=8, xmn=Damage_bbox_1[1], xmx= Damage_bbox_1[3], ymn=Damage_bbox_1[2], ymx=Damage_bbox_1[4], 
                crs="+proj=longlat +datum=WGS84") #each cell is 30 arcseconds x 30 arcseconds
  
  lenny = 1 #generate number of events according to a geometric distribution
  bbox = Damage_bbox_1
  miniDam <- DamageData
  lhazdat<-list(hazard_info=list(bbox=bbox,sdate=min(miniDam$sdate),fdate=max(miniDam$fdate),
                                 NumEvents=lenny,hazard="TC", I0=3, eventdates=rep(miniDam$sdate, lenny)))
  
  #for(i in 1:lenny){
  #  hazsdf <- simulateEvent(r, I0)
  #  if(is.null(hazsdf@data$wind_speed_mean)){
  #    next
  #  }
  #  eventCalc <- function(r, I0 = 3){ 
      # Input: 
      # - An empty raster field r (over the grid of interest)
      # Output:
      # - A spatial pixel data frame hazsdf containing the hazard intensity in each grid cell
      # Details:
      # - The earthquake intensity follows a gaussian kernel centered at (0,0)
      # - The maximum magnitude varies randomly in [5, 9.3] and the spread in [15, 25]
      # - The earthquake standard deviation in each cell is random uniform in [0.8,1.1]
      maxWind <- max(EOSDIS$wind_speed)
      r_1 <- setValues(r_1, EOSDIS$wind_speed)
      #sigma <- not sure what the spread of the maximum would be?
      #maxMag = runif(1, 4.5, 10)
      #sigma = runif(1, 20, 38)
      #r <- 
      #r <- setValues(r, spatialEco::gaussian.kernel(sigma=sigma, n=r@nrows)) 
      r_1 <- r_1 * (maxWind/r_1@data@max)
      sd <- setValues(r_1, sd(EOSDIS$wind_speed)) # ----------------------------------------- will have to change this LOOSEND
      names(r_1) <- 'wind_speed_mean'
      r_1$wind_speed_sd <- sd
      hazsdf <- as(r_1, 'SpatialPixelsDataFrame')
      hazsdf <- hazsdf[hazsdf$wind_speed_mean > 3,] ############################################### changed this to 4.5 from I0
      colnames(hazsdf@coords)<-rownames(hazsdf@bbox)<-c("Longitude","Latitude")
      hazsdf@coords
      return(hazsdf)
    }
    lhazdat <- append(lhazdat, new("HAZARD",
                                   obj=hazsdf,
                                   hazard="TC",
                                   dater=min(miniDam$sdate),
                                   I0=3,
                                   alertlevel="red", # ----------------------------------------- will have to change this LOOSEND, 
                                   # see GDACS TCs alertlevels if you want to include this, but not NB, don't think this is used in the model
                                   alertscore=0))
  #}
#  if (length(lhazdat)== 1){
#    print('Simulation failed: affected region under simulated event is too small')
#    return(NULL)
#  }
#  PopSim <- simulatePopulation(r)
#  GDPSim <- simulateGDP(r)
    
    #
    ### Population
    #
    
    vut_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_VUT_shp/gadm41_VUT_0.shp")),"Spatial")
    vut_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_VUT_shp/gadm41_VUT_1.shp")),"Spatial")
    vut_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_VUT_shp/gadm41_VUT_2.shp")),"Spatial")
    fji_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_FJI_shp/gadm41_FJI_0.shp")),"Spatial")
    fji_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_FJI_shp/gadm41_FJI_1.shp")),"Spatial")
    fji_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_FJI_shp/gadm41_FJI_2.shp")),"Spatial")
    ton_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_0.shp")),"Spatial")
    ton_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_1.shp")),"Spatial")
    ton_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_2.shp")),"Spatial")
    slb_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_SLB_shp/gadm41_SLB_0.shp")),"Spatial")
    slb_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_SLB_shp/gadm41_SLB_1.shp")),"Spatial")
    slb_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_SLB_shp/gadm41_SLB_2.shp")),"Spatial")
    
    vut_pop <- GetPopulationBbox(directory = dir, density = F, bbox = vut_ab_0@bbox, lowres = F, yr = "2020")
    ton_pop <- GetPopulationBbox(directory = dir, density = F, bbox = ton_ab_0@bbox, lowres = F, yr = "2020")
    slb_pop <- GetPopulationBbox(directory = dir, density = F, bbox = slb_ab_0@bbox, lowres = F, yr = "2020")
    
    #fji_pop <- GetPopulationBbox(directory = dir, density = F, bbox = fji_ab_0@bbox, lowres = F, yr = "2020")
    
    #fji_ab_0@bbox[1] <- 180
    #fji_ab_0@bbox[3] <- -180
    #fji_ab_0@bbox[4]
    #sum(population@data, na.rm = T)
    sum(slb_pop@data, na.rm = T)
    length(subset(vut_pop@data, vut_pop@data == "NA"))
    #table(vut_pop@data)
    #bbox <- fji_ab_0@bbox
    #yr <- "2020"
    
    # x is longitude, y is latitude
    # bbox<-c(min-longitude, min-latitude, max-longitude, max-latitude)
    Damage_bbox <- c(min(Damage$Longitude), min(Damage$Latitude), max(Damage$Longitude), max(Damage$Latitude)) # bbox for Damage
    
    fji_bbox_1st_half <- c(176.22, -21.0425, 180, -12.46172)
    fji_bbox_2nd_half <- c(-180, -21.0425, -176.22, -12.46172)
    fji_pop_1 <- GetPopulationBbox(directory = dir, density = F, bbox = fji_bbox_1st_half, lowres = F, yr = "2020")
    fji_pop_2 <- GetPopulationBbox(directory = dir, density = F, bbox = fji_bbox_2nd_half, lowres = F, yr = "2020")
    sum(fji_pop_1$Population, na.rm = T)
    sum(fji_pop_2$Population, na.rm = T) # only 28000 people here - maybe just leave 2nd half out?
    sum(vut_pop$Population, na.rm = T)
    sum(slb_pop$Population, na.rm = T)
    sum(ton_pop$Population, na.rm = T)
    
    vut_pop_raster <- raster(vut_pop)
    ton_pop_raster <- raster(ton_pop)
    slb_pop_raster <- raster(slb_pop)
    fji_pop_1_raster <- raster(fji_pop_1)
    fji_pop_2_raster <- raster(fji_pop_2)
    
    total_pop <- raster::merge(vut_pop_raster, ton_pop_raster, slb_pop_raster, fji_pop_1_raster, fji_pop_2_raster)
    total_pop <- as(total_pop, "SpatialPixelsDataFrame")
    sum(total_pop@data, na.rm = T)
    
    #
    ### GDP
    #
    
    GDP_fji_1 <-  GetKummu(dir = dir, bbox = c(176.22, -21.0425, 180, -12.46172))
    GDP_fji_2 <-  GetKummu(dir = dir, bbox = c(-180, -21.0425, -176.22, -12.46172))
    GDP_vut <- GetKummu(dir = dir, bbox = vut_ab_0@bbox)
    GDP_ton <- GetKummu(dir = dir, bbox = ton_ab_0@bbox)
    GDP_slb <- GetKummu(dir = dir, bbox = slb_ab_0@bbox)
    GDP_fji_1_raster <- raster(GDP_fji_1)
    GDP_fji_2_raster <- raster(GDP_fji_2)
    GDP_vut_raster <- raster(GDP_vut)
    GDP_slb_raster <- raster(GDP_slb)
    GDP_ton_raster <- raster(GDP_ton)
    total_GDP <- raster::merge(GDP_vut_raster, GDP_fji_1_raster, GDP_fji_2_raster, GDP_slb_raster, GDP_ton_raster)
    total_GDP <- as(total_GDP, "SpatialPixelsDataFrame")
    
    #
    ### Make the ODD object (old format of ODD for the time being while I get it to work)
    #
    DamageData$gmax <- rep(max(EOSDIS$wind_speed, length(miniDam$iso3))) #------------- should rename this max wind speed instead of gmax/or better
    # ----------------- the maximum hazard intensity?
    DamageData$buildDestroyed <- c(nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "FJI")),
                                   nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "SLB")),
                                   nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "TON")),
                                   nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "VUT")))
    #nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "VUT"))
    #nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "TON"))
    DamageData$qualifierBD <- rep("total", length(DamageData$iso3))
    miniDam <- DamageData
    ODD_harold <- new('ODD_new', lhazSDF=lhazdat, DamageData=DamageData, Pop = total_pop, GDP = total_GDP)
    
#  return(ODDSim)
#}


#lhazSDF_new <- list("hazard_info" = lhazdat, lhazSDF)

#lhazdat <- list(hazard_info = list(bbox = bbox, sdate = sdate, fdate = fdate, NumEvents = 1,
#                               hazard = "TC", I0 = 3, eventdates = sdate), hazsdf) # I got I0 = 3 from GetDisaster() function-ask Hamish on this

# Get the bbox for TC harold according to UNOSAT-UNITAR Damage
#Damage_bbox <- c(min(Damage$Longitude), min(Damage$Latitude), max(Damage$Longitude), max(Damage$Latitude)) # bbox for Damage
# x is longitude, y is latitude
# bbox<-c(min-longitude, min-latitude, max-longitude, max-latitude)
Damage_bbox_1 <- c(178.74634, -20.66485, 180, -14.93190)
Damage_bbox_2 <- c(-180, -20.66485, -178.74634, -14.93190)


# -------------------------------------------------- New initalize function ------------------------------------------------------- #

lhazSDF <- lhazdat

ExtractOSMnBuildings <- function(ODD, bbox){
  #Input: Bounding box of hazard
  #Output: Number of buildings in each of the population grid cells (denoted using ij)
  buildings<-GetOSMbuildingsODD(ODD, bbox) #retrieve number of buildings using OpenStreetMaps
  nBuildings = rasterize(buildings@coords, raster(ODD["Population"]), fun='count')@data@values #LOOSEEND: ensure same rasterization as population
  nBuildings[is.na(nBuildings)] = 0
  return(nBuildings)
}

ODD <- .Object
bbox <- .Object@bbox

GetOSMbuildingsODD<-function(ODD, bbox=NULL, minnum=50, plotty=F,timeout=60){
  
  if(is.null(bbox)) bbox<-ODD@bbox
  bbox <- c(bbox[1], bbox[2], bbox[3], bbox[4])
  if(bbox[1] < -180) bbox[1] <- -180
  
  area<-(bbox[4]-bbox[2])*(bbox[3]-bbox[1])
  # If the bounding box size is too large, separate into pixels
  # I know that an area of 0.01 normally returns decent results
  if(area>0.01){
    tbuildings<-tryCatch(ExtractOSMbuild(bbox),error=function(e) NULL)
    if(is.null(tbuildings)) {
      p<-ggplot(as.data.frame(ODD),aes(Longitude,Latitude))+stat_bin_2d(drop = F,binwidth = 0.1)
      pg<-(ggplot_build(p))$data[[1]]; rm(p)
      buildings<-data.frame()
      for(i in 1:nrow(pg)){
        bbox<-as.double(pg[i,c("xmin","ymin","xmax","ymax")])
        tbuildings<-tryCatch(ExtractOSMbuild(bbox,timeout=timeout),error=function(e) NULL)
        if(is.null(tbuildings)) next
        buildings%<>%rbind(tbuildings)
      }
    }
    else buildings <- tbuildings
  } else {
    # ensure bbox has an area of at least 0.01
    bbox<-expandBbox(bbox,0.01,scaling = F)
    buildings<-ExtractOSMbuild(bbox)
  }
  # i<-1
  # while((nrow(buildings)<minnum | sum(!is.na(buildings$building.levels))<minnum) & i<10){
  #   bbox%<>%expandBbox(1.1,scaling = T)
  #   buildings%<>%rbind(ExtractOSMbuild(bbox))
  #   buildings%<>%distinct()
  #   i<-i+1
  # }
  
  return(SpatialPointsDataFrame(coords = buildings[,c("Longitude","Latitude")],
                                data = buildings[,c("building.levels","area")],
                                proj4string = crs("+proj=longlat +datum=WGS84 +ellps=WGS84")))
  
  return(buildings)
  
}

# Initialisation of the ODD object
setMethod(f="initialize", signature="ODD_new",
          # definition=function(.Object,bbox,lhazSDF,dater=NULL,dir=directory,
          definition=function(.Object, lhazSDF=NULL, DamageData=NULL, dir="./", Model=list(
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
              .Object@IDPs<-DamageData[,c("sdate","gmax","qualifierDisp")]%>%
                transmute(date=sdate,IDPs=gmax,qualifier=qualifierDisp)
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
            
            #.Object@data$nBuildings <- ExtractOSMnBuildings(ODD = .Object, bbox=bbox)
            .Object@data$nBuildings <- ExtractOSMnBuildings(ODD = .Object, bbox=bbox)
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
            GDPfactors<-InterpGDPWB(iso3c,dater)
            for (iso in iso3c){
              indie<-.Object@data$ISO3C==iso & !is.na(.Object@data$ISO3C)
              .Object@data$Population[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
              .Object@data$GDP[indie]%<>%
                multiply_by(Popfactors$factor[Popfactors$iso3==iso])
            }
            
            print("Extract country indicators - INFORM:")
            # INFORM (Joint Research Center - JRC) data:
            INFORM<-InterpINFORMdata(Model$INFORM_vars,max(dater,as.Date("2014-10-22")),iso=iso3c)
            # World Income Database (WID) data:
            if(year==AsYear(Sys.Date())) year<-AsYear(Sys.Date())-1
            print("Extract country indicators - WID:")
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





setMethod(f="initialize", signature="ODD_new",
          definition=function(.Object, lhazSDF = NULL, DamageData = NULL, Pop = NULL, GDP = NULL, dir="./", Model=list(
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
                summarise(gmax=max(gmax),
                          qualifier=ifelse(all(is.na(gmax)), NA_character_, qualifierDisp[which.max(gmax)]),
                          mortality=max(mortality),
                          qualifierMort=qualifierMort[which.max(mortality)]) ################################ changed this line
              #                       buildDestroyed=max(buildDestroyed),  ############################## 
              #                        qualifierBD = ifelse(all(is.na(buildDestroyed), NA_character_, qualifierMort[which.max(mortality)])))
              .Object@IDPs<-DamageData[,c("sdate","gmax","qualifierDisp")]%>%
                transmute(date=sdate,IDPs=gmax,qualifier=qualifierDisp)
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
            
            obj<-total_pop
            .Object@data <- obj@data
            .Object@coords.nrs <-obj@coords.nrs
            .Object@grid <-obj@grid
            .Object@grid.index <-obj@grid.index
            .Object@coords <-obj@coords
            .Object@bbox <-obj@bbox
            .Object@proj4string <-crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            
            
            
            
            
            
            
            # ---------------------------------- need to change the below to get the actual nBuildings -----------------------#
            
            
            .Object@data$nBuildings <- round(runif(1,0.2,1) * PopSim$Population + rnorm(length(PopSim$Population),0,20))
            .Object@data$nBuildings[.Object@data$nBuildings < 0] <- 0
            
            
            
            
            
            
            # Including minshake polygon per hazard event using getcontour from adehabitatMA package
            # LOOSEEND
            .Object%<>%AddHazSDF(lhazSDF)
            
            # Extract empty indices to save time
            inds<-!is.na(.Object$Population)
            
            .Object$GDP <- total_GDP
            
            
            
            
            
            
            # ------------------------------- need to change the below to get the actual ISO3C ---------------------- #
            
            
            .Object@data$ISO3C<-NA_character_
            .Object@data$ISO3C[inds]<-'ABC'
            iso3c<-unique(.Object@data$ISO3C) ; iso3c<-iso3c[!is.na(iso3c)]
            
            # ------------------------------- need to change the below to interpolate population and GDP ------------ #
            
            # Skip the interpolation of population and GDP values for the simulations
            
            
            
            
            
            
            
            
            
            
            
            
            
            # Simulate INFORM (Joint Research Center - JRC) data:
            # Sample the indicators from a uniform distribution on [0,10]
            # Set some to 0 each with a probability of 0.1 (to mimic presence of zero values in real data).
            INFORM <- data.frame(iso3='ABC', 
                                 value=runif(length(Model$INFORM_vars),0,10),
                                 variable=Model$INFORM_vars)
            INFORM_missing <- rbernoulli(nrow(INFORM), 0.1)
            INFORM$value[INFORM_missing]<-0
            
            #Simulate World Income Database (WID) data
            #Simulate a cubic distribution with max given by max_inc (WID data appears close to cubic)
            max_inc <- runif(1,0.4,0.8)
            perc <- seq(0.1,0.9,0.1)
            cubic <- perc^3+runif(1,0.3,1)*perc^2
            WID <- data.frame(
              variable=Model$WID_perc,
              iso3='ABC',
              value=cubic*max_inc/max(cubic)
            )
            
            # Bind it all together!
            .Object@cIndies<-rbind(INFORM,WID)
            .Object@fIndies<-Model$fIndies
            
            linp<-rep(list(1.),length(unique(.Object@cIndies$iso3)))
            names(linp)<-unique(.Object@cIndies$iso3)
            .Object@modifier<-linp
            
            print("Checking ODDSim values")
            checkODD(.Object)
            
            return(.Object)
          }
)



# ------------------------------------------ Below is current functions not relating to TC Harold -----------------------------------#

mergeEMDAT_Helix <- function(EMDAT, CM){
  #cannot figure out a more efficient way of doing this using dplyr/joins
  for (i in 1:NROW(CM)){
    #find a match between the CM and the EMDAT data, using a three day window?
    match_indexes = which((EMDAT[,'iso3'] == CM[i,'iso3']) & 
                            (EMDAT[,'sdate'] > (CM[i, 'sdate'] - 3)) & 
                            (EMDAT[,'fdate'] < (ifelse(is.na(CM[i, 'fdate']), CM[i, 'sdate'] + 3, CM[i, 'fdate']))))
    
    if (length(match_indexes) > 0){
      if (length(match_indexes) > 1){ #if more than one match, take the closest date
        #EMDAT rows 286 and 287 cannot be distinguished with the exception of eventid (China // 2019-09-08)
        match_index = match_indexes[which.min(abs(EMDAT[match_indexes,'sdate']-CM[i, 'sdate']))]
      } else match_index = match_indexes
      EMDAT[match_index,'inHelix'] = TRUE
      #match_eventids = which(EMDAT[,'eventid'] == EMDAT[match_index,'eventid'])
      #EMDAT[match_eventids, 'eventid'] = CM[i,'eventid'] #change event id for all events in EMDAT with the same id
      
      CM[i,'mortality'] = EMDAT[match_index, 'mortality']
      CM[i, 'qualifierMort'] = EMDAT[match_index, 'qualifierMort']
      #need to think about which start date we use if they clash.
      CM[i,'fdate'] = if_else(CM[i,'sdate'] == EMDAT[match_index, 'sdate'], EMDAT[match_index, 'fdate'], CM[i,'fdate'])
    }
  }
  
  CM %<>% rbind(EMDAT %>% filter(inHelix == FALSE) %>% transmute(
    iso3 = iso3, 
    gmax = NA,
    qualifierDisp = NULL,
    hazard = hazard, 
    sdate = sdate, 
    fdate = fdate, 
    mortality = mortality,
    eventid = gsub("-", "", eventid) %>% as.numeric(),
    qualifierMort = qualifierMort, 
    qualifierDisp=NA))
  
  return(CM)
}

GetDisplacements_new<-function(haz, saved=F, reduce=T, GIDD=F, EMDAT=F, dir="./"){ #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ changed saved=T to saved=F
  # to bypass the ExtractCheckedDisp function above since I don't have access to the file in the original line, i.e.
  # DispData<-read.csv(paste0(dir,"IIDIPUS_Input/DispData_",haz,".csv"),na.strings = "-")
  # Also changed GIDD to False to avoid that step as the function GetGIDD isn't defined anywhere
  
  if(saved) return(ExtractCheckedDisp(dir))
  
  #CM<-GetHelix(haz=haz,reduce=reduce)
  
  Dispy <- readRDS(paste0(dir,"IIDIPUS_Input/DispData_EQ_V2.Rdata"))
  CM <- Dispy %>% transmute(
    iso3 = iso3, 
    gmax = gmax, 
    hazard = 'EQ',
    sdate = as.Date(sdate),
    fdate = as.Date(fdate),
    eventid = event_id,
    qualifierDisp = 'total')
  
  
  
  if(GIDD){ #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ changed this to false at start of the function: 
    # GetGIDD is not defined anywhere in the current GitHub repo
    # Ensure we don't access anything before 2018
    CM%<>%filter(AsYear(sdate)>=2017)
    # Extract older facts from GIDD (Global Internal Displacement Database)
    GIDD<-GetGIDD(dir,haz)
    # Combine Helix and GIDD data into one database
    CM<-MergeGIDD_Helix(GIDD,CM)
    # Ensure distinct events
    CM%<>%distinct
    # Remove what we have already extracted
    rm(GIDD)
  }
  if(EMDAT){
    startRow <- 7 #ignore metadata stored on the first 6 rows
    EMDAT <- read.xlsx(paste0(dir,"/Displacement_Data/emdat.xlsx"), sheetIndex = 1, startRow = startRow)
    EMDAT %<>% transmute(
      iso3 = ISO,
      mortality = ifelse(is.na(Total.Deaths), 0, Total.Deaths), #assume blank cells correspond to 0 deaths ??? 
      qualifierMort = 'total',
      hazard = ifelse(Disaster.Type == 'Earthquake', 'EQ', 'UK'), #Earthquake or Unknown
      sdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'),
      eventid = paste(Year, Seq, sep='-'), 
      fdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'), 
      inHelix = FALSE #tracks to avoid duplicates
    )
    
    CM <- mergeEMDAT_Helix(EMDAT, CM)
    rm(EMDAT)
    
  }
  
  return(CM)
}

# ----------------------------------- Version of GetDisplacements() that works --------------------------------- #

GetDisplacements<-function(haz, saved=F, reduce=T, GIDD=F, EMDAT=T, dir="./"){ #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ changed saved=T to saved=F
  # to bypass the ExtractCheckedDisp function above since I don't have access to the file in the original line, i.e.
  # DispData<-read.csv(paste0(dir,"IIDIPUS_Input/DispData_",haz,".csv"),na.strings = "-")
  # Also changed GIDD to False to avoid that step as the function GetGIDD isn't defined anywhere
  
  if(saved) return(ExtractCheckedDisp(dir))
  
  #CM<-GetHelix(haz=haz,reduce=reduce)
  
  Dispy <- readRDS(paste0(dir,"IIDIPUS_Input/DispData_EQ_V2.Rdata"))
  CM <- Dispy %>% transmute(
    iso3 = iso3, 
    gmax = gmax, 
    hazard = 'EQ',
    sdate = as.Date(sdate),
    fdate = as.Date(fdate),
    eventid = event_id,
    qualifierDisp = 'total')
  
  
  
  if(GIDD){ #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ changed this to false at start of the function: 
    # GetGIDD is not defined anywhere in the current GitHub repo
    # Ensure we don't access anything before 2018
    CM%<>%filter(AsYear(sdate)>=2017)
    # Extract older facts from GIDD (Global Internal Displacement Database)
    GIDD<-GetGIDD(dir,haz)
    # Combine Helix and GIDD data into one database
    CM<-MergeGIDD_Helix(GIDD,CM)
    # Ensure distinct events
    CM%<>%distinct
    # Remove what we have already extracted
    rm(GIDD)
  }
  if(EMDAT){
    startRow <- 7 #ignore metadata stored on the first 6 rows
    EMDAT <- read.xlsx(paste0(dir,"/Displacement_Data/emdat.xlsx"), startRow = startRow)
    EMDAT %<>% transmute(
      iso3 = ISO,
      mortality = ifelse(is.na(Total.Deaths), 0, Total.Deaths), #assume blank cells correspond to 0 deaths ??? 
      qualifierMort = 'total',
      hazard = ifelse(Disaster.Type == 'Earthquake', 'EQ', 'UK'), #Earthquake or Unknown
      sdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'),
      eventid = paste(Year, Seq, sep='-'), 
      fdate = as.Date(paste(Start.Day,Start.Month,Start.Year, sep='-'), format='%d-%m-%Y'), 
      inHelix = FALSE #tracks to avoid duplicates
    )
    
    CM <- mergeEMDAT_Helix(EMDAT, CM)
    rm(EMDAT)
    
  }
  
  return(CM)
}
df <- GetDisplacements(haz = "EQ")
ExtractData(haz = "EQ")
colnames(EMDAT)
colnames(CM)
# ----------------------------------- Old extract data function is below --------------------------------------- #

Damage_2 <- ExtractBDfiles(dir = dir,haz = "EQ")
Damage_2
Damage <- ExtractBDfiles(dir = dir, haz = "TC")
ExtractData(haz = "EQ")
ExtractData_test <- function(haz="EQ",dir="./",extractedData=F){
  
  if(extractedData) return(paste0(dir,"IIDIPUS_Input/ODDobjects/"))
  
  # Get the human displacement data from IDMC Helix & GIDD databases and other sources filtered by hazard
  DamageData<-GetDisplacements(haz, saved=F, GIDD=F, EMDAT=T)
  # Extract GDACS database on the event for further validation & alertscore benchmarking
  dfGDACS<-FilterGDACS(haz=haz,syear=min(AsYear(DamageData$sdate)),fyear=max(AsYear(DamageData$sdate)),red=T)
  # Extract all building damage points
  Damage<-ExtractBDfiles(dir = dir,haz = haz)
  # Per event, extract hazard & building damage objects (HAZARD & BD, resp.)
  path<-data.frame()
  for (ev in unique(DamageData$eventid)){
    # Subset displacement and disaster database objects
    miniDam<-DamageData%>%filter(eventid==-1) # changed this from eventid=ev to eventid=-1
    # Set some confining dates for the rest of the data to be assigned to this event
    maxdate<-miniDam$sdate-5
    if(is.na(miniDam$fdate)) mindate<-miniDam$sdate+3 else mindate<-miniDam$fdate+3
    # GDACS subset
    miniDACS<-dfGDACS%>%filter(iso3%in%unique(miniDam$iso3) & 
                                 sdate<mindate & sdate>maxdate)
    # Match displacement and hazard data and extract hazard maps
    # HazSDF includes SpatialPixelDataFrame object of hazmean & hazsd per date 
    # (list of each, bbox-cropped to remove M < minmag)
    #lhazSDF<-tryCatch(GetDisaster(miniDam,miniDACS),error=function(e) NULL)
    #lhazSDF<-tryCatch(GetDisaster(miniDam),error=function(e) NULL) # commented this out and trying to define lhazSDF as per Simulate.R
    
    # trying to figure out why GetDisaster(miniDam) is NULL
    if(!is.null(lhazSDF)) {
      print(paste0("Warning: no hazard data found for event ", unique(miniDam$iso3),
                   " ",unique(miniDam$hazard), " ", min(miniDam$sdate) ))
      next
    }
    
    # Create the ODD object:
    lhazSDF <- simulateODDSim(miniDam = miniDam, I0 = 4.5)
    lhazdat <- list(hazard_info=list(bbox=bbox,sdate=DamageData$sdate,fdate=DamageData$fdate,NumEvents=lenny,hazard="EQ",I0=4.5,eventdates=c()))
    lhazSDF_new <- list("hazard_info" = lhazdat, lhazSDF)
    ODDy<-tryCatch(new("ODD",lhazSDF=lhazSDF_new,DamageData=DamageData),error=function(e) NULL)
    if(is.null(ODDy)) {print(paste0("ODD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
    
    # Create a unique hazard event name
    namer<-paste0(ODDy@hazard,
                  str_remove_all(as.character.Date(min(ODDy@hazdates)),"-"),
                  unique(miniDam$iso3)[1],
                  "_",ODDy@eventid)
    # Save out objects to save on RAM
    ODDpath<-paste0(dir,"IIDIPUS_Input/ODDobjects/",namer)
    saveRDS(ODDy,ODDpath)
    
    HAZARDpath<-paste0(dir,"IIDIPUS_Input/HAZARDobjects/",namer)
    saveRDS(lhazSDF,HAZARDpath)
    rm(lhazSDF)
    
    ggsave(paste0(namer,".png"), plot=plotODDyBG(ODDy),path = paste0(directory,'Plots/IIDIPUS_BG/'),width = 8,height = 5)
    
    # Building damage subset
    miniDam<-Damage%>%filter(iso3%in%unique(miniDam$iso3) & 
                               sdate<mindate & sdate>maxdate)
    # Get building damage data and filter to matched hazard events
    BDpath=NA_character_
    if(nrow(miniDam)>0) {
      # Make building damage object BD
      BDy<- tryCatch(new("BD",Damage=miniDam,ODD=ODDy),error=function(e) NULL)
      if(is.null(BDy)) {print(paste0("BD FAIL: ",ev, " ",unique(miniDam$iso3)[1]," ", unique(miniDam$sdate)[1])) ;next}
      BDpath <-paste0(dir,"IIDIPUS_Input/BDobjects/",namer)
      # Save it out!
      saveRDS(BDy, BDpath)
    }
    
    # Path to file and ID
    path%<>%rbind(data.frame(ODDpath=ODDpath,
                             BDpath=BDpath,
                             eventid=ODDy@eventid))
    # Save some RAM
    rm(ODDy,BDy,miniDam)
  }
  
  return(path)
  
}
