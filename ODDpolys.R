################
### ODDpolys ###
################

### Functions:
# - convIso3Country: converts a country's ISO3 code to the country name
# - cleanValData: filters to the inputted country iso3 code, removes NA columns, orders the rows such that highest admin levels are at the 
# bottom, formats the date columns, adds a column describing the polygon which will be used in getPolys to extract data from OSM
# - getPolys: extracts polygons from OSM given the validation data spreadsheet
# - extractIndices: extracts the indices of the rows of the coordinates in the ODDpixels object that have points contained within the polygons
# that were extracted from getPolys.
# - initialize ODDpolys: fill the slots of ODDpolys with polyIndices (a list of vectors, one for each polygon, 
# with each vector containing the indices of the ODDpixels object whose points lie within the polygon), sourceInfo (data frame made up of source 
# date, source type, and the source itself), valInfo (dataframe with 2 columns - one details the type of validation data (i.e.
# mortality, displacement, b_damaged, or b_destroyed, and the other details the corresponding quantity of validation data, such that these two
# columns should be read in conjunction with each other e.g. mortality 2 means 2 deceased).

library(tidyverse)
library(magrittr)
library(OpenStreetMap)
library(osmdata)
library(sf)
library(raster)
library(devtools)
library(parallel)
#library(doParallel)
#library(foreach)
install_github('nathanvan/parallelsugar')
library(parallelsugar)

dir<-directory<-"C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/"
setwd(directory)

packred<-F
source('RCode/GetODDPackages.R')

ValData <- read.csv(paste0(dir, "Harold_data/news_sources/tc_harold_val_data_v6.csv"), header = TRUE , na.strings = c("","NA"))

convIso3Country<-function(iso3){
  countrycode::countrycode(sourcevar = iso3,
                           origin = "iso3c",
                           destination = "country.name", warn = F)
}

cleanValData <- function(ValData, ISO3){
  
  ### Input:
  # - ValData: a spreadsheet with validation data collected from scattered and sporadic sources, including online news outlets and Twitter.
  # Each row in this spreadsheet will have a non-zero entry for either building damage, mortality, or displacement.
  # - ISO3: the desired country's ISO3 code.
  ### Output:
  # - The original data filtered to rows relating to the given ISO3 code only, and only those rows with non-NA values in one of the following 
  # columns: mortality, displacement, b_damaged, b_destroyed. Also, the date columns are formatted.
  # Also, there are two columns added, valType and valQuantity, which tell us what type of validation data is in that row and what quantity this
  # type has.

  # Filtering to Vanuatu and non NA rows in mortality or displacement or building destruction or severe damage
  countryData <- ValData %>% filter(ValData$iso3 == ISO3)
  notnans <- which(!(is.na(countryData$mortality) & is.na(countryData$displacement) & 
                       is.na(countryData$b_destroyed) & is.na(countryData$b_damaged)))
  filteredData <- countryData[notnans, ]
  
  # Reorder rows to have the most granular admin level rows (highest admin numbers) at the bottom of the dataframe - this ensures that for 
  # polygons which share coordinates (e.g. a city within a county) that the predictions produced by DispX later will be broken down into
  # the most granular polygons possible (i.e. that the city coordinates are not absorbed into the county polygon)
  orderedData <- filteredData[order(filteredData$admin_4, filteredData$admin_3, filteredData$admin_2, filteredData$admin_1, na.last = FALSE), ]

  # Add a column which describes the polygon for each row (the format is "highest admin level, country")
  orderedData$polygon <- rep(NA, nrow(orderedData))
  for(i in 1:nrow(orderedData)){
    if(!(is.na(orderedData$admin_4[i]))){
      orderedData$polygon[i] <- paste0(orderedData$admin_4[i],", ",convIso3Country(orderedData$iso3[i]))
    }
    else if(!(is.na(orderedData$admin_3[i]))){
      orderedData$polygon[i] <- paste0(orderedData$admin_3[i],", ",convIso3Country(orderedData$iso3[i]))
    }
    else if(!(is.na(orderedData$admin_2[i]))){
      orderedData$polygon[i] <- paste0(orderedData$admin_2[i],", ",convIso3Country(orderedData$iso3[i]))
    }
    else {
      orderedData$polygon[i] <- paste0(orderedData$admin_1[i],", ",convIso3Country(orderedData$iso3[i]))
    }
  }
  
  # Date formatting
  orderedData$source_date %<>% as.Date(format = "%d/%m/%Y")
  orderedData$sdate %<>% as.Date(format = "%d/%m/%Y")
  orderedData$fdate %<>% as.Date(format = "%d/%m/%Y")
  
  return(orderedData)
  
}

cleanData <- cleanValData(ValData = ValData, ISO3 = "VUT")

getPolys <- function(cleanData, adminLevel = NULL){
  
  ### Input: 
  # - cleanData: the filtered, ordered, formatted validation data spreadsheet output from cleanValData, with a column added for the polygon
  # - adminLevel (NULL, 1, 2, 3, 4): if the user wants a specific administrative level, they can input the level 1, 2, 3, 4. 
  # Alternatively, NULL is the function default which extracts at the lowest level of granularity possible (highest admin level possible). 
  ### Output: 
  # - bindPolys: a spatial polygons data frame containing the polygons for each row in the ValData spreadsheet at the lowest level of granularity 
  # possible, i.e. if admin levels 1, 2, and 3 are given, the function will output the polygon corresponding to admin level 3. Note that this
  # is the default setting (corresponding to adminLevel = NULL), but the user can change to a specific admin level (1, 2, 3, or 4) if desired.
  ### Details: 
  # - The user is warned to manually check the data if two sources have similar figures with the source dates within one week of each other.
  # - The user is warned if multiple admin levels have the same name (e.g. New York, New York) 
  
  # Extract polygons loop
  polysList <- list()
  for (i in 1:nrow(cleanData)){
    
  # Warn the user if admin level names are the same    
    if(!(is.na(cleanData$admin_2[i])) & (cleanData$admin_1[i] == cleanData$admin_2[i]) |
       !(is.na(cleanData$admin_3[i])) & (cleanData$admin_2[i] == cleanData$admin_3[i]) |
       !(is.na(cleanData$admin_4[i])) & (cleanData$admin_3[i] == cleanData$admin_4[i]))
       {warning("Some administrative levels have the same name")}
    
    # Warn the user about potential double counting
      for (j in 2:nrow(cleanData)){
        
        if(!(is.na(cleanData$mortality[i])) & !(is.na(cleanData$mortality[j])) &
           (cleanData$mortality[i] == cleanData$mortality[j])
           & (!(is.na(cleanData$source_date[i])) & !(is.na(cleanData$source_date[i])) &
           (as.numeric(cleanData$source_date[i]) - as.numeric(cleanData$source_date[j]))) <= 7)
           {warning("Some sources for mortality figures are similar to each other - check for double counting")}
        
        if(!(is.na(cleanData$displacement[i])) & !(is.na(cleanData$displacement[j])) &
           (cleanData$displacement[i] == cleanData$displacement[j])
           & (!(is.na(cleanData$source_date[i])) & !(is.na(cleanData$source_date[i])) &
           (as.numeric(cleanData$source_date[i]) - as.numeric(cleanData$source_date[j]))) <= 7)
           {warning("Some sources for displacement figures are similar to each other - check for double counting")}
        
        if(!(is.na(cleanData$b_damaged[i])) & !(is.na(cleanData$b_damaged[j])) &
          (cleanData$b_damaged[i] == cleanData$b_damaged[j])
           & (!(is.na(cleanData$source_date[i])) & !(is.na(cleanData$source_date[i])) &
          (as.numeric(cleanData$source_date[i]) - as.numeric(cleanData$source_date[j]))) <= 7)
          {warning("Some sources for building damage (severe) figures are similar to each other - check for double counting")}    
        
        if(!(is.na(cleanData$b_destroyed[i])) & !(is.na(cleanData$b_destroyed[j])) &
        (cleanData$b_destroyed[i] == cleanData$b_destroyed[j])
         & (!(is.na(cleanData$source_date[i])) & !(is.na(cleanData$source_date[i])) &
        (as.numeric(cleanData$source_date[i]) - as.numeric(cleanData$source_date[j]))) <= 7)
        {warning("Some sources for building damage (destroyed) figures are similar to each other - check for double counting")}}
    
  # Extract polygons at the lowest level of granularity possible
    if(is.null(adminLevel) & !(is.na(cleanData$admin_4[i]))){
      polysList[[i]] <- getbb(cleanData$polygon[i], format_out = "sf_polygon")}
    else if(is.null(adminLevel) & !(is.na(cleanData$admin_3[i]))){
      polysList[[i]] <- getbb(cleanData$polygon[i], format_out = "sf_polygon")}
    else if(is.null(adminLevel) & !(is.na(cleanData$admin_2[i]))){
      polysList[[i]] <- getbb(cleanData$polygon[i], format_out = "sf_polygon")}
    else {
      polysList[[i]] <- getbb(cleanData$polygon[i], format_out = "sf_polygon")}
  }
  
  bindPolys <- dplyr::bind_rows(polysList)
  bindPolys %<>% as("Spatial")
  
  for(j in 1:length(bindPolys)){
    names(bindPolys@polygons)[j] <- cleanData$polygon[j]
  }
  
  return(bindPolys)

}

polys <- getPolys(cleanData = cleanData, adminLevel = NULL) 

inPoly<-function(poly, pop, iii = 1, sumFn = "sum"){
  
  if(any(class(pop) == "SpatialPointsDataFrame") | any(class(pop) == "SpatialPixelsDataFrame")){
    coords<-pop@coords
    data<-pop@data
  } else {
    coords<-as.data.frame(pop[,c("Longitude","Latitude")])
    data<-as.data.frame(pop)
  }
  
  insidepoly<-rep(FALSE,nrow(pop))
  
  for (i in 1:length(poly@Polygons)){
    # Get rid of values outside the bounding box first
    minipoly<-rep(FALSE,length(insidepoly))
    indies<-coords[,1]>=min(poly@Polygons[[i]]@coords[,1]) &
      coords[,1]<=max(poly@Polygons[[i]]@coords[,1]) &
      coords[,2]>=min(poly@Polygons[[i]]@coords[,2]) &
      coords[,2]<=max(poly@Polygons[[i]]@coords[,2])
    # Now we only need to calculate a few points that lie inside the polygon!
    minipoly[indies]<-sp::point.in.polygon(coords[indies,1],
                                           coords[indies,2],
                                           poly@Polygons[[i]]@coords[,1],
                                           poly@Polygons[[i]]@coords[,2])>0
    # Add to the total
    insidepoly<- insidepoly | minipoly
  }
  #outer<-match.fun(sumFn)(data[insidepoly,iii],na.rm=T)
  #return(list(vals=outer,indies=insidepoly))
  return(indies = which(insidepoly == "TRUE"))
}

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

#ODDpixels_2208 <- readRDS(paste0(dir, "IIDIPUS_Input/ODDobjects3/ODD_TC20200404VUT_7325"))

ODDpixels <- readRDS(paste0(dir,'IIDIPUS_Input/ODDobjects3/ODD_TC20200404VUT_7325'))

#
### Extract buildings from OpenStreetMap and convert to Longitude/Latitude from polygon
#

ExtractOSMbuild<-function(bbox,timeout=60){
  
  obj<-opq(bbox = c(bbox),timeout = timeout)%>%add_osm_feature("building") %>%
    opq_string()%>%osmdata_sf()
  obj<-obj$osm_polygons
  inds<-st_is_valid(obj$geometry); inds[is.na(inds)]<-FALSE
  obj<-obj[inds,]
  # obj%<>%st_as_sf()
  # st_crs(obj)<-st_crs("urn:ogc:def:crs:EPSG::4326")
  #obj%<>%dplyr::select(building.levels,geometry,name)
  obj%<>%dplyr::select(geometry,name)
  
  #  obj$building.levels%<>%as.numeric()
  # sf::sf_use_s2(FALSE)
  obj$area<-as.double(st_area(st_as_sf(obj$geometry)))
  # obj$area<-vapply(1:nrow(obj),function(i) st_area(st_as_sf(obj$geometry[i])),FUN.VALUE = numeric(1))
  obj$geometry%<>%st_centroid()
  obj$Longitude<-st_coordinates(obj$geometry)[,1]
  obj$Latitude<-st_coordinates(obj$geometry)[,2]
  obj%<>%as.data.frame%>%dplyr::select(-geometry)
  
  return(obj)
  
}

nbuildings <- ExtractOSMbuild(ODDpixels@bbox)
xy <- nbuildings[,c(3,4)]
nbuildings_spdf <- SpatialPointsDataFrame(coords = xy, data = nbuildings,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#
### Aggegrate the above to get how many buildings are in each square kilometre and make this a slot in the ODDpixels object
#

ParAggnBuildings<-function(ODDpixels,arrayz,ncores=3, funcy="sum",namer="nBuildings", napop=F){
  
  funcy<-match.fun(funcy)
  
  grid<-as.data.frame(ODDpixels@grid)
  ODDpixels@data$array<-NA
  arrayz<-data.frame(Longitude=arrayz@coords[,1],
                     Latitude=arrayz@coords[,2],
                     nBuildings=rep(1, length(arrayz$name)))
  
  # Create a function that takes in input i in 1:ncores that reduces FBpop (rename as tmp) as it goes
  # and returns a list of N/ncores values
  parAGG<-function(kk){
    
    if(napop) ijs<-which(!is.na(ODDpixels@data$Population)) else ijs<-1:nrow(ODDpixels@data)
    iiis<-floor((kk-1L)*length(ijs)/ncores+1L):floor((kk)*length(ijs)/ncores)
    if(kk==ncores) iiis<-floor((kk-1L)*length(ijs)/ncores+1L):length(ijs)
    ijs<-ijs[iiis]
    
    inds<-arrayz$Longitude< (min(ODDpixels@coords[ijs,1]) - grid$cellsize[1])&
      arrayz$Longitude>=(max(ODDpixels@coords[ijs,1]) + grid$cellsize[1])&
      arrayz$Latitude< (min(ODDpixels@coords[ijs,2]) - grid$cellsize[2])&
      arrayz$Latitude>=(max(ODDpixels@coords[ijs,2]) + grid$cellsize[2])
    
    tmp<-arrayz[!inds,]
    
    output<-rep(NA,length(ijs))
    i<-1
    for (z in ijs){
      
      inds<-tmp$Longitude< (ODDpixels@coords[z,1] + grid$cellsize[1])&
        tmp$Longitude>=(ODDpixels@coords[z,1] - grid$cellsize[1])&
        tmp$Latitude< (ODDpixels@coords[z,2] + grid$cellsize[2])&
        tmp$Latitude>=(ODDpixels@coords[z,2] - grid$cellsize[2])
      output[i]<-funcy(tmp$nBuildings[inds])
      
      tmp<-tmp[!inds,]
      i<-i+1
      
    }
    
    return(output)
    
  }
  
  if(napop) indies<-which(!is.na(ODDpixels@data$Population)) else indies<-1:nrow(ODDpixels@data)
  
  ODDpixels@data$array[indies]<-unlist(mclapply(1:ncores, FUN = parAGG, mc.cores = ncores))
  
  colnames(ODDpixels@data)[ncol(ODDpixels@data)]<-namer
  
  return(ODDpixels)
  
}

ODDpixels <- ParAggnBuildings(ODDpixels, nbuildings_spdf)

#
### Aggregating predictions per polygon
#

ODDpixels@data$polygon <- rep(NA, length(ODDpixels@data$nBuildings))
for(i in 1:length(ODDpolys@pixelsIndices)){
  ODDpixels@data$polygon[ODDpolys@pixelsIndices[[i]]] <- names(ODDpolys@polygons)[i]
}
saveRDS(ODDpixels, file = paste0(dir,'IIDIPUS_Input/ODDobjects_v4/TC20200404VUT_2208_agg_4'))

#
### Extract the indices of the pixels that lie within the polygons
#

extractIndices <- function(polys, ODDobject){
  
  ### Input:
  # - polys: a spatial polygons data frame which contains all the polygons extracted from the filtered validation data (the output of getPolys)
  # - ODDobject: either a spatial pixels data frame, previously of the ODD object structure in ODDobj.R, or a spatial points data frame, 
  # previously of the BD object structure in BDobj.R
  ### Output:
  # - inPolyInds: a list of vectors, one for each polygon, with each vector containing the indices of the ODDpixels/ODDpoints object whose 
  # points lie within the polygon
  
  inPolyInds <- list()
  for(k in 1:length(polys)){
    inPolyInds[[k]] <- inPoly(polys@polygons[[k]], pop = ODDobject@coords)
  }
  return(inPolyInds)
}

setClass("ODDpolys", 
         slots = c(pixelsIndices = "list",
#                   pointsIndices = "list",
                   sourceInfo = "data.frame",
                   valDF = "data.frame"
         ),
         contains = "SpatialPolygonsDataFrame")

setMethod(f = "initialize", signature = "ODDpolys",
          definition = function(.Object, polys, ODDpixels, ODDpoints, cleanData){
            
            ### Input:
            # polys: a spatial polygons data frame containing all the extracted polygons from the filtered validation data spreadsheet (the 
            # output of getPolys)
            # ODDpixels: a spatial pixels data frame, previously of the ODD object structure in ODDobj.R
            # cleanData: the filtered and formatted validation data spreadsheet output from filterValData
            ### Output:
            # ODDpolys: an ODDpolys object which has been initialized with all of its slots
            
            if(is.null(ODDpixels)) return(.Object)
            .Object@pixelsIndices <- extractIndices(polys, ODDpixels)
#            .Object@pointsIndices <- extractIndices(polys, ODDpoints)
            .Object@sourceInfo <- as.data.frame(cleanData[c("source_date", "source_type", "source")])
            
            valDF <- as.data.frame(cleanData %>%
                                             group_by(polygon) %>%
                                             summarise(b_destroyed = sum(b_destroyed, na.rm = T),
                                                       b_damaged = sum(b_damaged, na.rm = T),
                                                       mortality = sum(mortality, na.rm = T),
                                                       displacement = sum(displacement, na.rm = T)))
            b_totalDF <- ODDpixels@data %>%
              group_by(polygon) %>%
              summarise(b_total = sum(nBuildings, na.rm = T))
            mergedDF <- merge(valDF, b_totalDF, by = "polygon")
            mergedDF$b_unaffected <- mergedDF$b_total - mergedDF$b_destroyed - mergedDF$b_damaged
            .Object@valDF <- mergedDF
            
            .Object@data <- polys@data
            .Object@polygons <- polys@polygons
            .Object@bbox <- polys@bbox
            .Object@proj4string <- crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
            
            return(.Object)  
          
            }
)
ODDpixels <- readRDS(file = paste0(dir, "IIDIPUS_Input/ODDobjects_v4/TC20200404VUT_2208_agg_4"))
ODDpolys <- new("ODDpolys", polys = polys, ODDpixels = ODDpixels, cleanData = cleanData)

###################################### ------------------------------- ideas
# ideas:
# 1. allow for country names, not just admin level names
# 2. change getPolys around to allow the user to choose which admin levels they want to go with (could make the default the most granular),
# would be pretty easy to make this an input to the function, have done the hard part with how it is at the moment.
# 3. use lapply instead of the for loop - might be quicker - test with microbenchmark or simply just system.time()
# 4. if expanding to more than Vanuatu, will have to make sure each row in the spreadsheet doesn't have more than one non-zero entry in one of the
# columns: mortality, b_damaged, b_destroyed, displacement
# 5. Put descriptions for last two functions in (for what the input will be etc.)
# 6. change merge in DispX_new to be a quicker left join or something in dplyr
# 7. run everything through from start to finish and make sure it works (you've loaded all the necessary libraries etc.)
# 8. Allow the user to choose whether they aggregate the final results of DispX to admin 1 or nationally.
