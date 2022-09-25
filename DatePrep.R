##########################
### Data preprocessing ###
##########################

dir<-directory<-"C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/"
setwd(directory)

library(htmltools)
library(htmlwidgets)
library(sp)
library(maptools)
library(aim.analysis)
library(spdplyr)
library(raster)
library(sf)
library(magrittr)
library(akima)
library(stringr)
source('RCode/DispX_new (1).R')
source('RCode/BDXNEW (1).R')

BD <- readRDS(paste0(dir, "IIDIPUS_Input/BDobjects_v3/BD_TC20200404VUT_7325"))

# Read in the OpenStreetMap data from before TC harold
OldOSM <- st_read(paste0(dir,"Harold_data/hotosm_vut_buildings_gpkg/hotosm_vut_buildings.gpkg"))
OldOSM %<>% as("Spatial")

# Get the centroids of the polygons as the coordinates of the buildings
labptLON <- rep(NA, length(OldOSM@polygons))
labptLAT <- rep(NA, length(OldOSM@polygons))

for(i in 1:length(OldOSM@polygons)){
  labptLON[i] <- OldOSM@polygons[[i]]@Polygons[[1]]@labpt[1]
  labptLAT[i] <- OldOSM@polygons[[i]]@Polygons[[1]]@labpt[2]
}

labptDF <- data.frame(Longitude = labptLON, Latitude = labptLAT)
labptDF$grading <- rep(NA, nrow(labptDF))

# Interpolate hazard intensity values from the ODDpixels object
ODDpixels <- readRDS(paste0(dir, "IIDIPUS_Input/ODDobjects_v4/TC20200404VUT_2208_agg_4"))
ODDpixelsCOORDS <- as.data.frame(ODDpixels@coords)
hazMeanDF <- data.frame(ODDpixels@data[4:82])

labptInterp <- labptDF

for(i in 1:length(hazMeanDF)){
  notnans <- which(!is.na(hazMeanDF[,i]))
  layer <- interpp(ODDpixelsCOORDS$Longitude[notnans], ODDpixelsCOORDS$Latitude[notnans], z = hazMeanDF[notnans,i],
                   xo = labptDF$Longitude, yo = labptDF$Latitude, linear = T)
  layer<-c(layer$z)
  labptInterp[paste0("hazMean",i)]<-layer
}

# Interpolate population values from the ODDpixels object
PopulationDF <- data.frame(ODDpixels@data$Population)

for(i in 1:length(PopulationDF)){
  notnans <- which(!is.na(PopulationDF[,i]))
  layer <- interpp(ODDpixelsCOORDS$Longitude[notnans], ODDpixelsCOORDS$Latitude[notnans], z = PopulationDF[notnans,i],
                   xo = labptDF$Longitude, yo = labptDF$Latitude, linear = F)
  layer<-c(layer$z)
}

labptInterp$Population <- layer

# Remove any rows outside the area covered by Copernicus
OutsideBbox <- which(labptInterp$Longitude > BD@bbox[3] |
                       labptInterp$Longitude < BD@bbox[1] |
                       labptInterp$Latitude > BD@bbox[4] |
                       labptInterp$Latitude < BD@bbox[2])

labptInterp <- labptInterp[-OutsideBbox,]

# Randomly sample rows in this data frame to keep 12494 (i.e. same number as in the BD object) 
dim(BD@data)
IndicesToKeep <- sample(seq(1, nrow(labptInterp)), 12494, replace = FALSE)
labptInterp <- labptInterp[IndicesToKeep,]

# Randomly sample columns such that there is the same number of hazMean columns in the BD data as the labptInterp DF
SampledHazMeans <- sample(colnames(labptInterp[c(8:73)]), 54, replace = FALSE) # keep the same amount in BD object 
KeepHazMeans <- str_sort(SampledHazMeans, numeric = TRUE)
RightSizeDF <- labptInterp[KeepHazMeans]

# get the column names the same to ensure rbind works
hrange<-grep("hazMean",names(BD),value = T)
colnames(RightSizeDF) <- hrange
RightSizeDF$grading <- labptInterp$grading
RightSizeDF$Confidence <- labptInterp$Confidence
RightSizeDF$ISO3C <- labptInterp$ISO3C
RightSizeDF$Population <- labptInterp$Population
RightSizeDF$GDP <- labptInterp$GDP

# rbind the coordinates first
labptCOORDS <- data.frame(Longitude = labptInterp$Longitude,
                          Latitude = labptInterp$Latitude)
rbindCOORDS <- rbind(BD@coords, labptCOORDS)
rbindCOORDS <- as.matrix(rbindCOORDS)
BD@coords <- rbindCOORDS

# now rbind the grading, confidence, ..., hazMean values
BD@data <- rbind(BD@data, RightSizeDF)

# save the new object out
saveRDS(BD, file=paste0(dir,"BDnew"))
