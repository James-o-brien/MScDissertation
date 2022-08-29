# Extract Environment Variables
source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')

# We install one other package that is needed to check the S4 SpatialPolygon* objects:
remotes::install_github("nstauffer/aim.analysis")
install.packages("html")
library(htmltools)
library(htmlwidgets)
library(sp)
library(maptools)
library(aim.analysis)
library(spdplyr)
library(raster)
library(sf)
library(magrittr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# These scripts are to help out James through the start of his MSc   #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# 
# Layout of the examples:
# - Extract buildings using OSM scripts from GetOSM.R
# - Plot the buildings onto an OSM map
# - Extract schools using OSM scripts from GetOSM.R
# - Exercise: extract area of the polygons in the SpatialPolygonsDataFrame
# - Plot the polygon with the points on the same OSM map
# - Find which grid points lay within the area polygon
# - Find which (ungridded) points lay within the area
# - Copernicus building damage assessment data file extract area
# - Exercise: plot the area and the buildings assessed
# - Exercises: sub-national socio-economic indicators
#   (taken from Sub-national Human Development Index database- SHDI
#    which is from the Global Data Platform organisation)
# 
# Files required from Hamish:
# - Copernicus example
# - SHDI
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

source("RCode/GetOSM.R")
# Here is the function for OSM
ExtractOSMbuild <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- rbind(dplyr::select(obj$osm_polygons, building.levels, geometry, amenity, name),
             dplyr::select(obj$osm_multipolygons, building.levels, geometry, amenity, name))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}
# Let's first check out what the object returned looks like
# Using Oxford as an example with the bounding box in the following format:
# bbox<-c(min-longitude, min-latitude, max-longitude, max-latitude)
bbox <- c(-1.313772, 51.710991, -1.178414, 51.803223)
obj <- ExtractOSMbuild(bbox, earlyret=T)
head(obj$osm_multipolygons)
# Let's check out the types of buildings we extracted:
table(c(obj$osm_polygons$amenity, obj$osm_multipolygons$amenity))
# Notice that not everything is a building!

# Now let's let the function do it's entire job
obj <- ExtractOSMbuild(bbox, earlyret=F)

# In order to plot, let's first builéd the 'basemap' (the background map you lay your plots onto)
mad_map <- get_map(bbox, source = "stamen", maptype = "toner", zoom = 9)
# Plot it all:
ggmap(mad_map) + geom_polygon(data = obj, aes(x=long, lat, group = group))
# Find schools only
schools <- obj %>% filter(amenity == "school")
ggmap(mad_map) + geom_polygon(data = schools, aes(x=long, lat, group = group))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# EXERCISE: plot only schools and hospitals

schools_hosps <- obj %>% filter(amenity == "school" | amenity == "doctors" | amenity == "clinic")
length(subset(obj$amenity, obj$amenity == "hospital")) # note that medical centres with doctors for outpatient care only should be tagged clinic, while individual doctors' office as doctors
ggmap(mad_map) + geom_polygon(data = schools_hosps, aes(x=long, lat, group = group))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# EXERCISE: create an extra column in the obj dataframe - the area of each of the polygons

#obj <- spTransform(obj, CRS("+proj=aea +lat_1=24 +lat_2=-33 +lat_0=0 +lon_0=24 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
#obj$area_2 <- rgeos::gArea(obj)
obj$area <- area(obj)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
?getbb
??getbb

Ireland <- getbb("Ireland", format_out = "sf_polygon")
Vanuatu <- getbb("vanuatu", format_out = "sf_polygon")

# Now let's extract the shapefile (sf) (polygon) of the administrative boundaries of Oxford:
Oxford <- getbb("Oxford, UK", format_out = "sf_polygon")
# Now let's convert the object to the S4 objects that we depend heavily on in ODDRIN:
Oxford %<>% as("Spatial")
Oxford@polygons
# Scripts to calculate whether a point lies within a polygon or not. This is useful because we can see what buildings are within a given admin boundary
inPoly <- function(poly, pop, iii=1, sumFn=NULL){
  
  insidepoly<-rep(FALSE, nrow(pop)) # repeat FALSE the number of rows in pop times
  
  if(class(pop) %in% c("SpatialPixelsDataFrame", "SpatialPointsDataFrame")){ # if the class of pop is one of these
    coords <- pop@coords  
  } else if(class(pop) == "SpatialPolygonsDataFrame"){
    coords <- coordinates(pop) # what is difference between pop@coords and coordinates(pop)? the sp article suggests they are the same but coordinates(x) is recommended
  } else stop("Please try inputting S4 sp spatial objects")
  
  for (i in 1:length(poly@Polygons)){
    insidepoly <- insidepoly | sp::point.in.polygon(coords[,1], # not sure what's going on here.
                                                   coords[,2],
                                                   poly@Polygons[[i]]@coords[,1],
                                                   poly@Polygons[[i]]@coords[,2])>0
  }
  
  if(!is.null(sumFn)) return(match.fun(sumFn)(pop@data[insidepoly, iii], na.rm=T))
  return(insidepoly)
  
}

vals <- sapply(1:length(Oxford@polygons), function(j) inPoly(Oxford@polygons[[j]], obj))
# How many buildings were actually contained in the administrative region of Oxford?
100*sum(!vals)/length(vals) # this gives the percentage
sum(!vals) # this gives the number of buildings within Oxford

# Get admin boundaries for Luganville:

vanuatu_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_VUT_shp/gadm40_VUT_0.shp")),"Spatial")
vanuatu_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_VUT_shp/gadm40_VUT_1.shp")),"Spatial")
vanuatu_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_VUT_shp/gadm40_VUT_2.shp")),"Spatial")

obj_3 <- ExtractOSMbuild(vanuatu_ab_2@bbox, earlyret = F)
vanuatu_ab_1
sanma <- vanuatu_ab_1 %>% filter(vanuatu_ab_1$NAME_1 == "Sanma")
sanma@bbox
ExtractOSMbuild(sanma@bbox)
df <- ExtractOSMbuild_2(sanma@bbox)
table(c(obj$osm_polygons$amenity, obj$osm_multipolygons$amenity))
table(c(df$amenity))

luganville <- vanuatu_ab_2 %>% filter(vanuatu_ab_2$NAME_2 == "Luganville")
luganville@bbox
coordinates(vanuatu_ab_2)
vanuatu_ab_2@bbox
vanuatu_ab_2$NAME_2
schools <- obj %>% filter(amenity == "school")

library(GADMTools)
gadm_getBbox(vanuatu_ab_2)
gadm_sp.loadCountries("gadm40_VUT_2.shp",
                      level = 2,
                      basefile = "C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/Harold_data/admin_boundaries/gadm40_VUT_shp/shp"
                      )

vanuatu_ab_2 <- asvanuatu_ab_2
vanuatu_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_VUT_shp/gadm40_VUT_0.shp")),"Spatial")
file.exists("C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/Harold_data/admin_boundaries/gadm40_VUT_shp/shp/gadm40_VUT_2.shp")
vanuatu_ab_2 <- gadm_sf_import_shp(dir = "C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/Harold_data/admin_boundaries/gadm40_VUT_shp/shp/", 
                                   name = "gadm40_VUT_2",
                                   level = 2,
                                   del = c("DCODE", "NAME3", "SDCODE"),
                                   renamed = c("ISO" = "COUNTRY",
                                                          "NAME_0" = "COUNTRY_LO",
                                                                          "NAME_1" = "NAME1",
                                                                          "NAME_2" = "NAME2"),
                                                              keepall = FALSE)
?as
install.packages("GADMTools")
library(GADMTools)
v_ab_2_2 <- gadm_sf_import_shp(dir = "C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/Harold_data/admin_boundaries/gadm40_VUT_shp/gadm40_VUT_0.shp",
                               level = 0)
coordinates(vanuatu_ab_2)
vanuatu_ab_2@bbox
obj_3 <-  ExtractOSMbuild(vanuatu_ab_2)

fiji_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_FJI_shp/gadm40_FJI_0.shp")),"Spatial")
fiji_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_FJI_shp/gadm40_FJI_1.shp")),"Spatial")
fiji_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm40_FJI_shp/gadm40_FJI_2.shp")),"Spatial")

tonga_ab_0 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_0.shp")),"Spatial")
tonga_ab_1 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_1.shp")),"Spatial")
tonga_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_TON_shp/gadm41_TON_2.shp")),"Spatial")

# Extract Copernicus polygon file:
CopPoly <- as(sf::st_read(paste0(dir,"COPERNICUS_Damage/EQ20210814HTI/EMSR536_AOI01_GRA_PRODUCT_r1_VECTORS_v1_vector/EMSR536_AOI01_GRA_PRODUCT_areaOfInterestA_r1_v1.shp")),"Spatial")
CopPoints <- as(sf::st_read(paste0(dir,"COPERNICUS_Damage/EQ20210814HTI/EMSR536_AOI05_GRA_PRODUCT_r1_VECTORS_v2_vector/EMSR536_AOI05_GRA_PRODUCT_builtUpP_r1_v2.shp")),"Spatial")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# EXERCISE: 1. extract all the OSM buildings within CopPoly@bbox
# 2. Then use inpoly to reduce to only the OSM buildings that are inside the polygons of CopPoly
# 3. Then try to find which buildings may have been damaged and tagged in the CopPoints object

obj_2 <- ExtractOSMbuild(CopPoly@bbox, earlyret=T)
names(obj_2$osm_multipolygons)
names(obj_2$osm_polygons)
# Let's check out the types of buildings we extracted:
table(c(obj_2$osm_polygons$amenity, obj_2$osm_multipolygons$amenity))
# Notice that not everything is a building!

ExtractOSMbuild_2 <- function(bbox, timeout=60, earlyret=F){
# Note this is an altered function of the ExtractOSMbuild since, in this case, obj$multipolygons does not contain the columns building.levels,
# amenity or name. So I use dplyr::bind_rows which adds NAs when there's no match for the column. Might be an idea to write a more general function
# where the user inputs the column names they wish to have in the final object and this is passed into the function. 
    
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- dplyr::bind_rows(dplyr::select(obj$osm_polygons, building.levels, geometry, amenity, name),
               dplyr::select(obj$osm_multipolygons, geometry))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}
df <- ExtractOSMbuild_2(luganville@bbox)
table(c(obj$osm_polygons$amenity))


obj_4$osm_polygons
# Now let's let the function do it's entire job
obj_2 <- ExtractOSMbuild_2(CopPoly@bbox, earlyret=F)
(obj_2$amenity)
# In order to plot, let's first builéd the 'basemap' (the background map you lay your plots onto)
mad_map <- get_map(CopPoly@bbox, source = "stamen", maptype = "toner", zoom = 9)
# Plot it all:
ggmap(mad_map) + geom_polygon(data = obj_2, aes(x=long, lat, group = group))
# Find fuel only
fuel <- obj_2 %>% filter(amenity == "fuel")
length(fuel)
ggmap(mad_map) + geom_polygon(data = fuel, aes(x=long, lat, group = group))

# 2.

vals_2 <- sapply(1:length(CopPoly@polygons), function(j) inPoly(CopPoly@polygons[[j]], obj_2))
# How many buildings were actually contained in the administrative region of CopPoly?
100*sum(!vals_2)/length(vals_2)

# 3.

CopPoints

# Not sure what this means here: this is a 9x9 tibble and tells us the damage grading for each of the 9 buildings; not sure how this ties into previous exercises?

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
library(magrittr)
SHDI <- xlsx::read.xlsx(paste0(dir, "Data_misc/SHDI-Database_VUT.xlsx"), sheetName = "Sheet1") # Human Development Index Database for Vanuatu
SHDI  %<>%  filter(year == max(SHDI$year) & level == "Subnat") # Filtering for subnational level 
# Now extract the population data for Vanuatu (source - WorldPop constrained):
pop <- raster::raster(paste0(dir, "Data_misc/vut_ppp_2020_UNadj_constrained.tif")) %>% as("SpatialPixelsDataFrame") 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# EXERCISE: 4. extract all the OSM area polygons for each 'region' of SHDI
# 5. Then find the population values in 'pop' that lie within each region
# 6. Finally, make a new column in the pop object that stores each of the socio-economic variables

SHDI

# 4. 

# Now let's extract the shapefile (sf) (polygon) of the administrative boundaries of the different regions of Vanuatu:
?getbb
Tafea <- getbb("Tafea, Vanuatu", format_out = "sf_polygon")

Shefa <- getbb("Shefa, Vanuatu", format_out = "sf_polygon")
Malampa <- getbb("Malampa, Vanuatu", format_out = "sf_polygon")
Penama <- getbb("Penama, Vanuatu", format_out = "sf_polygon")
Sanma <- getbb("Sanma, Vanuatu", format_out = "sf_polygon")
Torba <- getbb("Torba, Vanuatu", format_out = "sf_polygon")
Port_Vila <- getbb("Port Vila", format_out = "sf_polygon") # No polygonal boundary for Port_Vila
Luganville <- getbb("Luganville, Sanma, Vanuatu", format_out = "sf_polygon") # No polygonal boundary for Luganville

# Now let's convert the object to the S4 objects that we depend heavily on in ODDRIN:
Tafea %<>% as("Spatial")
Shefa %<>% as("Spatial")
Malampa %<>% as("Spatial")
Penama %<>% as("Spatial")
Sanma %<>% as("Spatial")
Torba %<>% as("Spatial")

# 5.

pop@coords
Sanma@bbox
Tafea$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 168.35376 & pop@coords[ ,1] <= 170.44998 & pop@coords[, 2] >= -20.46274 & pop@coords[, 2] <= -18.27167))
Shefa$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 167.46060 & pop@coords[ ,1] <= 169.3133 & pop@coords[, 2] >= -18.35001 & pop@coords[, 2] <= -16.5000))
Malampa$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 166.73228 & pop@coords[ ,1] <= 168.73395 & pop@coords[, 2] >= -16.97992 & pop@coords[, 2] <= -15.76333))
Penama$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 167.430 & pop@coords[ ,1] <= 168.5515 & pop@coords[, 2] >= -16.095 & pop@coords[, 2] <= -14.4750))
Sanma$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 166.3355 & pop@coords[ ,1] <= 167.55333 & pop@coords[, 2] >= -15.9432 & pop@coords[, 2] <= -14.13796))
Torba$pop <- sum(subset(pop@data, pop@coords[ ,1] >= 166.3377 & pop@coords[ ,1] <= 168.29893 & pop@coords[, 2] >= -14.6600 & pop@coords[, 2] <= -12.87138))
sum(pop@data)
sum(Tafea$pop, Shefa$pop, Malampa$pop, Penama$pop, Sanma$pop, Torba$pop) # Not sure above is correct since the sense checks don't pass (i.e. from looking at population figures in SDHI,
# googling population figures, and adding up the populations of each individual sub region (this comes out to be more than the total population in pop)). Some of this could be due to not
# being able to get bbox for Port Vila and Luganville. I think Luganville is contained within Sanma? But Port Vila is the capital, so this is surprising.

# Adding a region column whose value is based on the bbox for each of the regions. 

pop$region <- ifelse((pop@coords[ ,1] >= 168.35376 & pop@coords[ ,1] <= 170.44998 & pop@coords[, 2] >= -20.46274 & pop@coords[, 2] <= -18.27167),
                     "Tafea",
                     ifelse((pop@coords[ ,1] >= 167.46060 & pop@coords[ ,1] <= 169.3133 & pop@coords[, 2] >= -18.35001 & pop@coords[, 2] <= -16.5000),
                            "Shefa",
                     ifelse((pop@coords[ ,1] >= 166.73228 & pop@coords[ ,1] <= 168.73395 & pop@coords[, 2] >= -16.97992 & pop@coords[, 2] <= -15.76333),
                            "Malampa",
                     ifelse((pop@coords[ ,1] >= 167.430 & pop@coords[ ,1] <= 168.5515 & pop@coords[, 2] >= -16.095 & pop@coords[, 2] <= -14.4750),
                            "Penama",
                     ifelse((pop@coords[ ,1] >= 166.3355 & pop@coords[ ,1] <= 167.55333 & pop@coords[, 2] >= -15.9432 & pop@coords[, 2] <= -14.13796),
                            "Sanma",
                     ifelse((pop@coords[ ,1] >= 166.3377 & pop@coords[ ,1] <= 168.29893 & pop@coords[, 2] >= -14.6600 & pop@coords[, 2] <= -12.87138),
                            "Torba",
                            "NA")))))
                     )

# Alternative to above: each row will have the region and population
# within that bbox. Can just sum like so to get the total population in each region:
sum(subset(pop@data$vut_ppp_2020_UNadj_constrained, pop@data$region == "Tafea"))

# 6.

# Making a new column in pop that stores the socio-economic variables:

pop_3 <- merge(pop, SHDI, by.x = "region")

#$$$$$$$$$$$$$$$$$$$$$ Misc deleted from above $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

head(pop@data)
max(pop@coords[, 1])
length(subset(pop, pop@coords[ ,1] >= 168.35376 & pop@coords[ ,1] <= 170.44998 & pop@coords[, 2] >= -20.46274 & pop@coords[, 2] <= -18.27167))

SHDI
pop_2 <- merge(SHDI, pop)
pop_2 %<>% as("SpatialPixelsDataFrame")


# In order to plot, let's first builéd the 'basemap' (the background map you lay your plots onto)
mad_map <- get_map(bbox, source = "stamen", maptype = "toner", zoom = 9)
# Plot it all:
ggmap(mad_map) + geom_polygon(data = obj, aes(x=long, lat, group = group))
# Find schools only
schools <- obj %>% filter(amenity == "school")
ggmap(mad_map) + geom_polygon(data = schools, aes(x=long, lat, group = group))


# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #
# From here, I attempt to make some plots of Vanuatu, Fiji, The Solomon Islands, and Tonga, and highlight buildings damaged on these plots #
# I also calculate the number of buildings damaged in the admin boundaries of each of the countries.                                       #
# I start with trying to plot the buildings damaged in espiritu santo and colour according to population
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ #

# Extract Copernicus polygon file:
cop_poly_es <- as(sf::st_read(paste0(dir,"Harold_data/copernicus/espiritu_santo/vector_package/EMSR434_AOI05_GRA_PRODUCT_r2_VECTORS_v1_vector/EMSR434_AOI05_GRA_PRODUCT_areaOfInterestA_r2_v1.shp")), "Spatial")
cop_points_es <- as(sf::st_read(paste0(dir,"Harold_data/copernicus/espiritu_santo/vector_package/EMSR434_AOI05_GRA_PRODUCT_r2_VECTORS_v1_vector/EMSR434_AOI05_GRA_PRODUCT_builtUpP_r2_v1.shp")), "Spatial")



obj_4 <- ExtractOSMbuild_2(cop_poly_es@bbox, earlyret=F)
names(obj_4$osm_multipolygons)
names(obj_4$osm_polygons)
# Let's check out the types of buildings we extracted:
table(c(obj_4$osm_polygons$name))
# Notice that not everything is a building!
obj_4$osm_polygons$name
obj_4$osm_multipolygons

ExtractOSMbuild_3 <- function(bbox, timeout=60, earlyret=F){
  # Note this is an altered function of the ExtractOSMbuild since, in this case, obj$multipolygons does not contain the columns building.levels,
  # amenity or name. So I use dplyr::bind_rows which adds NAs when there's no match for the column. Might be an idea to write a more general function
  # where the user inputs the column names they wish to have in the final object and this is passed into the function. 
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- as(dplyr::select(obj$osm_polygons, osm_id, name, building, height, source, geometry), "Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

obj_4 <- ExtractOSMbuild_3(cop_poly_es@bbox, earlyret=F)
obj_4$osm_polygons
# In order to plot, let's first builéd the 'basemap' (the background map you lay your plots onto)
mad_map <- get_map(obj_4@bbox, source = "google", maptype = "satellite", zoom = 14)
?get_map
plot(obj_4, add = TRUE)
# Plot it all:
ggmap(mad_map) + geom_polygon(data = obj_4, aes(x=long, lat, group = group))
# Find schools only
schools <- obj %>% filter(amenity == "school")
ggmap(mad_map) + geom_polygon(data = schools, aes(x=long, lat, group = group))

install.packages("tmap") # install the CRAN version
library(tmap)
vignette("tmap-nutshell")
# A couple of basic plots show the package's intuitive syntax and attractive default parameters.
??qtm
qtm(shp = lnd, fill = "Partic_Per", fill.palette = "-Blues") # not shown
qtm(shp = lnd, fill = c("Partic_Per", "Pop_2001"), fill.palette = "Blues", ncol = 2)

library(ggplot2)
p <- ggplot(obj_4@data, aes(Partic_Per, Pop_2001))
cop_points_es[1,]@bbox
cop_points_es@bbox
?get_googlemap
install.packages("pacman")
library(pacman)
install.packages("BiocManager")
library(BiocManager)
install.packages("rstudioapi")
library(rstudioapi)

??register_google
ggmap::register_google(key = "AIzaSyCH0baB0F6OzcyvHFrBWm9zI-tqwI51rYQ")
pacman::p_load(ggmap, osmdata, get_googlemap)
??get_google
library(ggmap)
map <- get_googlemap(centre = c(lon = -16.4903, lat = 167.7968), zoom = )
?get_googlemap
ggmap(get_stamenmap(bbox = obj_4@bbox), extent = 'normal', zoom = 20, maptype = 'watercolor')
?ggmap
map <- get_stamenmap(bbox = obj_4@bbox)

# contour overlay
ggmap(get_map(maptype = "satellite"), extent = "device") +
  stat_density2d(aes(x = lon, y = lat, colour = class), data = as.data.frame(cop_points_es), bins = 5)


# adding additional content
library(grid)
baylor <- get_map("one bear place, waco, texas", zoom = 15, maptype = "satellite")
ggmap(baylor)



?get_googlemap
obj_4@bbox
coordinates(obj_4)
cop_points_es@coords
coordinates(cop_poly_es)
ggmap(map, extent = 'device')

# From here, we should add a 




# ---------------------------- Data Visualisation using Copernicus_Damage and UNOSAT_Damage for TC Harold ----------------------#

source("RCode/GetOSM.R")
library(ggplot2)
library(ggtext)
library(RColorBrewer)

# Here is the function for OSM
ExtractOSMbuild <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- rbind(dplyr::select(obj$osm_polygons, building.levels, geometry, amenity, name),
               dplyr::select(obj$osm_multipolygons, building.levels, geometry, amenity, name))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

ExtractOSMbuild_no_bls <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- rbind(dplyr::select(obj$osm_polygons, geometry, amenity, name),
               dplyr::select(obj$osm_multipolygons, geometry, amenity, name))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

ExtractOSMbuild_no_bls_am <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- rbind(dplyr::select(obj$osm_polygons, geometry, name),
               dplyr::select(obj$osm_multipolygons, geometry, name))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

ExtractOSMbuild_no_bls_am_name <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- rbind(dplyr::select(obj$osm_polygons, geometry),
               dplyr::select(obj$osm_multipolygons, geometry))%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

ExtractOSMbuild_no_bls_am_name_no_multipolygons <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- dplyr::select(obj$osm_polygons, geometry)%>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}


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

# x is longitude, y is latitude
# bbox<-c(min-longitude, min-latitude, max-longitude, max-latitude)
Damage_bbox <- c(min(Damage$Longitude), min(Damage$Latitude), max(Damage$Longitude), max(Damage$Latitude)) # bbox for Damage

sanma <- subset(vut_ab_2, vut_ab_2@data$NAME_1 == "Sanma")
sanma_obj <- ExtractOSMbuild_no_bls_am_name(sanma@bbox, earlyret = F)
mad_map <- get_map(sanma@bbox, source = "stamen", maptype = "watercolor", zoom = 9)
ggmap(mad_map) + geom_polygon(data = sanma_obj, aes(x=long, lat, group = group)) 
sanma_damage <- subset(Damage, (Damage$Longitude >= sanma@bbox[1] & Damage$Longitude <= sanma@bbox[3] &
                                  Damage$Latitude >= sanma@bbox[2] & Damage$Latitude <= sanma@bbox[4]))
sanma_sf <- st_as_sf(sanma_damage, coords = 2:3) # Converting sanma_damage to an S3 object
nrow(subset(Damage, Damage$grading == "destroyed" & Damage$iso3 == "VUT"))

# ----------------------------- Plots of 3 worst affected regions in terms of building damage from Harold ----------------------------#

#--------------------------------# 
### Luganville, Sanma, Vanuatu ###
#--------------------------------#

display.brewer.pal(n = 11, name = "RdYlBu")
cbPalette <- brewer.pal(name = "RdYlBu", n = 11) # Defining the colour blind palette to be used in the plots
theme_update(plot.title = element_text(hjust = 0.5)) # Centering plot titles in ggplot
sanma_sf_map <- get_stamenmap(
    bbox = unname(st_bbox(sanma_sf)),
    zoom = 12, maptype = 'terrain', source = 'stamen'
  ) %>% ggmap()
p <- sanma_sf_map + 
  geom_sf(
    data = sanma_sf, 
    #aes(size = copper), 
    color = "black", bg = cbPalette[2], pch = 23, alpha = 1, lwd = 2, # pch 2 is triangles, pch 3 is cross hairs, pch 4 is cross hairs, pch 5 is diamonds, pch 6 is upside down triangles 
    show.legend = 'point', inherit.aes = F
  ) 
p + 
  labs(title = "Buildings <b style='color:#D73027'>damaged</b> in Sanma, Vanuatu",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))


