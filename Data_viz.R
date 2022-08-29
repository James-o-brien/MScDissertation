# ---------------------------- Data Visualisation using Copernicus_Damage and UNOSAT_Damage for TC Harold ----------------------#

# Extract Environment Variables
source('RCode/GetEnv.R')
# Download and install the necessary packages:
source('RCode/GetODDPackages.R')
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
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(ggmap)
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

Damage

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
################################ Wrap the below in a function and give to Hamish/Max?


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
    color = cbPalette[2], bg = "black", pch = 23, alpha = 1, lwd = 2, # pch 2 is triangles, pch 3 is cross hairs, pch 4 is cross hairs, pch 5 is diamonds, pch 6 is upside down triangles 
    show.legend = 'point', inherit.aes = F # alpha above is [0,1] and defines the transparency of the colours, 0 being more transparent
  ) 
p + 
  labs(title = "Buildings <b style='color:#D73027'>damaged</b> in Sanma, Vanuatu",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))


# if either of the two other worst affected areas have destroyed buildings, you can colour these differently on the plot and in the title
# e.g. Buildings *damaged* and *destroyed* in Melampa, etc.

# ------------------------------------------------------ GDP Plots -----------------------------------------------------------------#
# can't find a nice way to illustrate the GDP for TC Harold, come back to just a general visualisation over say Europe when you have time

plotGDP(GDP_fji_1)
plotGDP(GDP_slb)
GDP_vut
plotGDP(GDP_vut)
GDP_fji_1_df <- as.data.frame(GDP_fji_1)
GDP_fji_2_df <- as.data.frame(GDP_fji_2)
GDP_ton_df <- as.data.frame(GDP_ton)
GDP_slb_df <- as.data.frame(GDP_slb)
GDP_vut_df <- as.data.frame(GDP_vut)
totalGDP <- merge(GDP_fji_1_df, GDP_fji_2_df)
totalGDP$X2015

EuropeGDP <- st_as_sf(EuropeGDP, coords = EuropeGDP@coords)
EuropeGDP@coords

europe_map <- get_stamenmap(
  bbox = c(-10.582033, 36.002883, 32.132810, 59.526237),
  zoom = 8, maptype = 'terrain', source = 'stamen'
) %>% ggmap()
europe_map

p2 <- europe_map +
  geom_sf(
    data = EuropeGDP, 
#    aes(EuropeGDP$X2015), 
#    color = "red", bg = "black", 
    #pch = 23, alpha = 1, lwd = 2, # pch 2 is triangles, pch 3 is cross hairs, pch 4 is cross hairs, pch 5 is diamonds, pch 6 is upside down triangles 
    #show.legend = 'point', 
    inherit.aes = F # alpha above is [0,1] and defines the transparency of the colours, 0 being more transparent
) 
p2
  

p <- sanma_sf_map + 
  geom_sf(
    data = sanma_sf, 
    #aes(size = copper), 
    color = cbPalette[2], bg = "black", pch = 23, alpha = 1, lwd = 2, # pch 2 is triangles, pch 3 is cross hairs, pch 4 is cross hairs, pch 5 is diamonds, pch 6 is upside down triangles 
    show.legend = 'point', inherit.aes = F # alpha above is [0,1] and defines the transparency of the colours, 0 being more transparent
  ) 



plotGDP<-function(GDP,zoom=7){
  mad_map <- get_stamenmap(GDP@bbox, source = "stamen", maptype = "terrain", zoom=zoom)
  p<-ggmap(mad_map) + xlab("Longitude") + ylab("Latitude")
  p+geom_contour_filled(data = as.data.frame(GDP),
                        mapping = aes(x,y,z=X2015),
                        inherit.aes = F,
                        alpha=0.8)+ 
    labs(fill = "GDP-PPP [USD-2015]")
}
?geom_contour_filled
vut_ab_0@bbox
fji_ab_0@bbox
slb_ab_0@bbox
ton_ab_0@bbox

EuropeGDP <- GetKummu(dir = dir, bbox = c(-10.582033, 36.002883, 32.132810, 59.526237))

plotGDP(EuropeGDP)

# http://bboxfinder.com/#36.002883,-10.582033,59.526237,32.132810

# http://bboxfinder.com/#44.331184,-10.784178,59.401145,15.319338
totalGDP <- GetKummu(dir = dir, bbox = c(155.39250, -22.34972, 180, -4.44522))
totalGDP2 <- GetKummu(dir = dir, bbox = c(-180, -22.34972, -173.7350, -4.44522))
plotGDP(totalGDP)
plotGDP(totalGDP2)
GetKummu<-function(dir,bbox=NULL,yr=2015L){
  
  iii<-yr-1989L
  
  # file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_PPP_30arcsec_v3.nc")
  file<-paste0(dir,"Demography_Data/SocioEconomic/KUMMU/GDP_per_capita_PPP_1990_2015_v2.nc")
  GDP<-brick(file,varname="GDP_per_capita_PPP")
  GDP<-GDP[[iii]]
  
  if(!is.null(bbox)) {
    e <- as(raster::extent(c(bbox[c(1,3,2,4)])), 'SpatialPolygons')
    crs(e) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    GDP%<>%raster::crop(e)
  }
  
  GDP%<>%as('SpatialPixelsDataFrame')
  
  return(GDP)
  
}









# ---------------------------------------- Visualising the different data types in ODDRIN ----------------------------------------- #

points <- data.frame(long = c(76.75,77.25,76.75,77.25,77.75),
                     lat = c(11.25,10.75,10.25,10.25,10.25)) %>% st_as_sf(coords=c('long','lat'), crs=4326)
cellsize = 0.15
g2 = st_make_grid(
  st_as_sfc(
    st_bbox(points) + 
      c(-cellsize/2, -cellsize/2,
        cellsize/2, cellsize/2)),
  what="polygons", cellsize=cellsize)
par(mfrow=c(1,1))

par(mar = c(0.1, 0.1, 2, 0.5), cex.main = 1.5, cex.lab = 2.5, cex.sub =2.5) 
plot(g2, main = "Spatial data: pixel versus point", font.main =18, cex.lab = 2.5) 
plot(points, pch = 3, add = T)
legend(x = "topright",
#       inset = c(-0.45, 0), # You will need to fine-tune the first
       # value depending on the windows size
       legend = c("Pixels", "Points"), 
#       lty = c(1, 2),
#       col = c(2, 3),
#       lwd = 2,
       pch = c(0, 3),
       xpd = TRUE,
        bty = "n")
labs(title = "Buildings <b style='color:#D73027'>damaged</b> in Sanma, Vanuatu",
       x = "Longitude", y = "Latitude") +
theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))








x = c(76.75,77.25,76.75,77.25,77.75)
y = c(11.25,10.75,10.25,10.25,10.25)
polygon(long, lat)

plot(1, 1, col = "white", xlab = "X", ylab = "Y")
polygon(x = c(0.7, 1.3, 1.2, 0.8),                           # X-Coordinates of polygon
        y = c(0.6, 0.8, 1.4, 1),                             # Y-Coordinates of polygon
        col = "#1b98e0")   

plot(x,y, type ="n")
pts = cbind(x, y)
polygon(pts[order(atan2(x-mean(x),y-mean(y))),])

x <- c(0.66, 0.26, 0.90, 0.06, 0.94, 0.37)
y <- c(0.99, 0.20, 0.38, 0.77, 0.71, 0.17)

xnew <- x[order(Arg(scale(x) + scale(y) * 1i))]
ynew <- y[order(Arg(scale(x) + scale(y) * 1i))]

plot(xnew, ynew, type = "n")
polygon(xnew ,ynew)


# -------------- Europe examples


library(rgdal) #v1.5-28
library(rgeos) #v.0.5-9
library(ggplot2) # 3.3.5
library(rworldmap) #plot worldmap v.1.3-6
library(dplyr) #v.1.0.7
library(sf)

#Create dataframe of coordinates that fall in Europe
coord <- data.frame(cbind(runif(1000,-15,45),runif(1000,30,75)))
colnames(coord) <- c("long","lat")

#Exlude ocean points following this post
#https://gis.stackexchange.com/questions/181586/remove-points-which-are-out-of-shapefile-or-raster-extent
URL <- "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip"
fil <- basename(URL)
if (!file.exists(fil)) download.file(URL, fil)
fils <- unzip(fil)
oceans <- readOGR(grep("shp$", fils, value=TRUE), "ne_110m_ocean",
                  stringsAsFactors=FALSE, verbose=FALSE)


europe_coord <- data.frame(long = coord$long,
                           lat = coord$lat)

coordinates(europe_coord) <- ~long+lat
proj4string(europe_coord) <- CRS(proj4string(oceans))

ocean_points <- over(europe_coord, oceans)

#Add ocean points to dataset
coord$ocean <- ocean_points$featurecla
#Exlude ocean points
europe_land <- coord %>% filter(is.na(ocean))

#Load worldmap
world <- map_data("world")
#Plot europe spatial data
ggplot() + geom_map(data = world, map = world,
                    aes(long, lat, map_id = region), color = "white", 
                    fill = "lightgray", size = 0.1) +
  geom_point(data = europe_land,aes(long, lat),
             alpha = 0.7, size = 0.05) + ylim(0,70) +
  coord_sf(xlim = c(-15, 45), ylim = c(30, 75), expand = FALSE)

bb <-
  c(
    "xmin" = 4.005,
    "xmax" = 12.005,
    "ymin" = 45.008,
    "ymax" = 51.005
  ) %>%
  sf::st_bbox() %>%
  sf::st_as_sfc() %>%
  sf::st_as_sf(crs = 4326) %>%
  sf::st_transform(crs = 4326)

europe_land_sf <- europe_land %>% 
  st_as_sf(coords = c("long", "lat"), dim = "XY") %>% 
  st_set_crs(4326)

bb_extraction <- st_intersection(bb, europe_land_sf)


ggplot() +
  geom_map(
    data = world,
    map = world,
    aes(long, lat, map_id = region),
    color = "white",
    fill = "lightgray",
    size = 0.1
  ) +
  geom_point(data = europe_land,
             aes(long, lat),
             alpha = 0.7,
             size = 0.05) + ylim(0, 70) +
  geom_sf(data = bb, colour = "black", fill = NA) +
  geom_sf(data = bb2) +
  coord_sf(xlim = c(-15, 45),
           ylim = c(30, 75),
           expand = FALSE)
?geom_sf  
?geom_polygon()

polygon(x = c(4, 12),                           # X-Coordinates of polygon
        y = c(45, 51),                             # Y-Coordinates of polygon
        col = "#1b98e0")   

bb2 <- c(4,12,45,51)
bb2 <- c(
  "xmin" = 4,
  "xmax" = 12,
  "ymin" = 45,
  "ymax" = 51
)

# ---------------------------- California polygon

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

library(ggspatial)
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = counties, color = "black", fill = cbPalette[3]) + # 6 is the winner
  coord_sf(xlim = c(-132.4, -112.1), ylim = c(28, 49), expand = FALSE) +
  scale_fill_viridis_c(trans = "sqrt", alpha = .4) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  #  coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE) +
  #  xlab("Longitude") + ylab("Latitude") +
  #  ggtitle("Observation Sites", subtitle = "(2 sites in Palm Beach County, Florida)") +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  labs(title = "Polygonal spatial data: <b style='color:#F46D43'>counties</b> of California, US",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))



###########################################################################################
################# Max vs hourly max wind speed and cumulative wind speed ##################
###########################################################################################

library(ggtext)
library(RColorBrewer)
library(ggplot2)
require(gridExtra)

BDobject <- readRDS(file = paste0(dir, "IIDIPUS_Input/BDobjects_v3/BD_TC20200404VUT_7325"))

display.brewer.pal(n = 11, name = "RdYlBu")
cbPalette <- brewer.pal(name = "RdYlBu", n = 11) # Defining the colour blind palette to be used in the plots

# Plot the relative frequencies of gradings
p <- ggplot(data.frame(BDobject@data$grading), aes(x=BDobject@data$grading)) +
  geom_bar(fill = cbPalette[9]) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal()

p +
  labs(title = "Frequency of damage gradings",
       x = "Grading", y = "Frequency") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

# Create a dataframe of all the hazmeans as columns and all the buildings as rows
hazMeanDF <- data.frame(BDobject@data[6:59])

#####################################################################
##################### Total maximum wind speed ######################
#####################################################################

#
### Get a dataframe with 1. the maximum wind speeds over each row in the hazMean dataframe and 2. the grading for each row
#

maxwindspeedDF <- data.frame(maxwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(maxwindspeedDF$maxwindspeed)){
  maxwindspeedDF$maxwindspeed[i] <- max(hazMeanDF[i,], na.rm = T)
}

#
### Violin plot for each level of grading (first all of the gradings on the same plot)
#

theme_update(plot.title = element_text(hjust = 0.5)) # Centering plot titles in ggplot

p <- ggplot(maxwindspeedDF, aes(x=grading, y=maxwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2])
p +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Max wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Now each with its own individual plot compared to the complement
#

# Damaged

Damaged <- ggplot(maxwindspeedDF, aes(x=grading=="Damaged", y=maxwindspeed, fill=grading=="Damaged")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==Damaged",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# destroyed

destroyed <- ggplot(maxwindspeedDF, aes(x=grading=="destroyed", y=maxwindspeed, fill=grading=="destroyed")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==destroyed",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# moderate

moderate <- ggplot(maxwindspeedDF, aes(x=grading=="moderate", y=maxwindspeed, fill=grading=="moderate")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==moderate",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# possible

possible <- ggplot(maxwindspeedDF, aes(x=grading=="possible", y=maxwindspeed, fill=grading=="possible")) +
  theme(legend.position="none") +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==possible",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

grid.arrange(Damaged, destroyed, moderate, possible, ncol=2, nrow=2)

#####################################################################
################## Cumulative maximum wind speed ####################
#####################################################################

#
### Get the cumulative max wind speed rather than the maximum of the maximum 
# i.e. sum each row in the hazMean dataframe rather than taking the maximum of each row

cltwindspeedDF <- data.frame(cltwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(cltwindspeedDF$cltwindspeed)){
  cltwindspeedDF$cltwindspeed[i] <- sum(hazMeanDF[i,], na.rm = T)
}

p <- ggplot(cltwindspeedDF, aes(x=grading, y=cltwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.3, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + stat_summary(fun=mean, colour = cbPalette[2])

p +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Cumulative wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Each individually
#

# Damaged

Damaged <- ggplot(cltwindspeedDF, aes(x=grading=="Damaged", y=cltwindspeed, fill=grading=="Damaged")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==Damaged",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# destroyed

destroyed <- ggplot(cltwindspeedDF, aes(x=grading=="destroyed", y=cltwindspeed, fill=grading=="destroyed")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==destroyed",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# moderate

moderate <- ggplot(cltwindspeedDF, aes(x=grading=="moderate", y=cltwindspeed, fill=grading=="moderate")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==moderate",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# possible

possible <- ggplot(cltwindspeedDF, aes(x=grading=="possible", y=cltwindspeed, fill=grading=="possible")) +
  theme(legend.position="none") +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==possible",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

grid.arrange(Damaged, destroyed, moderate, possible, ncol=2, nrow=2)
