##########################
### Data Visualisation ###
##########################

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
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggtext)
library(gridExtra)
library(ggpubr)

dir<-directory<-"C:/Users/james/OneDrive/Documents/Oxford/Dissertation/Code/ODDRIN/"
setwd(directory)

source('RCode/GetSatDamage.R')

###############################################################
### Building damage/destruction visualisation for TC Harold ###
###############################################################

#
### Define the functions which extract buildings from OpenStreetMap
#

ExtractOSMBuildNew <- function(bbox, timeout=60, earlyret=F){
  
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

ExtractOSMBuildNewMulti <- function(bbox, timeout=60, earlyret=F){
  
  # This extract the OSM data and converts it into a convenient form for us
  obj <- opq(bbox = bbox, timeout = timeout) %>% add_osm_feature("building") %>% # opq builds an overpass query, which is read-only API that serves up custom selected parts of the OSM map data
    # add_osm_feature adds a feature to an overpass query
    opq_string() %>% osmdata_sf() # opq_string converts obj to a character string query to be submitted to the overpass API. osmdata_sf returns an OSM Overpass query as an osmdata object in sf format
  # Output for study of the general OSM object
  if(earlyret) return(obj)
  # This next step extracts only what has polygon geometry (generally buildings)
  # (Notice that the S4 objects are really easy to combine together - using rbind only)
  obj <- dplyr::select(obj$osm_polygons, geometry) %>%
    as("Spatial")
  projection(obj) <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  # Make sure the polygons are valid in R! (usually fine, but I'm paranoid...)
  err <- tryCatch(repair_geometry(obj, verbose=T), error=function(e) NA)
  if(is.na(err)) warning("Some elements of the OSM object have faulty polygons")
  
  return(obj)
  
}

#
### Load the administrative units of Vanuatu and the Copernicus/UNOSAT data 
#

vut_ab_2 <- as(sf::st_read(paste0(dir,"Harold_data/admin_boundaries/gadm41_VUT_shp/gadm41_VUT_2.shp")),"Spatial")
Damage <- ExtractBDfiles(dir = dir, haz = "TC")
Damage_bbox <- c(min(Damage$Longitude), min(Damage$Latitude), max(Damage$Longitude), max(Damage$Latitude))

#
### Sanma, Vanuatu 
#

cbPalette <- brewer.pal(name = "RdYlBu", n = 11) # Defining the colour blind palette to be used in the plots
theme_update(plot.title = element_text(hjust = 0.5)) # Centering plot titles in ggplot

sanma <- subset(vut_ab_2, vut_ab_2@data$NAME_1 == "Sanma")
sanma_obj <- ExtractOSMBuildNew(sanma@bbox, earlyret = F)
sanma_damage <- subset(Damage, (Damage$Longitude >= sanma@bbox[1] & Damage$Longitude <= sanma@bbox[3] &
                                   Damage$Latitude >= sanma@bbox[2] & Damage$Latitude <= sanma@bbox[4]))
sanma_sf <- st_as_sf(sanma_damage, coords = 2:3) # Converting sanma_damage to an S3 object
sanma_sf_Damaged <- subset(sanma_sf, sanma_sf$grading=="Damaged")
sanma_sf_destroyed <- subset(sanma_sf, sanma_sf$grading=="destroyed")

sanma_sf_map <- get_stamenmap(
  bbox = unname(st_bbox(sanma_sf)),
  zoom = 12, maptype = 'terrain', source = 'stamen'
) %>% ggmap()

sanma_plot <- sanma_sf_map +
  geom_sf(
    data = sanma_sf_Damaged, 
    #aes(size = copper),
    color = cbPalette[4], bg = "black", pch = 23, alpha = 1, lwd = 2,
    show.legend = 'point', inherit.aes = F 
  ) +
  geom_sf(
    data = sanma_sf_destroyed,
    #aes(size = copper),
    color = cbPalette[2], bg = "black", pch = 23, alpha = 1, lwd = 2, 
    show.legend = 'point', inherit.aes = F 
  ) + 
  labs(title = 
         "Buildings <b style='color:#FF8C00'>damaged</b> and <b style='color:#D73027'>destroyed</b> in Sanma, Vanuatu",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18, hjust = 0.5, margin = margin(0,0,20,0)),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))

#######################################################
### Visualising the different data types in ODDRIN  ###
#######################################################

#
### Grids and points
#

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

par("mar"=c(3,4,4,2))
par(mar = c(0.1, 0.1, 4, 0.5), cex.main = 2, cex.lab = 2.5, cex.sub =2.5) 
plot(g2, main = "Gridded spatial data vs Point spatial data", font.main =18, cex.lab = 2.5) 
plot(points, pch = 3, add = T)
legend(x = "right",
       legend = c("Gridded data", "Point data"), 
       pch = c(0, 3),
       xpd = TRUE,
       cex = 1.5,
        bty = "n")

#
### California polygon
#

theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
counties <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
counties <- subset(counties, grepl("california", counties$ID))
ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  geom_sf(data = counties, color = "black", fill = cbPalette[3]) +
  coord_sf(xlim = c(-132.4, -112.1), ylim = c(28, 49), expand = FALSE) +
  scale_fill_viridis_c(trans = "sqrt", alpha = .4) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  labs(title = "Polygonal spatial data: <b style='color:#F46D43'>counties</b> of California, US",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 24, hjust= 0.5, margin = margin(0,0,30,0)), 
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))

###########################################################################################
################# Max vs hourly max wind speed and cumulative wind speed ##################
###########################################################################################

BD <- readRDS(file = paste0(dir, "IIDIPUS_Input/BDobjects_v3/BD_TC20200404VUT_7325"))
bluesPalette <- brewer.pal(name = "Blues", n = 9)

# Plot the relative frequencies of gradings
p <- ggplot(data.frame(BD@data$grading), aes(x=BD@data$grading)) +
  geom_bar(fill = cbPalette[8], colour = "black") + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() +
  scale_x_discrete(breaks=c("Damaged","destroyed","moderate", "possible"),
                   labels=c("Damaged","Destroyed", "Moderate", "Possible")) +
  labs(title = "Frequency of building damage gradings",
       x = "Grading", y = "Frequency") +
  theme(plot.title = element_text(lineheight = 1.1, size = 24, hjust = 0.5, margin = margin(0,0,20,0)),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"),
                axis.text = element_text(size=12),
                axis.title.x = element_text(size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
                axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)))

#
### Total maximum wind speed 
#

# Create a dataframe of all the hazmeans as columns and all the buildings as rows
hazMeanDF <- data.frame(BDobject@data[6:59])

# Max wind speed DF
maxwindspeedDF <- data.frame(maxwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(maxwindspeedDF$maxwindspeed)){
  maxwindspeedDF$maxwindspeed[i] <- max(hazMeanDF[i,], na.rm = T)
}

#
### Violin plot for each level of grading
#

p <- ggplot(maxwindspeedDF, aes(x=grading, y=maxwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Max wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Cumulative maximum wind speed 
#

cltwindspeedDF <- data.frame(cltwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(cltwindspeedDF$cltwindspeed)){
  cltwindspeedDF$cltwindspeed[i] <- sum(hazMeanDF[i,], na.rm = T)
}

p <- ggplot(cltwindspeedDF, aes(x=grading, y=cltwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.3, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + stat_summary(fun=mean, colour = cbPalette[2]) +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Cumulative wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Box plots instead of violin
#

MaxWindSpeed <- ggplot(maxwindspeedDF, aes(x=grading, y=maxwindspeed)) +
  geom_boxplot(fill = cbPalette[8], colour = "black", position=position_dodge(0.8)) + 
  theme_minimal() +
  scale_x_discrete(breaks=c("Damaged","destroyed","moderate", "possible"),
                     labels=c("Damaged","Destroyed", "Moderate", "Possible")) +
  theme(legend.position="none") +
  labs(x = "Grading", y = "Total maximum wind speed - log(m/s)") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=12, margin = margin(t = 0, r = 10, b = 0, l = 0)))

CltWindSpeed <- ggplot(cltwindspeedDF, aes(x=grading, y=cltwindspeed)) +
  geom_boxplot(fill = cbPalette[8], colour = "black", position=position_dodge(0.8)) + 
  theme_minimal() +
  scale_x_discrete(breaks=c("Damaged","destroyed","moderate", "possible"),
                   labels=c("Damaged","Destroyed", "Moderate", "Possible")) +
  theme(legend.position="none") +
  labs(x = "Grading", y = "Cumulative maximum wind speed - log(m/s)") +
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.text = element_text(size=12),
        axis.title.x = element_text(size=12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=12, margin = margin(t = 0, r = 10, b = 0, l = 0)))

plot <- ggarrange(MaxWindSpeed, CltWindSpeed, ncol=2, nrow=1)
annotate_figure(plot, top = text_grob("Building damage gradings: Maximum vs Cumulative Wind Speed", size = 18)) +
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines"))


######################################
### Plot Displacement over Vanuatu ###
######################################

mad_map <- get_stamenmap(ODDy@bbox,source = "stamen",maptype = "terrain",zoom=12)
q <- ggmap(mad_map)
p <- q + geom_raster(data = as.data.frame(ODDy), aes(Longitude, Latitude, fill = Disp),
                   alpha = 0.5, interpolate = T, inherit.aes = FALSE) + coord_cartesian() +
  scale_fill_gradient2(low = "blue", mid="blue", high = "red",trans = "log",
                       labels = c(1,8,60),
                       breaks = c(1,8,60),
                       na.value = "transparent") +
  labs(fill = "Displaced Population",
    title = "Predicted displacement across Vanuatu",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18, hjust = 0.5, margin = margin(0,0,30,0)),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines"))
