#
### 1. Spatial Pixels Data Frame (formerly the ODD object in ODDobj.R)
#
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

#
### 2. Spatial Points Data Frame (formerly the BD object in BDObj.R) 
#
setClass("ODDpoints", 
         slots = c(hazard = "character",
                   cIndies = "data.frame",
                   fIndies = "list",
                   I0 = "numeric",
                   hazdates = "Date",
                   eventid = "numeric",
                   coefs = "numeric",
                   buildingsfile = "character",
                   modifier = "list"),
         contains = "SpatialPointsDataFrame")

#
### 3. Spatial Polygons Data Frame (not formerly allowed for in ODDRIN). 
#
setClass("ODDpolys", 
         slots = c(polyIndices = "numeric",
                   sourceInfo = "character",
                   buildingsDam = "numeric"),
         contains = "SpatialPolygonsDataFrame")

#
### 5. Full validation object
#
setClassUnion("ODDobjs", c("ODDpixels", "ODDpoints", "ODDpolys"))
setClass("ODD", 
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
                   modifier = "list",
                   coefs = "numeric",
                   buildingsfile = "character"
         ),
         contains = "ODDobjs")

# Qs:
# 1. What type of slots would go into 3. or how could I find this out by myself?
# 2. What type of data structure is 4. and what type of slots would go into it or how could I find this out by myself?

#### Add the ones from slots of ODD_spatial_points into ODD, and change DispX (take a copy of it into a new file) to have @data[["Population"]] 
# instead of $Population, etc. $ is not a feature of a general S4 object, but implemented by the Spatial class in the sp package.
# Then, look at including the data from the Philippines, how to do this. If the example runs ok after all this, look into the gridded background 
# data with aggregated impact data.