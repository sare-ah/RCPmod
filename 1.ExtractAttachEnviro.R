#######################################################################################################################
# Extract and attach enviromental data
# 
# Objective:  Get environmental variables at sampling points and attach to dataframe of species observations
#
# Author:     Katie Gale
#             Katie.Gale@dfo-mpo.gc.ca
# Date:       August 7, 2018
######################################################################################################################

library(rgdal)
library(raster)

################################
# Import rasters for analysis  #
################################
setwd("D:/Documents/!GIS Files/EnvironmentalRasters/Nearshore/Rasters/")
rasfiles<-list.files(getwd(), pattern = "(*.)tif$",recursive=T)
ras<-stack(rasfiles[c(1:3, 5:26)])

##################################################
# Bring in spatialized points with species data 
##################################################
setwd("D:/Documents/Projects/WorldClass_DiveSurveys/RCP/")
sp<-readOGR("./Data/SpatialPointsWithSpecies/SpatialPointsWithSpecies.shp")
sp<-spTransform(sp, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"))
proj4string(ras)==proj4string(sp)

########################################################
# Extract enviromental data from rasters to points
########################################################
spEnv = data.frame(sp,raster::extract(ras, sp)) ##add .parameters to each record
summary(is.na(spEnv)) 

spEnv_complete<-spEnv[complete.cases(spEnv[!names(spEnv) %in% names(sp)]),] #remove cases where environmental data aren't available
summary(is.na(spEnv_complete)) 
write.csv(spEnv_complete,"./Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviro.csv", row.names=F)
