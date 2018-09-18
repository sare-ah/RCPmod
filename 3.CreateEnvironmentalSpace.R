#######################################################################################################################
# Regions of Common Profile - environmental space
# 
# Objective:  Convert and prepare rasters for prediction from RCP Model. This code should only need to be run once for the enviro vars of interest 
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

##### 
# Set up
##### 
# Set working directory - move back to parent directory
setwd('..')
dir<-getwd()

outdir <- file.path("./Models", Sys.Date() )  # Set this to Models/Tday's Date/...
suppressWarnings( dir.create( "../",outdir, recursive = TRUE ) )

#Location of rasters
rasLoc<-c("D:/Documents/!GIS Files/EnvironmentalRasters/Nearshore/Rasters/")

# Load helper functions
source("Scripts/RCP_Helper_Functions.R")

# Make packages available
UsePackages( pkgs=c("dplyr","RCPmod", "raster", "rasterVis", "tidyr","corrplot") ) 

##### 
# Transform input (site x species+env) data first
####
# Load species data to get only the environmental variables of interest 
spEnv<- read.csv("./Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv", header=T, stringsAsFactors=F)
spEnv<-spEnv[spEnv$area=="NCC",] #confirm NCC only
spEnv <- spEnv[!names(spEnv)%in% c( "TransDepth","fcode","BoP1","BoP2","area","surv","coords.x1","coords.x2","optional")]
species <- names(spEnv)[2:160] #159 species
enviro.variables<-names(spEnv)[162:ncol(spEnv)] #11 enviro
enviro.variables<-c("bathy","fetchSum","FBPI")

id_vars="SourceKey" # Unique identifier for each record
sample_vars="Survey" # Field or fields that describe potential sampling effects

# Transform environmental variables in input data into orthogonal polynomial terms
spEnv_poly_list <- poly_data(poly_vars=enviro.variables, degree=c(rep(1, length(enviro.variables))), 
                      id_vars=id_vars,sample_vars=sample_vars, 
                      species_vars=species, data=spEnv)
spEnv_poly_data <- spEnv_poly_list$rcp_data

saveRDS(spEnv_poly_data, file=paste0( outdir, "/spEnv_poly_data.rds"))
saveRDS(spEnv_poly_list, file=paste0( outdir, "/spEnv_poly_list.rds"))

##### 
# Transform environmental space (raster stack) 
####
# Load rasters in as raster stack then go back to project working directory
setwd(rasLoc)
ras<-stack(paste0(enviro.variables,".tif"))
setwd(dir)

# Convert rasters to dataframe
envRasterPoints <- as.data.frame(rasterToPoints(ras))
envRasterPoints <- na.omit( envRasterPoints )
saveRDS(envRasterPoints, file=paste0(outdir, "/envRasterPoints.rds"))

# Transform using stored polys, predict and plot results
envRasterPoints_poly <- poly_pred_space(envRasterPoints, spEnv_poly_list$poly_output,sampling_vals=NULL,sampling_name=NULL, sampling_factor_levels =NULL)
saveRDS(envRasterPoints_poly, file=paste0(outdir, "/envRasterPoints_poly.rds"))

#Transform environmental space again, but with sampling levels
#Right now this makes outputs IDENTICAL to the poly_pred_space, and identical for each factor level, differing only by an  an extra column with the survey name. Not sure what the point is.

sampling_factor_levels<-unique(spEnv$Survey)

#for (i in 1:length(sampling_factor_levels)){
envRasterPoints_poly_level <- poly_pred_space(envRasterPoints, spEnv_poly_list$poly_output,sampling_vals=sampling_factor_levels[1],
                                 sampling_name="Survey", sampling_factor_levels =sampling_factor_levels)

saveRDS(envRasterPoints_poly_level, file=paste0(outdir, "/envRasterPoints_poly_level_",sampling_factor_levels[1],".rds" ))
#}


