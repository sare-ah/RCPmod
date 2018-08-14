#######################################################################################################################
# Regions of Common Profile - environmental space
# 
# Objective:  Convert and prepare rasters for prediction from RCP Model. This code should only need to be run once for the enviro vars of interest 
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

#set up
# Set working directory - move back to parent directory
setwd('..')

outdir <- file.path("../Models", Sys.Date() )  # Set this to Models/Tday's Date/...
suppressWarnings( dir.create( "../",outdir, recursive = TRUE ) )

# Load helper functions
# =====================
# contains additional functions for installing packages, data transformations, plotting etc.
source("../Scripts/RCP_Helper_Functions.R")

# Make packages available
UsePackages( pkgs=c("dplyr","RCPmod", "raster", "rasterVis", "tidyr","corrplot") ) 

# Load species data to get only the environmental variables of interest 
spEnv<- read.csv("../Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv", header=T, stringsAsFactors=F)
spEnv<-spEnv[spEnv$area=="NCC",] #confirm NCC only
covariates.species <- spEnv[!names(spEnv)%in% c( "TransDepth","fcode","BoP1","BoP2","area","surv","coords.x1","coords.x2","optional")]
species <- names(covariates.species)[2:160] #159 species
enviro.variables<-names(covariates.species)[162:170] #9 enviro

id_vars="SourceKey" # Unique identifier for each record
sample_vars="Survey" # Field or fields that describe potential sampling effects

# Load rasters & create dataframe of prediction space --- NEED to set environmental layers up to be read in as a raster brick (?I think?)
setwd("D:/Documents/!GIS Files/EnvironmentalRasters/Nearshore/Rasters/")
rasSub<-stack(paste0(enviro.variables,".tif"))

# Set working directory - move back to parent directory
setwd('..')

# Convert rasters to dataframe and log transform depth
pred_space_rcp <- as.data.frame(rasterToPoints(rasSub))
pred_space_rcp <- na.omit( pred_space_rcp )
saveRDS(pred_space_rcp, file=paste0( outdir, "/pred_space_rcp.rds"))
pred_space_rcp<-readRDS(paste0(outdir, "/pred_space_rcp.rds"))

# Generate datafile with orthogonal polynomial terms

# ??poly_data() --- creates a list of two objects
# 1. Generate orthogonal polynomials for RCPmod input
# 2. Save polynomial bases to transform prediction space
# poly_data(predictor vars, poly degrees, vars not transformed, sampling level, species, data)
rcp_poly <- poly_data(poly_vars=enviro.variables, degree=c(rep(2, length(enviro.variables))), 
                      id_vars=id_vars,sample_vars=sample_vars, 
                      species_vars=species, data=covariates.species)

rcp_data <- rcp_poly$rcp_data

saveRDS(rcp_poly, file=paste0( outdir, "/rcp_poly.rds"))
saveRDS(rcp_data, file=paste0( outdir, "/rcp_data.rds"))


# ??poly_pred_space()  
# Transform using stored polys, predict and plot results
# Note: only accomodates one sampling term at the moment
# poly_pred_space(enviro.pts, orthogonal polynomials, sampling_vals, sampling_name, sampling_factor_levels)
rcp_poly_pred <- poly_pred_space(pred_space_rcp, rcp_poly$poly_output, 
                                 sampling_vals=NULL, 
                                 sampling_name=NULL, sampling_factor_levels =NULL)

saveRDS(rcp_poly_pred, file=paste0(outdir, "/rcp_poly_pred.rds"))

#Transform environmental space again, but with sampling levels
#Right now this makes outputs IDENTICAL to the poly_pred_space, and identical for each factor level, differing only by an  an extra column with the survey name.
#What is the use of this??
sampling_factor_levels<-unique(covariates.species$Survey)

for (i in 1:length(sampling_factor_levels)){
rcp_poly_pred <- poly_pred_space(pred_space_rcp, rcp_poly$poly_output, 
                                 sampling_vals=sampling_factor_levels[i], 
                                 sampling_name="Survey", sampling_factor_levels =sampling_factor_levels)

saveRDS(rcp_poly_pred, file=paste0(outdir, "/rcp_poly_pred_level_",sampling_factor_levels[i],".rds" ))
}


