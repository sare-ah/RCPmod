#######################################################################################################################
# Regions of Common Profile
# 
# Objective:  Determine the regions of common profile from benthic habitat mapping dive survey data
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

# start fresh
rm(list=ls())

# Check which version of R is being used and reset if necessary
Sys.getenv("R_ARCH")
# The message returned will tell you which version of R is being used
# "/i386" 32 bit R --- which is necessary to grab data from MS Access database
# "/64"   64 bit R
# To reset: Select Tools menu | Global Options... | R Version: | Change
# Then you will have to open and close R for the changes to take effect

# Set working directory - move back to parent directory
setwd('..')

# Load helper functions
# =====================
# contains additional functions for installing packages, data transformations, plotting etc.
source("./Scripts/RCP_Helper_Functions.R")

# Make packages available
UsePackages( pkgs=c("dplyr","RCPmod", "raster", "rasterVis", "tidyr", "png") ) 

# Controls
# ========
enviro.variables <- c("Long_MP", "log_depth", "caisom_floor_temperature" )  # CHANGE to match environmental variables used for BHM
id_vars="HaulIndex"
sample_vars="Season"

distribution <- "NegBin"  # Change to "Bernoulli" for P/A data
gen.start.val <- "random2"
n.sites <- 181 # Number of sites in study area

# Load the data
# =============
# Fish example code has species and covariates in the same csv file
covariates.species <- read.csv("./Sample_data/SubAntFish_bioenv.csv", header=T, stringsAsFactors=F)
species <- names(covariates.species)[9:23]

# Generate datafile with orthogonal polynomial terms
rcp_env_vars <- c( "Long_MP", "log_depth", "caisom_floor_temperature" )  # CHANGE to match environmental variables used for BHM

# ??poly_data() --- creates a list of two objects
# 1. Generate orthogonal polynomials for RCPmod input
# 2. Save polynomial bases to transform prediction space
# poly_data(predictor vars, poly degrees, vars not transformed, sampling level, species, data)
rcp_poly <- poly_data(poly_vars=rcp_env_vars, degree=c(2,2,2), 
                    id_vars="HaulIndex",sample_vars="Season", 
                    species_vars=species, data=covariates.species)
rcp_data <- rcp_poly$rcp_data

# Load rasters & create dataframe of prediction space --- NEED to set environmental layers up to be read in as a raster brick (?I think?)
pred_masked <- brick( "./Sample_data/pred_masked" )
plot(pred_masked)

# Convert rasters to dataframe and log transform depth
pred_space_rcp <- as.data.frame(rasterToPoints(
  subset(pred_masked, c( "Long_MP", "bathymetry", "caisom_floor_temperature" ))))
pred_space_rcp <- na.omit( pred_space_rcp )
pred_space_rcp$log_depth <- log( pred_space_rcp$bathymetry* -1 )

# ??poly_pred_space()  --- I AM CONFUSED BY THIS FUNCTION
# Transform using stored polys, predict and plot results
# Note: only accomodates one sampling term at the moment
# poly_pred_space(enviro.pts, orthogonal polynomials, sampling_vals, sampling_nam, sampling_factor_levels)
rcp_poly_pred <- poly_pred_space(pred_space_rcp, rcp_poly$poly_output, 
                               sampling_vals="Autumn/Winter", 
                               sampling_name="Season", sampling_factor_levels = c("Autumn/Winter","Spring","summer"))

# Create RCP formula
form <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",paste(names(rcp_data)[18:23], collapse="+")))
form

# Run RCPs --- with sampling effect
# =================================
# Season/Year as sampling effect ---- CHANGE to survey (month/year) for BHM data
# This bit  will take a while on a single core machine
nstarts <- 100  # editted this line from nstarts <- 1000
max.nRCP <- 6   # editted this line from max.nRCP <- 6
nRCPs_samp <- list() 

# for( ii in 1:max.nRCP) # This for loop throws an error when nRCP equals 2. Function is aborted. Confusing b/c the example starts at nRCP=1
#    nRCPs_samp[[ii]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=ii,
#                                         inits="random2", nstart=nstarts, dist="NegBin", mc.cores=1)

# Populate list manually
nRCPs_samp[[1]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=1,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[2]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=2,
                                    inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1) # throws an error!
nRCPs_samp[[3]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=3,
                                    inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[4]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=4,
                                    inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[5]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=5,
                                    inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[6]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=6,
                                    inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)

saveRDS(nRCPs_samp, file="nRCPs_samp.rds")

# Determine best model
# ====================
# Name list elements & remove elements with null values (model failed to run properly for nRCP=2)
names(nRCPs_samp) <- seq_along(nRCPs_samp)
nRCPs_samp[sapply(nRCPs_samp, is.null)] <- NULL

# Get BICs
RCPsamp_BICs <- sapply( nRCPs_samp, function(x) sapply( x, function(y) y$BIC))

# Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
RCPsamp_minPosteriorSites <- cbind( n.sites, sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
RCPsamp_ObviouslyBad <- RCPsamp_minPosteriorSites < 2
RCPsamp_BICs[RCPsamp_ObviouslyBad] <- NA

# Plot minimum BIC for each nRCP
RCPsamp_minBICs <- apply(RCPsamp_BICs, 2, min, na.rm=TRUE)

plot( names(nRCPs_samp), RCPsamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep(names(nRCPs_samp), each=nrow( RCPsamp_BICs)), RCPsamp_BICs, pch=20)

# Choose best model from above 
# set variables to match the lowest BIC with a number of RCPs
nRCP_best <- as.numeric(which.min(RCPsamp_minBICs))
inits_best <- unlist( nRCPs_samp[[nRCP_best]][[3]]$coef)

# Rerun best model (for full output)
# ==================================
control <- list( optimise=FALSE, quiet=FALSE)
RCPsamp_fin <- regimix(form.RCP=form, form.spp=~Season, 
                     nRCP=nRCP_best, data=rcp_data, dist="NegBin", inits = inits_best, control=control)

# Clean-up workspace
rm(RCPsamp_BICs,RCPsamp_minPosteriorSites, RCPsamp_ObviouslyBad, RCPsamp_minBICs, RCPsamp_goodun, control)

# Plot model diagnostics
# ======================
# To do:  Add code to save figures to pdf files
# Residual plots
plot.regimix(RCPsamp_fin, type="RQR", fitted.scale="log")

# Cooks Distance Plots
# Takes a while...
tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,3,4,5,6,7,8,9,10,20,30,40,50), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # original
#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # simplified, but missing horizontal line
plot( tmp, minWidth=2, ncuts=111 ) 

saveRDS(tmp, file="Fish.Cooks.Dist.rds")
tmp <- readRDS("Fish.Cooks.Dist.rds")

# Examine dispersion parameter for negative Binomial  --- would this change with a Bernoulli distribution?
par(mfrow=c(1,1))
hist(RCPsamp_fin$coefs$disp, xlab="Dispersion Parameter", 
     main="Negative Binomial Model", col="grey", cex.main=0.8, cex=0.8, cex.lab=0.8 )

# Generate bootstrap estimates of parameters
# ==========================================
# Warning --- bootstrap is slooow
# rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=1000, mc.cores=1) # original code
rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=100, mc.cores=1) # simplified

# Evaluate sampling factor effect
# ===============================
# ??Sp_abund_all()
# Calculate average, SD and CI of species abundances in each  RCP
# Sp_abund_all(bootstrap object) --- FUNCTION IS FOR FISH EXAMPLE --- NEED TO REWRITE!!!
RCP_abund_samp<-Sp_abund_all(rcpsamp_boots)

# Get autumn values and make pretty
aut_samp <- as.data.frame(matrix(data=paste0(sprintf("%.2f", round(RCP_abund_samp$autumn$mean,2)), " (", sprintf("%.2f",round(RCP_abund_samp$autumn$sd,2)), ")"),
                               ncol=3, nrow=15, byrow=TRUE))
names(aut_samp) <- paste0("RCP", 1:3)
rownames(aut_samp) <- gsub("."," ", species, fixed=TRUE)

# ??sampling_dotplot2()
# Plot of sampling factor effects
sampling_dotplot2(RCPsamp_fin,rcpsamp_boots,legend_fact=c("Spring", "Summer"), col=c("black", "red"), lty=c(1,2))

# Spatial Predictions
# ===================
RCPsamp_SpPreds<-predict.regimix(object=RCPsamp_fin, object2=rcpsamp_boots, newdata=rcp_poly_pred)
# ??predict_maps2_SDF2()
# Plot RCP predictions
predict_maps2_SDF2(RCPsamp_SpPreds, pred_space=pred_space_rcp, pred_crop=pred_masked, nRCP=3)

# To do:
# Is there a way to export predictions to a shapefile?

