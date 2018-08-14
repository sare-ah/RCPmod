#######################################################################################################################
# Regions of Common Profile
# 
# Objective:  Determine the regions of common profile from benthic habitat mapping dive survey data
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

# TO DO:
# Add code to export spatial predictions as a raster layer
# Add code to build species profiles of each RCP (probability of species being a member of an RCP)
# Write documentation for this code

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
source("../Scripts/RCP_Helper_Functions.R")

# Make packages available
UsePackages( pkgs=c("dplyr","RCPmod", "raster", "rasterVis", "tidyr","corrplot") ) 

# Controls
# ========
# REQUIRED
id_vars="SourceKey" # Unique identifier for each record
sample_vars="Survey" # Field or fields that describe potential sampling effects
nstarts <- 10  # Editted this line from nstarts <- 1000
#max.nRCP <- 6   # Not used when manually testing number of RCPs
distribution <- "Bernoulli"  # Change to "Bernoulli" for P/A data or NegBin for abund
gen.start.val <- "random2"

# OPTIONAL
min.prevalence <- TRUE # specify T/F to subset data based on species prevalence
species.n <- 30  # Minimum species prevalence
subset.data <- T # Specify T/F to subset data or use all communities 
subset.size <- 200 # Specifcy random subset size --- originally set to 1000


# Directories  
# ===========
# Subfolder for output 
outdir <- file.path("../Models", Sys.Date() )  # Set this to Models/Tday's Date/...
suppressWarnings( dir.create( outdir, recursive = TRUE ) )

# Load the data
# =============
# Fish example code has species and covariates in the same csv file
spEnv<- read.csv("../Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv", header=T, stringsAsFactors=F)
spEnv<-spEnv[spEnv$area=="NCC",] #confirm NCC only
covariates.species <- spEnv[!names(spEnv)%in% c( "TransDepth","fcode","BoP1","BoP2","area","surv","coords.x1","coords.x2","optional")]

siteNo <- covariates.species[,id_vars] 
species <- names(covariates.species)[2:160] #159 species
enviro.variables<-names(covariates.species)[162:170] #9 enviro

n.sites <- nrow(covariates.species) # Number of sites in study area
#n.abund <- 2 # Define where species data starts in covariates.species, test with names(covariates.species)[n.abund] ??

rcp_data <- readRDS(paste0(outdir, "/rcp_data.rds"))

# Conditional processing to subset sites used - ensure sitename column is correctly specified in controls
if (subset.data) {
  # specify the  subset to use in the analysis
  set.seed(subset.size)
  sample.sites <- sample(siteNo, subset.size)
  rcp_data <- rcp_data[rcp_data[,id_vars] %in% sample.sites,]
  print(paste0("Successfully subsetted [",subset.size,"] random sites"))
} else {
  print("No subsetting performed")
}

# Record site order --- we may need this later if we are subsetting the data 
site.names <- rcp_data[,id_vars]

# Conditional processing to subset number of species used based on species prevalence
if (min.prevalence) {
  # Determine in decreasing order the total count of each species within the study area
  species.count <- data.frame(count=sort(colSums(rcp_data[,names(rcp_data) %in% species]), decreasing=T))
  species.count$species <- row.names(species.count)
  # Select the species to be dropped and not used in the model
  dropped.species <- species.count$species[species.count$count < species.n]
  # Remove species columns from covariates.species that are not prevalent enough - otherwise poly_data() fails!
  rcp_data <- rcp_data[ , !(names(rcp_data) %in% dropped.species)]
  # Remove species names from species vector - otherwise poly_data() fails!
  species <- setdiff(species, dropped.species)
}

length(species) #38 species
length(enviro.variables) #9 enviro
length(site.names) #200

## In paper they use forward selection to select variables. Do we need to do that???


#filename <- paste( "..", outdir, "Enviro.Layers.pdf", sep = "/" )
# pdf( file=filename, width = 6, height = 4 )
#   plot(ras)
# dev.off()


# Create RCP formula
form <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",c(paste(c(paste0(enviro.variables, "1"), paste0(enviro.variables, "2")), collapse="+"))))


# Run RCPs --- WITH sampling effect
# =================================
# Season/Year as sampling effect ---- CHANGE to survey (month/year) for BHM data
# This bit  will take a while on a single core machine
nRCPs_samp <- list() 

# for( ii in 1:max.nRCP) # This for loop throws an error when nRCP equals 2. Function is aborted. Confusing b/c the example starts at nRCP=1
#    nRCPs_samp[[ii]] <- regimix.multifit(form.RCP=form, form.spp= ~ Season, data=rcp_data, nRCP=ii,
#                                         inits="random2", nstart=nstarts, dist="NegBin", mc.cores=1)

# Populate list manually
#nRCPs_samp[[1]] <- regimix.multifit(form.RCP=form, form.spp= ~ sample_vars, data=rcp_data, nRCP=1,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
#nRCPs_samp[[2]] <- regimix.multifit(form.RCP=form, form.spp= ~ sample_vars, data=rcp_data, nRCP=2,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1) # throws an error!
nRCPs_samp[[3]] <- regimix.multifit(form.RCP=form, form.spp= paste0("~ ", sample_vars), data=rcp_data, nRCP=3,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[4]] <- regimix.multifit(form.RCP=form, form.spp= paste0("~ ", sample_vars), data=rcp_data, nRCP=4,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[5]] <- regimix.multifit(form.RCP=form, form.spp= paste0("~ ", sample_vars), data=rcp_data, nRCP=5,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[6]] <- regimix.multifit(form.RCP=form, form.spp= paste0("~ ", sample_vars), data=rcp_data, nRCP=6,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_samp[[7]] <- regimix.multifit(form.RCP=form, form.spp= paste0("~ ", sample_vars), data=rcp_data, nRCP=7,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)


saveRDS(nRCPs_samp, file=paste0(outdir,"/nRCPs_samp.rds"))
nRCPs_samp <- readRDS(paste0(outdir,"/nRCPs_samp.rds"))

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

pdf( file= paste0(outdir, "/nRCPs_Samp.pdf"), width = 4, height = 4 )
plot( names(nRCPs_samp), RCPsamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep(names(nRCPs_samp), each=nrow( RCPsamp_BICs)), RCPsamp_BICs, pch=20)
dev.off()

# Choose best model from above with the lowest BIC
minBICs <- as.data.frame(RCPsamp_minBICs) # Create dataframe to manipulate b/c I can't capture index value of a list
minBICs$RCP <- as.numeric( row.names(minBICs) ) # Create column for number of RCPs
nRCP_best <- minBICs$RCP[which.min(minBICs$RCPsamp_minBICs)]

# set variables to match the number of RCPs with the lowest BIC 
#inits_best <- unlist( nRCPs_samp[[nRCP_best]][[3]]$coef) #didn't use this - ??. using this causes model to fail from here on. Flag

# Rerun best model (for full output)
# ==================================
control <- list( optimise=FALSE, quiet=FALSE)
RCPsamp_fin <- regimix(form.RCP=form, form.spp=paste0("~ ", sample_vars), 
                       data=rcp_data, dist=distribution, nRCP=nRCP_best, inits = "random2", control=control)

saveRDS(RCPsamp_fin, file=paste0(outdir,"/nRCPsamp_fin.rds"))
RCPsamp_fin <- readRDS(paste0(outdir,"/nRCPsamp_fin.rds"))

# Clean-up workspace
#rm(RCPsamp_BICs,RCPsamp_minPosteriorSites, RCPsamp_ObviouslyBad, RCPsamp_minBICs, control, minBICs)

# Plot model diagnostics
# ======================
# To do:  Add code to save figures to pdf files
# Residual plots
pdf( file=paste0(outdir, "/Residuals_Samp.pdf"), width = 6, height = 4 )
plot.regimix(RCPsamp_fin, type="RQR", fitted.scale="log")
dev.off()


# Cooks Distance Plots #????
# Takes a while...
#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=NULL, mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # original
tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=3, doPlot=FALSE) # quick run (3 runs instead of 200)

#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # simplified, but missing horizontal line
saveRDS( tmp, file=paste0(outdir,"/Fish.Cooks.Dist.rds" ))

#tmp <- readRDS("Fish.Cooks.Dist.rds")
pdf( file=paste0(outdir, "/Cooks.Distance.Samp.pdf"), width = 6, height = 4 )
plot( tmp, minWidth=2, ncuts=111 ) 
dev.off()

# Examine dispersion parameter for negative Binomial  --- would this change with a Bernoulli distribution?
# #No "disp" parameter in bernoulli output
# par( mfrow=c(1,1) )
# filename <- paste( "..", outdir, "Dispersion.Parameter.pdf", sep = "/" )
# pdf( file=filename, width = 4, height = 4)
# hist(RCPsamp_fin$coefs$disp, xlab="Dispersion Parameter", 
#      main="Negative Binomial Model", col="grey", cex.main=0.8, cex=0.8, cex.lab=0.8 )
# dev.off()

# Generate bootstrap estimates of parameters
# ==========================================
# Warning --- bootstrap is slooow
# rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=1000, mc.cores=1) # original code
rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=10, mc.cores=1) # simplified
saveRDS( rcpsamp_boots, file=paste0(outdir,"/RCPSamp_boots.rds" ))
#rcpsamp_boots <- readRDS(paste0(outdir,"/RCPSamp_boots.rds"))


# Evaluate sampling factor effect
# ===============================
# Calculate average, SD and CI of species abundances in each  RCP
# Sp_abund_gen(bootstrap object, input species file with levels)
RCP_abund_samp<-Sp_abund_gen(rcpsamp_boots, covariates.species)

#Plot results
png(paste0(outdir,"/RCP_abund_sp.png"), height=20, width=20, res=300, units="cm", pointsize = 12)
sp_abund_plot(RCP_abund_samp,"Survey")
dev.off()

# ??sampling_dotplot2()
# Plot of sampling factor effects (relative to base case)
pdf( file=paste0(outdir, "/Sampling.Factor.Effects.pdf"), width = 8, height = 8)
sampling_dotplot2(RCPsamp_fin,rcpsamp_boots,legend_fact=c("PAC 2014-29", "PAC 2014-58","PAC 2015-52"), col=c("black", "red", "blue"), lty=c(1,2,3))
dev.off()


# Spatial Predictions
# ===================

#Load environmental space
#rcp_poly_pred<-readRDS(paste0(outdir,"/rcp_poly_pred.rds"))
rcp_poly_pred_lev<-readRDS(paste0(outdir,"/rcp_poly_pred_level_P1458.rds"))
rcpsamp_boots<-readRDS(paste0(outdir,"/rcpsamp_boots.rds"))
nRCPsamp_fin<-readRDS(paste0(outdir,"/nRCPsamp_fin.rds"))


RCPsamp_SpPreds<-predict.regimix(object=nRCPsamp_fin, object2=rcpsamp_boots, newdata=rcp_poly_pred_lev)


# ??predict_maps2_SDF2()
# Plot RCP predictions
filename <- paste( ".", outdir, "Spatial.Preds.Samp.pdf", sep = "/" )
pdf( file=filename, width = 6, height = 4 )
predict_maps2_SDF2(RCPsamp_SpPreds, pred_space=pred_space_rcp, pred_crop=pred_masked, nRCP=nRCP_best)
dev.off()

# Clean-up variables...
rm(nRCP_best)

#================================================================================================
#################################################################################################
#================================================================================================

# Run RCPs --- WITHOUT sampling effect
# ====================================
nRCPs_NoSamp <- list()

# for( ii in 1:max.nRCP){
#   nRCPs_NoSamp[[ii]] = tryCatch({regimix.multifit(form.RCP=form, data=rcp_data, nRCP=ii,
#                                                   inits="random2", nstart=nstarts, dist="NegBin", mc.cores=1)},
#                                 error=function(ii){return(NA)})
# }

# Populate list manually
nRCPs_NoSamp[[1]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=1,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_NoSamp[[2]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=2,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1) # throws an error!
nRCPs_NoSamp[[3]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=3,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_NoSamp[[4]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=4,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_NoSamp[[5]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=5,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
nRCPs_NoSamp[[6]] <- regimix.multifit(form.RCP=form, data=rcp_data, nRCP=6,
                                      inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)

saveRDS(nRCPs_NoSamp, file="nRCPs_NoSamp.rds")

# Determine best model
# ====================
# Name list elements & remove elements with null values (model failed to run properly for nRCP=2)
names(nRCPs_NoSamp) <- seq_along(nRCPs_NoSamp)
nRCPs_NoSamp[sapply(nRCPs_NoSamp, is.null)] <- NULL

# Get BICs
RCPNoSamp_BICs <- sapply( nRCPs_NoSamp, function(x) sapply( x, function(y) y$BIC))

# Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
RCPNoSamp_minPosteriorSites <- cbind( n.sites, sapply( nRCPs_NoSamp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
RCPNoSamp_ObviouslyBad <- RCPNoSamp_minPosteriorSites < 2
RCPNoSamp_BICs[RCPNoSamp_ObviouslyBad] <- NA

# Plot minimum BIC for each nRCP
RCPNoSamp_minBICs <- apply( RCPNoSamp_BICs, 2, min, na.rm=TRUE )
filename <- paste( ".", outdir, "nRCPs_NoSamp.pdf", sep = "/" )
pdf(file=filename, width = 6, height = 4)
plot( names(nRCPs_NoSamp), RCPNoSamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep(names(nRCPs_NoSamp), each=nrow( RCPNoSamp_BICs)), RCPNoSamp_BICs, pch=20)
dev.off()

# Choose best model from above with the lowest BIC
minBICs <- as.data.frame(RCPNoSamp_minBICs) # Create dataframe to manipulate b/c I can't capture index value of a list
minBICs$RCP <- as.numeric( row.names(minBICs) ) # Create column for number of RCPs
nRCP_best <- minBICs$RCP[which.min(minBICs$RCPNoSamp_minBICs)]

# set variables to match the number of RCPs with the lowest BIC 
inits_best <- unlist( nRCPs_NoSamp[[nRCP_best]][[3]]$coef)

# Rerun best model (for full output)
# ==================================
control <- list( optimise=FALSE, quiet=FALSE)
RCPNoSamp_fin<-regimix(form.RCP=form, 
                       nRCP=nRCP_best, data=rcp_data, dist=distribution, inits = inits_best, control=control)

# Clean-up workspace
rm(RCPNoSamp_BICs,RCPNoSamp_minPosteriorSites, RCPNoSamp_ObviouslyBad, RCPNoSamp_minBICs, minBICs)

# Plot model diagnostics
# ======================
# Residual plot
filename <- paste( ".", outdir, "Residuals_NoSamp.pdf", sep = "/" )
pdf(file=filename, width = 6, height = 4)
plot.regimix(RCPNoSamp_fin, type="RQR", fitted.scale="log") 
dev.off()

# Cooks Distance Plots
# will take a while to run
tmp <- stability.regimix(RCPNoSamp_fin, oosSizeRange=c(1,2,3,4,5,6,7,8,9,10,20,30,40,50), mc.cores=1, times=RCPsamp_fin$n)
#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # simplified, but missing horizontal line
saveRDS(tmp, file="Fish.Cooks.Dist.NoSamp.rds")
#tmp <- readRDS("Fish.Cooks.Dist.NoSamp.rds")
filename <- paste( ".", outdir, "Cooks.Distance.NoSamp.pdf", sep = "/" )
pdf(file=filename, width = 6, height = 4)
plot( tmp, minWidth=2, ncuts=111 ) 
dev.off()

# Generate bootstrap estimates of parameters
# ==========================================
# Warning --- bootstrap is slooow
rcpNoSamp_boots <- regiboot( RCPNoSamp_fin, type="BayesBoot", nboot=1000, mc.cores=1 )
saveRDS( rcpNoSamp_boots, file="RCPNoSamp_boots.rds" )
#rcpNoSamp_boots <- readRDS("RCPNoSamp_boots.rds")


# Spatial Predictions
# ===================
RCPNoSamp_SpPreds <- predict.regimix( object=RCPNoSamp_fin, object2=rcpNoSamp_boots, newdata=rcp_poly_pred )
# ??predict_maps2_SDF2()
# Plot RCP predictions
# function(predictions,     #output from predict.regimix
# pred_space,      #dataframe containing coordinates for the prediction space
# pred_crop,       #raster of extent of prediction space (used in rasterize function)
# nRCP,            #the number of RCPs
# my.ylim=NULL,    #y limit to plot
# my.xlim=NULL,    #x limit to plot
# my.asp=1)        #aspect for plotting

filename <- paste( ".", outdir, "Spatial.Preds.NoSamp.pdf", sep = "/" )
pdf( file=filename, width = 4, height = 4 )
predict_maps2_SDF2( RCPNoSamp_SpPreds, pred_space=pred_space_rcp, pred_crop=pred_masked, nRCP=nRCP_best )
dev.off()

