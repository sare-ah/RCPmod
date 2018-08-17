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
# Add code to export spatial predictions as a raster layer (KG: rasterize function as Hill wrote fails. I was able to output it as a spatialpointsdataframe, which you can view in R.)
# Add code to build species profiles of each RCP (probability of species being a member of an RCP) (KG: Done, assuming you can interpret the mean "prevalence" (alpha) of the each species as its probability)
# Write documentation for this code

#### Outstanding questions
##########################
# Should we need to set up variable selection for the model? In the paper they do forward variable selection wrapped with selection of the number of RCPs. I've started a script to do this, and Jess suggested running it in parallel to reduce run time. Still need to be thoughtful about iterations. 18 enviro vars (9 x 2 polynomials) is 262143 models. Not doing the polynomials is only 511. 
# Do we actually need to calculate orthogonal polynomials? Hill says it's to "avoid convergence problems". Not sure what that means.  
# I don't think the model outputs tells us "the most important environmental variables for distinguishing clusters". Doing a BIC comparison for models with/without each variable might give this. 
# What is actually being bootstrapped in regiboot? Is it model inputs bootstrapped and run through the single model, or inputs bootstrapped to get a new model? 
# I don't think the results tell us "effect of Survey". So we could add this to the model selection/BIC list to see how much of an impact including Survey has on the model BIC?


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

# Directories  
# ===========
# Subfolder for output 
outdir <- file.path("../Models", Sys.Date() )  # Set this to Models/Tday's Date/...
suppressWarnings( dir.create( outdir, recursive = TRUE ) )

setwd(outdir)

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


# Load the data
# =============
# Fish example code has species and covariates in the same csv file
spEnv<- read.csv("../Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv", header=T, stringsAsFactors=F)
spEnv<-spEnv[spEnv$area=="NCC",] #confirm NCC only
covariates.species <- spEnv[!names(spEnv)%in% c( "TransDepth","fcode","BoP1","BoP2","area","surv","coords.x1","coords.x2","optional")]

# Load some outputs of 3.create environmental space
rcp_data <- readRDS("rcp_data.rds") #data.frame of species data with polynomial-transformed environmental variables
rasSub<- readRDS("rasSub.rds")      #raster stack of predictor variables

siteNo <- covariates.species[,id_vars] 
species <- names(covariates.species)[2:160] #159 species
enviro.variables<-names(covariates.species)[162:170] #9 enviro

n.sites <- nrow(covariates.species) # Number of sites in study area

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

saveRDS(rcp_data, "rcp_data_subsetted.rds")


## In paper they use forward selection to select variables. Do we need to do that???

# pdf( file=paste0( outdir, "/Enviro.Layers.pdf"), width = 6, height = 4 )
#   plot(ras)
# dev.off()


# Create RCP formula
form <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",c(paste(c(paste0(enviro.variables, "1"), paste0(enviro.variables, "2")), collapse="+"))))


# Run RCPs --- WITH sampling effect
# =================================
# Using "Survey" as sampling effect (defined above)
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

saveRDS(nRCPs_samp, file="nRCPs_samp.rds")
nRCPs_samp <- readRDS("nRCPs_samp.rds")

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

pdf( file= "nRCPs_Samp.pdf", width = 4, height = 4 )
plot( names(nRCPs_samp), RCPsamp_minBICs, type='b', ylab="BIC", xlab="nRCP", pch=20)
points( rep(names(nRCPs_samp), each=nrow( RCPsamp_BICs)), RCPsamp_BICs, pch=20)
dev.off()

# Choose best model from above with the lowest BIC
minBICs <- as.data.frame(RCPsamp_minBICs) # Create dataframe to manipulate b/c I can't capture index value of a list
minBICs$RCP <- as.numeric( row.names(minBICs) ) # Create column for number of RCPs
nRCP_best <- minBICs$RCP[which.min(minBICs$RCPsamp_minBICs)]

# set variables to match the number of RCPs with the lowest BIC 
#inits_best <- unlist( nRCPs_samp[[nRCP_best]][[3]]$coef) #FLAG - using these defined coefficients causes the model to fail running. Changed inits to "random2" instead but might need to be fixed.

# Rerun best model (for full output)
# ==================================
control <- list( optimise=FALSE, quiet=FALSE)
RCPsamp_fin <- regimix(form.RCP=form, form.spp=paste0("~ ", sample_vars), 
                       data=rcp_data, dist=distribution, nRCP=nRCP_best, inits = "random2", control=control)

saveRDS(RCPsamp_fin, file="nRCPsamp_fin.rds")
RCPsamp_fin <- readRDS("nRCPsamp_fin.rds")

# Clean-up workspace
#rm(RCPsamp_BICs,RCPsamp_minPosteriorSites, RCPsamp_ObviouslyBad, RCPsamp_minBICs, control, minBICs)

# Plot model diagnostics
# ======================
# To do:  Add code to save figures to pdf files
# Residual plots
pdf( file=paste0("Residuals_Samp.pdf"), width = 6, height = 4 )
plot.regimix(RCPsamp_fin, type="RQR", fitted.scale="log")
dev.off()


# Cooks Distance Plots #????
# Takes a while...
#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=NULL, mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # original
tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=3, doPlot=FALSE) # quick run (3 runs instead of 200)

#tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # simplified, but missing horizontal line
saveRDS( tmp, file="Fish.Cooks.Dist.rds" )

#tmp <- readRDS("Fish.Cooks.Dist.rds")
pdf( file="Cooks.Distance.Samp.pdf", width = 6, height = 4 )
plot( tmp, minWidth=2, ncuts=111 ) 
dev.off()

# Examine dispersion parameter for negative Binomial  --- would this change with a Bernoulli distribution?
# FLAG: doesn't work because there's no "disp" parameter in bernoulli output
# par( mfrow=c(1,1) )
# pdf( file="Dispersion.Parameter.pdf", width = 4, height = 4)
# hist(RCPsamp_fin$coefs$disp, xlab="Dispersion Parameter", 
#      main="Negative Binomial Model", col="grey", cex.main=0.8, cex=0.8, cex.lab=0.8 )
# dev.off()

# Generate bootstrap estimates of parameters
# ==========================================
# Warning --- bootstrap is slooow
# rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=1000, mc.cores=1) # original code
rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=10, mc.cores=1) # simplified
saveRDS( rcpsamp_boots, file="RCPSamp_boots.rds" )
rcpsamp_boots <- readRDS("RCPSamp_boots.rds")


# Evaluate sampling factor effect
# ===============================
# Calculate average, SD and CI of species abundances in each RCP
# Sp_abund_gen(bootstrap object, input species file with levels)
RCP_abund_samp<-Sp_abund_gen(rcpsamp_boots, covariates.species)

#Plot results
png("RCP_abund_sp.png", height=20, width=20, res=300, units="cm", pointsize = 12)
sp_abund_plot(RCP_abund_samp,"Survey")
dev.off()

# ??sampling_dotplot2()
# Plot of sampling factor effects (relative to base case)
# Flag - could be generalized to make sure the factor levels are right. Not sure here because the "base case" isn't listed. 
pdf( file="Sampling.Factor.Effects.pdf", width = 8, height = 8)
sampling_dotplot2(RCPsamp_fin,rcpsamp_boots,legend_fact=c("PAC 2014-29", "PAC 2014-58","PAC 2015-52"), col=c("black", "red", "blue"), lty=c(1,2,3))
dev.off()


#What do the beta values coming out of regiboot object tells us? "RCPs dependence on covariates". Is looking at them useful? 


# Spatial Predictions
# ===================

#Load environmental space
#rcp_poly_pred<-readRDS(paste0(outdir,"/rcp_poly_pred.rds")
rcp_poly_pred_lev<-readRDS("rcp_poly_pred_level_P1362.rds")
rcpsamp_boots<-readRDS("rcpsamp_boots.rds")
nRCPsamp_fin<-readRDS("nRCPsamp_fin.rds")
pred_space_rcp<-readRDS("pred_space_rcp.rds")

#This command maxes out my ram, but works on GIS computer. Looks like Fig 4 in hill paper. 
#This is a probability of prediction for every raster cell in the rasSub stack. 
#Predict.regimix takes "newdata" of env+sp variables and calcualates probability of RCP membership. 
#it does this for the single model output (object) and the boot output (object2)
RCPsamp_SpPreds<-predict.regimix(object=nRCPsamp_fin, object2=rcpsamp_boots, newdata=rcp_poly_pred_lev)

RCPsamp_SpPreds<-readRDS("RCPsamp_SpPreds.rds")

df<-SpatialPointsDataFrame(coords= pred_space_rcp[,1:2], data=as.data.frame(RCPsamp_SpPreds$bootPreds))
rgdal::writeOGR(df, "/spatial", layer="bootPreds", driver="ESRI Shapefile")

dbf<-foreign::read.dbf("spatial/bootPreds.dbf")
dbf$probOfAssign<-apply(dbf, 1, FUN=max)
dbf$topRCP<-apply(dbf,1,function(x) match(max(x),x)) 
foreign::write.dbf(dbf,"spatial/bootPreds.dbf")

# ??predict_maps2_SDF2()
# Plot RCP predictions
nRCP_best=3

pdf( file=paste0("Spatial.Preds.Samp.pdf" ), width = 6, height = 4 )
predict_maps2_SDF2(RCPsamp_SpPreds, pred_space=pred_space_rcp, pred_crop=rasSub, nRCP=nRCP_best)
dev.off()

## giving predict.regimix the input sampling sites (here, rcp_data) gives probability of RCP membership for each site. 
rcp_data<-readRDS("rcp_data.rds")
RCPsamp_SpPreds_intialSites<-predict.regimix(object=nRCPsamp_fin, object2=rcpsamp_boots, newdata=rcp_data)

#the outputs are both point predictions (from model?) and bootstrap predictions
#take the point predictions and plot them - it looks like Figure 3 in Hill paper
pt_preds<-as.data.frame(RCPsamp_SpPreds_intialSites[[1]])
rcp_data_melt<-melt(covariates.species[,!names(covariates.species)%in% species], id.vars=c("SourceKey","Survey"))
summary(rcp_data_melt$SourceKey==rcp_data$SourceKey) #order of points is the same. attach the spatial predictions
rcp_data_wPred<-cbind.data.frame(rcp_data_melt, pt_preds)
rcp_data_wPred_melt<-melt(rcp_data_wPred, id.vars=c("SourceKey","Survey", "variable","value"))
names(rcp_data_wPred_melt)[3:6]<-c("enviro","enviro_val","RCP","RCP_prob")

#Plot - matches Figure 3 in the paper. 
png("envirospace.png",res=400, width=50, height=20, units="cm")
ggplot(data=rcp_data_wPred_melt, aes(x=enviro_val, y=RCP_prob, col=RCP))+geom_point(alpha=0.2)+facet_grid(RCP~enviro, scales="free_x")+theme_bw()
dev.off()



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

pdf( file=paste0( outdir, "/Spatial.Preds.NoSamp.pdf" ), width = 4, height = 4 )
predict_maps2_SDF2( RCPNoSamp_SpPreds, pred_space=pred_space_rcp, pred_crop=ras, nRCP=nRCP_best )
dev.off()

