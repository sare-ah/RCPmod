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
# Add code to build species profiles of each RCP (probability of species being a member of an RCP) (KG: Done, assuming you can interpret the mean "prevalence" (alpha) of the each species as its probability)
# Write documentation for this code

# Questions / Notes / In progress
##########################
# In the paper they do forward variable selection wrapped with selection of the number of RCPs. 
# I wrote "LoopModelSelection.R", which compares BICs of models run with different combinations of env vars, nRCPs, and survey/not survey to get: 
#  1) Contribution or "importance" of individual variables
#  2) The "best" model (lowest BIC), both env var terms and nRCPs. 
#  3) We're going to use the survey factor, but this code will let us compare the model fits with and without so we can report. 

# To reduce the number of env vars that are being used, 2 degree and higher orthogonal polynomials were not calculated.  1 degree orthogonal polynomials were calculated because: Hill says it's to "avoid convergence problems", and regimix throws an error with the raw data, because it's not all same order of magnitude. 

# The model outputs don't give "the most important environmental variables for distinguishing clusters" or the "effect of survey". But, the model selection script will let us compare the model fits w/ and w/o variables or survey effects and get those metrics. 

# Looks like regiboot is actually bootstrapping the "coefficients" which are the prior probabilities. 

# Using the initial priors from regimix.multifit causes very poor fit in the output mode. Using "random2" start position results in good-looking diagnostic plots. 

# One of the diagnostics, the dispersion parameter for negative binomial, doesn't work for Bernoulli distribution because there's no "disp" parameter in bernoulli output

# Question: What do the beta values coming out of regiboot object tells us? "RCPs dependence on covariates". Is looking at them useful? 

# Note: pursue difference between distance to substrate and substrate. Can the model take factors? 

# Could do 60/40 validation e.g., hold one survey back. Predicted vs observed would give a kappa statistic. 

# Try predicting mdoel onto different survey years, to see how survey predict works.  

#Do cluster analysis and random forest to compare

##########
# BEFORE RUNNING THIS SCRIPT
# Run: 
#   1. 0.LoopModelSelection
#   2. 1.ExtractAttachEnviro
#   3. 2.RemoveCorrelatedEnviros
#   4. 3.CreateEnvironmentalSpace

## Using code "LoopModelSelection.R", loop through combinations of variables and RCPs to determine:
## a) The best combination of variables (measured by BIC). Because of computational limitations, this was limited to 2-, 3-, 10-, and 11-factor combinations of the env vars
## b) The contribution or "importance" of individual variables. By comparing the 11-factor model to each 10-factor model, we can get the change in BIC when each factor is dropped. 
##    The highest change can be interpreted as the most important variable for the model. 
## c) The best number of RCPs (groups) for the data. 
## NOTE that because of computational and time limitations, the model selection was run on a subset (600 or 30%) of sites


# ## JUST TO NOTE: 
# # There are a small number of points whose depths don't match what was recorded in fcode. This stems from the spatialized points tool in the bottom patches. Jess Nephin has done a new tool to spatialize that we could use. 
#look<-read.csv("./Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv")
#look<-look[,!names(look)%in% species]
#plot(look$SourceKey,look$bathy)
#points(look$SourceKey[look$bathy<0],look$bathy[look$bathy<0], col="red")
#head(look[look$bathy<0,])
#table(look$fcode[look$bathy<0]) #mostly fine. 291 intertidal, 70 is 0-5, 12 is 5-10, 1 is 10-20 
#points(look$SourceKey[look$bathy<0&!look$fcode%in% c("ITD", "0to5")],look$bathy[look$bathy<0&!look$fcode%in% c("ITD", "0to5")], col="blue", cex=2)
#
#look[look$bathy<0&!look$fcode%in% c("ITD", "0to5"),]


#####################
# start fresh
rm(list=ls())

# Set working directory - move back to parent directory
setwd('..')

# Load helper functions
# =====================
# contains additional functions for installing packages, data transformations, plotting etc.
source("./Scripts/RCP_Helper_Functions.R")

# Make packages available
UsePackages( pkgs=c("dplyr","RCPmod", "raster", "rasterVis", "tidyr","corrplot","reshape") ) 

# Directories  
# ===========
# Subfolder for output 
outdir <- file.path("./Models", Sys.Date() )  # Set this to Models/Tday's Date/...
#outdir <- file.path("./Models/2018-09-14/" )
#suppressWarnings( dir.create( outdir, recursive = TRUE ) )

# Controls
# ========
# REQUIRED
id_vars="SourceKey" # Unique identifier for each record
sample_vars="Survey" # Field or fields that describe potential sampling effects
nstarts <- 10  # Editted this line from nstarts <- 1000
distribution <- "Bernoulli"  # Change to "Bernoulli" for P/A data or NegBin for abund
gen.start.val <- "random2"

# OPTIONAL
min.prevalence <- TRUE # specify T/F to subset data based on species prevalence
species.n <- 40  # Minimum species prevalence
subset.data <- F # Specify T/F to subset data or use all communities 
subset.size <- 600 # Specifcy random subset size --- originally set to 1000

# Load the data from previous scripts
# =============


# Load the input data
# =============
spEnv_poly_data<-readRDS(paste0(outdir,"/spEnv_poly_data.rds")) #this is the input data with transformed env vars

species <- names(spEnv_poly_data)[3:161] #159 species
enviro.variables<-names(spEnv_poly_data)[162:ncol(spEnv_poly_data)] #11 enviro

n.sites <- nrow(spEnv_poly_data) # Number of sites in study area
site.names <- spEnv_poly_data[,id_vars] # Record site order - we may need this later if we are subsetting the data 

# Conditional processing to subset sites used - ensure sitename column is correctly specified in controls
if (subset.data) {
  # specify the  subset to use in the analysis
  set.seed(subset.size)
  sample.sites <- sample(site.names, subset.size)
  spEnv_poly_data <- spEnv_poly_data[spEnv_poly_data[,id_vars] %in% sample.sites,]
  print(paste0("Successfully subsetted [",subset.size,"] random sites"))
  site.names<-sample.sites
} else {
  print("No subsetting performed")
}

# Conditional processing to subset number of species used based on species prevalence
if (min.prevalence) {
  # Determine in decreasing order the total count of each species within the study area
  species.count <- data.frame(count=sort(colSums(spEnv_poly_data[,names(spEnv_poly_data) %in% species]), decreasing=T))
  species.count$species <- row.names(species.count)
  # Select the species to be dropped and not used in the model
  dropped.species <- species.count$species[species.count$count < species.n]
  # Remove species columns from covariates.species that are not prevalent enough - otherwise poly_data() fails!
  spEnv_poly_data <- spEnv_poly_data[ , !(names(spEnv_poly_data) %in% dropped.species)]
  # Remove species names from species vector - otherwise poly_data() fails!
  species <- setdiff(species, dropped.species)
}

length(species) #115 species
length(enviro.variables) #3 enviro
length(site.names) #1994

# Create RCP model using outputs of LoopModelSelection
varsFromBestModel<-gsub("1","",c(enviro.variables))
nRCP_best<-6

form <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",c(paste(c(paste0(varsFromBestModel, "1")), collapse="+"))))

# Run best model (for full output)
# ==================================
control <- list( optimise=FALSE, quiet=FALSE)

# #With "random2" start
RCPmod_surv <- regimix(form.RCP=form, form.spp=paste0("~ ", sample_vars),data=spEnv_poly_data, dist=distribution, nRCP=nRCP_best, inits = "random2", control=control)

saveRDS(RCPmod_surv, file=paste0(outdir,"/RCPmod_surv.rds"))


# Plot model diagnostics
# ======================
# Residual plots - Species Level
png( file=paste0(outdir,"/Diagnostics/Residuals_Samp_random2.png"), width = 15, height = 10, units="cm", res=200)
plot.regimix(RCPmod_surv, type="RQR", fitted.scale="log")
dev.off()

# Residual plots - Site Level
png( file=paste0(outdir,"/Diagnostics/Residuals_Samp_Site_random2.png"), width = 15, height = 10, units="cm", res=200)
plot.regimix(RCPmod_surv, type="deviance", fitted.scale="log")
dev.off()


# # Cooks Distance Plots
# # Not 100% sure how to interpret
# #tmp <- stability.regimix(RCPmod_surv, oosSizeRange=NULL, mc.cores=1, times=RCPmod_surv$n, doPlot=FALSE) # original
# tmp <- stability.regimix(RCPmod_surv, oosSizeRange=seq( from=1,to=RCPmod_surv$n%/%5,length=5), mc.cores=1, times=RCPmod_surv$n, doPlot=FALSE) # quick run (3 runs instead of 200) #this took 2.5 days to run halfway. Killed it to get to the more interesting results
# #tmp <- stability.regimix(RCPsamp_fin, oosSizeRange=c(1,2,5,10,20), mc.cores=1, times=RCPsamp_fin$n, doPlot=FALSE) # simplified, but missing horizontal line
# 
# saveRDS( tmp, file=paste0(outdir,"/Fish.Cooks.Dist_quick.rds" ))
# 
# png( file=paste0(outdir,"/Cooks.Distance.Samp.png"),width = 20, height = 10, units="cm", res=200)
# plot( tmp, minWidth=2, ncuts=111 )
# dev.off()

# Generate bootstrap estimates of parameters
# ==========================================
# Warning --- bootstrap is slooow
# rcpsamp_boots<-regiboot(RCPsamp_fin, type="BayesBoot", nboot=1000, mc.cores=1) # original code
RCP_boot<-regiboot(RCPmod_surv, type="BayesBoot", nboot=10, mc.cores=1) # simplified
saveRDS(RCP_boot, file=paste0(outdir,"/RCP_boot.rds" ))


# Evaluate sampling factor effect
# ===============================
# Calculate average, SD and CI of species abundances in each RCP, and plot results
# Sp_abund_gen(bootstrap object, input species file with levels)
RCP_sp_abund<-sp_abund_gen(RCP_boot, spEnv_poly_data)
saveRDS(RCP_sp_abund, paste0(outdir, "/RCP_sp_abund.rds"))

#sp_abund_plot(sp_abund, bySampleEffect, sampleEffectName, sortByRCP, refRCP=NULL)
png(paste0(outdir,"/Diagnostics/RCP_abund_sp_surv.png"), height=20, width=20, res=300, units="cm", pointsize = 12)
sp_abund_plot(RCP_sp_abund,bySampleEffect=TRUE,"Survey", sortByRCP = F)
dev.off()

png(paste0(outdir,"/Diagnostics/RCP_abund_sp_nosurv.png"), height=12, width=20, res=300, units="cm", pointsize = 12)
sp_abund_plot(RCP_sp_abund,bySampleEffect=FALSE, sortByRCP = T)
dev.off()

# Plot of sampling factor effects (relative to base case)
# Flag - code could be generalized to make sure the factor levels are right. Not 100% sure here because the "base case" isn't listed. 
png( file=paste0(outdir,"/Diagnostics/Sampling.Factor.Effects.png"), height=25, width=20, res=300, units="cm", pointsize = 12)
sampling_dotplot2(RCPmod_surv,RCP_boot,legend_fact=c("PAC 2014-29", "PAC 2014-58","PAC 2015-52"), col=c("black", "red", "blue"), lty=c(1,2,3))
dev.off()


# Spatial Predictions
# ===================

#Load environmental space
#rcp_poly_pred<-readRDS(paste0(outdir,"/rcp_poly_pred.rds"))
envRasterPoints_poly_level<-readRDS(paste0(outdir,"/envRasterPoints_poly_level_P1552.rds"))
envRasterPoints<-readRDS(paste0(outdir,"/envRasterPoints.rds"))
envRasterPoints_poly<-readRDS(paste0(outdir,"/envRasterPoints_poly.rds"))

# #This command maxes out my ram, but works on GIS computer. Looks like Fig 4 in hill paper. 
# #This is a probability of prediction for every raster cell in the rasSub stack. 
# #Predict.regimix takes "newdata" of env+sp variables and calcualates probability of RCP membership. 
# #it does this for the single model output (object) and the boot output (object2)
# RCPsamp_preds<-readRDS("D:/Documents/Projects/WorldClass_DiveSurveys/RCP/Models/2018-08-13/RCPsamp_preds.rds")
RCPsamp_preds1<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=envRasterPoints_poly_level[1:500000,])
saveRDS(RCPsamp_preds1, paste0(outdir,"/RCPsamp_preds1.rds"))
RCPsamp_preds2<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=envRasterPoints_poly_level[500001:1000000,])
saveRDS(RCPsamp_preds2, paste0(outdir,"/RCPsamp_preds2.rds"))
RCPsamp_preds3<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=envRasterPoints_poly_level[1000001:5000000,])
saveRDS(RCPsamp_preds3, paste0(outdir,"/RCPsamp_preds3.rds"))
RCPsamp_preds4<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=envRasterPoints_poly_level[5000001:10000000,])
saveRDS(RCPsamp_preds4, paste0(outdir,"/RCPsamp_preds4.rds"))
RCPsamp_preds5<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=envRasterPoints_poly_level[10000001:nrow(envRasterPoints_poly_level),])
saveRDS(RCPsamp_preds5, paste0(outdir,"/RCPsamp_preds5.rds"))

RCPsamp_preds<-RCPsamp_preds1
for (i in 1:3){
RCPsamp_preds[[i]]<-rbind(RCPsamp_preds1[[i]],RCPsamp_preds2[[i]],RCPsamp_preds3[[i]],RCPsamp_preds4[[i]],RCPsamp_preds5[[i]])
}
lower<-rbind(RCPsamp_preds1[[4]][,,1],RCPsamp_preds2[[4]][,,1],RCPsamp_preds3[[4]][,,1],RCPsamp_preds4[[4]][,,1],RCPsamp_preds5[[4]][,,1])
upper<-rbind(RCPsamp_preds1[[4]][,,2],RCPsamp_preds2[[4]][,,2],RCPsamp_preds3[[4]][,,2],RCPsamp_preds4[[4]][,,2],RCPsamp_preds5[[4]][,,2])
RCPsamp_preds[[4]][,,1]<-NULL
RCPsamp_preds[[4]][,,2]<-NULL
RCPsamp_preds[[4]][,,1]<-lower
RCPsamp_preds[[4]][,,2]<-upper
# CIs not working but boot and pt preds are fine.

saveRDS(RCPsamp_preds, paste0(outdir,"/RCPsamp_preds.rds"))

## Create output rasters of predicted RCPs at each raster cell
predictClasses(envRasterPoints, RCPsamp_preds, paste0(outdir,"/RCP_Output_Probabilities.tif"),"prob")
predictClasses(envRasterPoints, RCPsamp_preds, paste0(outdir,"/RCP_Output_Classes.tif"),"classes")

##### 
# Get response of each RCP to environmental variables (Figure 3 in Hill paper)
## giving predict.regimix the input sampling sites (here, rcp_data) gives probability of RCP membership for each site. 
RCP_preds_intialSites<-predict.regimix(object=RCPmod_surv, object2=RCP_boot, newdata=spEnv_poly_data)

#Get the point predictions from predict.regimix
## NOTE: Should this be the bootstrap predictions?
pt_preds<-as.data.frame(RCP_preds_intialSites[[1]])

#Bring in original environmental data that wasn't transformed
#Using the non-orthogonalized data looks identical to orthogonalized data, but the units are more meaningful on the plots. 
spEnvOrig<-read.csv("./Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv")
spEnvOrig<-spEnvOrig[,names(spEnvOrig)%in%c("SourceKey",gsub("1","",enviro.variables))]
spEnv_poly_data_wOrig<-merge(spEnv_poly_data, spEnvOrig, by="SourceKey")

#Reorganize data to have the value of each input point's environmental variables and the probability it falls into each RCP
#Survey doesn't matter here, because each of the input points only has one associated survey
#This is taking the input points and calculating the probability of falling in an RCP, given the model, given that sites' env space  
spEnv_melt<-melt(spEnv_poly_data_wOrig[,!names(spEnv_poly_data_wOrig)%in% species], id.vars=c("SourceKey","Survey"))
summary(spEnv_melt$SourceKey==spEnv_melt$SourceKey) #order of points is the same. attach the spatial predictions
spEnv_wPreds<-cbind.data.frame(spEnv_melt, pt_preds)
spEnv_wPreds_melt<-melt(spEnv_wPreds, id.vars=c("SourceKey","Survey", "variable","value"))
names(spEnv_wPreds_melt)[3:6]<-c("enviro","enviro_val","RCP","RCP_prob")
spEnv_wPreds_melt<-spEnv_wPreds_melt[spEnv_wPreds_melt$enviro%in%gsub("1","",enviro.variables), ] #remove orthogonal values

#Plot - matches Figure 3 in the paper. 
png(paste0(outdir,"/Diagnostics/envirospace_unpoly.png"),res=400, width=50, height=20, units="cm")
ggplot(data=spEnv_wPreds_melt, aes(x=enviro_val, y=RCP_prob, col=RCP))+geom_point(alpha=0.2)+facet_grid(RCP~enviro, scales="free_x")+theme_bw()
dev.off()

## Get "indicator species" for each RCP (needs to be from predict to get variation)
inds<-RCPindSp(spEnv_poly_data, RCP_preds_intialSites, indvalCutoff=0.15, pValCutoff=0.05, id_var="SourceKey", species=species)
inds

save.image(paste0(outdir,"/RCP_Run_bestFormWithSurv_2018.09.17.Rdata")) 

###
# what else is there to do? 
# finalize the model selection
# ? convergence
# - compare 1 vs 2 or 5% species cutoff? 



