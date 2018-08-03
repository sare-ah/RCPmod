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

# Load the data
# =============
# veg example code has species and covariates in the same csv file
covariates.species <- read.csv("C:/Users/daviessa/Documents/R/Exercises/RCPmod_Foster_etal_2017/rcp-survey-artifacts/NSWVeg_covariates_species.csv", header=T, stringsAsFactors=F)

# Parameratise models
# ===================


