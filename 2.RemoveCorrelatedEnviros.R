#######################################################################################################################
# Choose environmental variables for RCP
# 
# Objective:  Remove highly correlated environmental variables
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

library(corrplot)

# Set working directory - move back to parent directory
setwd('..')

outdir<-("Models/2018-08-22/")

#Load species data with attached environmental data
spEnv<- read.csv("../Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviro.csv", header=T, stringsAsFactors=F)
spEnv<-spEnv[spEnv$area=="NCC",] #confirm NCC only

env <- names(spEnv)[171:ncol(spEnv)]

#Look at all variables correlated
corAll <- cor(spEnv[names(spEnv) %in% env], use = "pairwise.complete.obs") 
summary(corAll[upper.tri(corAll)])
corAll<-round(corAll, 2)

png(paste0(outdir,"corrPlot.png"), res=300, height=15, width=15, units="cm")
corrplot(corAll, method="square")
dev.off()

#Look at plot and find high correlations 
#Strong correlations
corTab<-data.frame(corAll)
corTab$name<-rownames(corTab)

corTab$name[abs(corTab$MBPI)>0.5] #Medium, Broad, and Fine BPI are highly correlated. Just keep Fine.
corTab$name[abs(corTab$tidalcurr_meanSummer)>0.5] #All the tidalcurrents and circ are highly correlated. Just keep summer tidal. 
corTab$name[abs(corTab$salt_meanSummer)>0.5] #All the salinities are pretty correlated, and most are correlated with temp summer and winter. Exception is min salinity.  Just keep summer temp and min salinity. 
corTab$name[abs(corTab$fetchSum)>0.5] #Sum of fetch and minimum fetch are correlated. Just keep sum of fetch. 
corTab$name[abs(corTab$dist_sand)>0.5] #mixed and sand are highly correlated. Remove mixed.

#Choose ones to keep and check new correlation plot
toKeep<-c("ArcRug","bathy","fetchSum","FBPI","salt_min","Slope","temp_meanSummer","tidalcurr_meanSummer", "dist_mud","dist_sand","dist_rock") #11 variables

corLim <- cor(spEnv[names(spEnv) %in% toKeep], use = "pairwise.complete.obs") 
summary(corLim[upper.tri(corLim)])

png(paste0(outdir,"corrPlotLim.png"), res=300, height=15, width=15, units="cm")
corrplot(corLim, method="square")
dev.off()

#Write new species-enviro table
toDrop<-env[!env %in% toKeep]
spEnv2<-spEnv[,!names(spEnv) %in% c(toDrop, "optional")]
names(spEnv2)

write.csv(spEnv2,"../Data/SpeciesMatrix_withEnviro/SpeciesByDepthSpatialEnviroRemCor.csv", row.names=F)
