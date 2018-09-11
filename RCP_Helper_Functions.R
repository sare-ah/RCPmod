#########################################################################################
### Helper Functions for RCPmod used in Foster, Hill & Lyons 2016 JRSS paper          ###
### Written by N. Hill July 2015.                                                     ###
### Works with RCPmod version 2.88                                                    ###
#########################################################################################

#### load packages----
# Install missing packages and load required packages (if required)
UsePackages <- function( pkgs, update=FALSE, locn="http://cran.rstudio.com/" ) {
  # Identify missing (i.e., not yet installed) packages
  newPkgs <- pkgs[!(pkgs %in% installed.packages( )[, "Package"])]
  # Install missing packages if required
  if( length(newPkgs) )  install.packages( newPkgs, repos=locn )
  # Loop over all packages
  for( i in 1:length(pkgs) ) {
    # Load required packages using 'library'
    eval( parse(text=paste("library(", pkgs[i], ")", sep="")) )
  }  # End i loop over package names
  # Update packages if requested
  if( update ) update.packages( ask=FALSE )
}  # End UsePackages function

#### poly_data----
# Generate orthogonal polynomials for RCPmod input. Useful to avoid convergence problems
# Save polynomial bases to transform prediction space

poly_data<-function(poly_vars,      #vector of predictor variable names to create orthogonal polynomials. Matched to column names in 'data'
                    degree,         #vector of same length as poly_vars specifying the polynomial degree for each predictor variable in poly_vars
                    id_vars,        #vector of ID variable names (Not transformed). Matched to column names in 'data'
                    sample_vars=NULL, #vector of sampling level variable names (e.g. Gear type). Matched to column names in 'data'
                    species_vars,    #vector of response variables names. Matched to column names in 'data'
                    offset=NULL,     #name of offset variable.  Matched to column names in 'data'
                    data, ...)
  {
  store_polys<-list()
  for(i in 1:length(poly_vars)){
    store_polys[[i]]<-poly(data[,poly_vars[i]], degree=degree[i])
    dimnames(store_polys[[i]])[[2]]<-paste0(poly_vars[i],seq(1:degree[i]))
  }
  names(store_polys)<-poly_vars
  
  rcp_data<-na.omit(cbind(subset(data, select= c(id_vars, sample_vars, offset, species_vars)),do.call(cbind, store_polys)))
  
  return(list(rcp_data=rcp_data, poly_output=store_polys))
}


### poly_pred_space----
# Transforming prediction space using same basis as orthogonal polynomials used to build models
# Note: only accomodates one sampling term at the moment
# Note: offset isn't actually used in predict function that predicts RCP membership, but will keep as it might be useful to predict expected abundance of species at site.

poly_pred_space<-function(pred_space,             #dataframe containing variables (ras to points)
                          poly_output,            #extracted list of stored polynomial attribtutes from 'poly_data' function
                          offset_val=NULL,        #an offset value. Possibly mean of offset used in model building. Will be logged within function
                          offset_name=NULL,       #name of offset used in RCP model
                          sampling_vals=NULL,     #level of sampling factor for prediction 
                          sampling_name=NULL,     #name of sampling factor used in RCP model
                          sampling_factor_levels=NULL )  #levels of sampling factor used in RCP model
  {
  
  # transform predictors using saved orthogonal polynomial attributes
  pred_polys<-list()
  vars<-names(poly_output)
  for( i in 1: length(vars)){
    pred_polys[[i]]<- predict( poly_output[[i]], pred_space[, names(pred_space) %in% vars[i]]) 
    dimnames(pred_polys[[i]])[[2]]<-dimnames(poly_output[[i]])[[2]]
  }
  pred_polys_df<-as.data.frame(do.call(cbind, pred_polys)) #data frame of polynomial-ized environmental data for footprint of raster stack
  
  #create offset term
  if(!is.null(offset_val)){
    pred_polys_df$offset<-log(offset_val)
    names(pred_polys_df)[ncol(pred_polys_df)]<-paste0("log(", offset_name, ")")
  }
  
  # create sampling variable. 
  # only accommodates one sampling factor
  if(!is.null(sampling_vals)){
    reps<- length(sampling_vals)
    pred_polys_df$sampling<-factor(sampling_vals,levels=sampling_factor_levels)
    names(pred_polys_df)[ncol(pred_polys_df)]<-sampling_name
  }
  
  return(pred_polys_df)
}



#### Calculate average SD and CI of species abundances in each RCP in each level of sampling factor.
# Updated by K. Gale to be generalized to any input/any sampling levels 
# Requires:
# boot_obj= regiboot object
# covariates.species= the species input table that contains a field for the sampling level (sample_vars, defined at start of RCP Mod)
# Right now cannot do a version without the sampling level

Sp_abund_gen<-function(boot_obj, covariates.species) {
  
  require(tidyr)
  invlogit<-function(x){exp(x)/(1+exp(x))}
  
  levels<-data.frame(sample_var=unique(covariates.species[,sample_vars]))
  levels$sample_var_short<-tolower(gsub("[[:punct:]]|\\s","",levels$sample_var)) #drop all special characters or spaces from level name. 
  levels$sample_var_short<-gsub("\\s","",levels$sample_var_short)
  results_levels<-paste0("res_",levels$sample_var_short)
  
  #set up coefficient extraction
  taus<-grepl("tau",dimnames(boot_obj)[[2]]) #for each species, difference from the mean (alpha) to mean for each RCP, for the first n-1 RCPs. 
  #The last RCP means (tau) are calculated by the "sum to zero" constraint (negative sum of taus for the other RCPs, so sum(tau)=0 for each species )
  alphas<-grepl("alpha",dimnames(boot_obj)[[2]]) #mean species prevalence overall in the dataset, on the logit scale
  gammas<-grepl("gamma",dimnames(boot_obj)[[2]]) #Specifies the species dependence on the covariates
  
  results<-rep( list(list()), length(results_levels) ) 
  
  for (n in 1:length(results_levels)){
    names(results)[n]<-results_levels[n]
  }
  
  #run loop for each row in boot object
  for(i in 1:dim(boot_obj)[1]){
    
    #extract and reformat coeficients 
    #tau - "regional profiles"
    temp_tau<-data.frame(value=boot_obj[i,taus])
    rownames(temp_tau)<-gsub("^A_","A-",rownames(temp_tau)) #This is specific to the dive dataset - not super generalized
    rownames(temp_tau)<-gsub("^I_","I-",rownames(temp_tau))
    temp_tau$species<-sapply(strsplit(rownames(temp_tau),"_"), "[", 1)
    temp_tau$coef<-sapply(strsplit(rownames(temp_tau),"_"), "[", 3)
    
    #"Base Case" of the model. 
    tau_base<-spread(temp_tau[,1:3], species, value) #convert to wide-mode table by "coef"
    tau_base<-rbind( tau_base[,-1], -colSums( tau_base[,-1])) #add a row of the negative sum of tau for RCP (row), such that the sum of tau for each species across the RCPs is 0.
    
    #alpha- OK as is. Mean of species prevalence (logit) overall. 
    temp_alphas<-boot_obj[i,alphas]
    names(temp_alphas)<-gsub("^A_","A-",names(temp_alphas)) #not super generalized
    names(temp_alphas)<-gsub("^I_","I-",names(temp_alphas))
    
    
    #gamma. Species dependence on covariates.
    temp_gamma<-data.frame(value=boot_obj[i,gammas])
    rownames(temp_gamma)<-gsub("^A_","A-",rownames(temp_gamma)) #not super generalized
    rownames(temp_gamma)<-gsub("^I_","I-",rownames(temp_gamma))
    
    temp_gamma$species<-sapply(strsplit(rownames(temp_gamma),"_"), "[", 1)
    temp_gamma$coef<-sapply(strsplit(rownames(temp_gamma),"_"), "[", 3)
    temp_gamma$sample_var<-gsub(paste(sample_vars),"",temp_gamma$coef)
    temp_gamma<-merge(temp_gamma, levels, by="sample_var")
    
    splitCoef<-unique(temp_gamma$sample_var_short)
    gamma_list<-list()
    for (n in 1:length(splitCoef)){
      gamma_list[[n]]<-temp_gamma[temp_gamma$sample_var_short==splitCoef[n],"value"]
      names(gamma_list)[[n]]<-paste0("gamma_",unique(temp_gamma[temp_gamma$sample_var_short==splitCoef[n],"sample_var_short"]))
    }
    
    #Identify Base Case -  the one that isn't in the gamma list
    base<-results_levels[!results_levels %in% gsub("gamma", "res",names(gamma_list))]
    notbase<-results_levels[results_levels %in% gsub("gamma", "res",names(gamma_list))]
    
    #calulate base case values (the level of the model for which there is no offset)
    lps <- sweep( tau_base, 2, temp_alphas, "+") #add alpha (sp mean prevalence) to base taus (difference from mean in each RCP) 
    
    #Add base case prevalence to results list
    
    
    results[names(results)==base][[1]][[i]] <- as.matrix(round(invlogit( lps),2)) #make pretty. This is the species prevalence in each RCP in base. ###Changed this and the following from exp(lps) to exp(lps)/(1-exp(lps)) for the p/a dataset -- are all the alpha outputs on the logit scale? 
    
    
    #calculate values for each sampling level
    for (x in 1:length(notbase)){
     #add each level's gamma (species dependence on level) to all base case prevalences 
      dat<-gamma_list[names(gamma_list)==gsub("res","gamma",notbase[x])][[1]]
      results[names(results)==notbase[x]][[1]][[i]] <- as.matrix(round(invlogit(sweep( lps, 2, gamma_list[names(gamma_list)==gsub("res","gamma",notbase[x])][[1]],"+")),2)) #make pretty. This is the species prevalence in each RCP in spring.
     }
    }
  
  #compile summaries of results, summarized across boot objects. These are mean & summary stats for species prevalences in each RCP in each level
  results_summary<-rep( list(list()), length(results_levels) ) 
  for (n in 1:length(results_levels)){
    names(results_summary)[n]<-results_levels[n]
  }
  
  for (n in 1:length(results_levels)){
    results_summary[names(results_summary)==results_levels[n]][[1]]<-list(mean=apply(simplify2array(results[names(results)==results_levels[n]][[1]]), c(1,2), mean),
                                                                          sd=apply(simplify2array(results[names(results)==results_levels[n]][[1]]), c(1,2), sd),
                                                                          lower=apply(simplify2array(results[names(results)==results_levels[n]][[1]]), c(1,2), function(x) quantile(x, probs=0.025)),
                                                                          upper=apply(simplify2array(results[names(results)==results_levels[n]][[1]]), c(1,2), function(x) quantile(x, probs=0.975)))
  }
  
  last_list<-list()
  for (i in 1:length(results_levels)){
    test<-data.frame(t(do.call("cbind.data.frame",results_summary[[i]])))
    test$species<-gsub("^sd\\.|^lower\\.|^upper\\.|^mean\\.","",row.names(test)) 
    test$stat<-substr(row.names(test),1,2) #should be an easier way I just can't remeber right now
    test$stat<-ifelse(test$stat=="me","mean",ifelse(test$stat=="lo","lower",ifelse(test$stat=="up","upper",ifelse(test$stat=="sd","sd",""))))
    test$sample_var<-gsub("res_","",names(results_summary)[i])
    last_list[[i]]<-test
  }
  final<-do.call("rbind",last_list)  
  
  final_melt<-reshape::melt(final, id.vars=c("species","sample_var","stat"))
  names(final_melt)[4]<-"RCP"
  final_melt$value<-round(final_melt$value,2)
  final_cast<-reshape::cast(final_melt, species+RCP+sample_var~stat)
  final_cast$RCP<-gsub("X","",final_cast$RCP)
  
  return(final_cast)}

#### 
# Plot output of sp_abund_gen to get predicted prevalence of each species in each RCP. 
# Option bySampleEffect specifies whether results are broken down by Sampling Effect (E.g., survey)
# If bySampleEffect=T, mean(Prevalence) from bootstrap are averaged across Sample Effect Levels (mean + CIs of bootstrap mean shown)
# if bySampleEffect=F, mean(Prevalence) are used straight from bootstrap
# sortByRCP = TRUE will sort by a specified RCP, or by the first one if refRCP is blank
# sp_abund: output of sp_abund_gen

sp_abund_plot<-function(sp_abund, bySampleEffect, sampleEffectName, sortByRCP, refRCP=NULL){
require(ggplot2)
  require(plyr)
  if(bySampleEffect==F){
    sp_abund<-ddply(sp_abund, c("species", "RCP"), summarize,meanv=mean(mean), lowerv=quantile(mean, probs=0.025),upperv=quantile(mean, probs=0.975))
      
    if (sortByRCP==T){
        if (is.null(refRCP)==T){  
          forOrder<-sp_abund[sp_abund$RCP==sp_abund$RCP[1],]
              } else {
          forOrder<-sp_abund[sp_abund$RCP==sp_abund$RCP[refRCP],]  
                        }
      forOrder$species<-factor(forOrder$species, levels = forOrder$species[order(-forOrder$meanv)])
      sp_abund$species <- factor(sp_abund$species, levels = levels(forOrder$species))
          } else {}
        
    ggplot(data=sp_abund, aes(x=species, y=meanv, color=RCP))+
      facet_grid(. ~ RCP)+geom_point()+coord_flip()+
      geom_linerange(aes(ymin=lowerv,ymax=upperv))+
      theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_x_discrete(limits=rev(levels(factor(sp_abund$species))))+
      scale_y_continuous(expand = c(0.0, 0))+xlab(c("Species"))+ylab(c("Mean Prevalence"))+
      theme(legend.position="none")
    
  }else {
    
    if (sortByRCP==T){
      if (is.null(refRCP)==T){  
        forOrder<-sp_abund[sp_abund$RCP==sp_abund$RCP[1]&sp_abund$sample_var==sp_abund$sample_var[1],]
      } else {
        forOrder<-sp_abund[sp_abund$RCP==sp_abund$RCP[refRCP]&sp_abund$sample_var==sp_abund$sample_var[1],]  
      }
      forOrder$species<-factor(forOrder$species, levels = forOrder$species[order(-forOrder$mean)])
      sp_abund$species <- factor(sp_abund$species, levels = levels(forOrder$species))
    } else {}
  
    ggplot(data=sp_abund, aes(x=species, y=mean, color=sample_var))+
    facet_grid(sample_var ~ RCP)+geom_point()+coord_flip()+
   geom_linerange(aes(ymin=lower,ymax=upper))+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_x_discrete(limits=rev(levels(factor(sp_abund$species))))+
    scale_y_continuous(expand = c(0.0, 0))+xlab(c("Species"))+ylab(c("Mean Prevalence"))+
   labs(color=sampleEffectName)
  }
}

## Indicator Species

# require(labdsv)
# require(reshape2)
# 
# 
# 
# sp_abund1<-data.frame(sp_abund[,c("species","mean","RCP")])
# sp_melt<- melt(sp_abund1, id.vars=c("species","RCP"), measure.vars=c("mean")) 
# sp_cast<-cast(sp_melt, RCP~species, fun="mean")
# 
# ind<-indval(sp_cast[2:ncol(sp_cast)], sp_cast$RCP)
# 
# indvals<-melt(data.frame(species=row.names(ind$indval),ind$indval))
# abund<-melt(data.frame(species=row.names(ind$relabu),ind$relabu))
# maxcl<-melt(data.frame(species=names(ind$maxcls),ind$maxcls))
# 
# data.frame

#### sampling_dotplot---
# Produce dotplot of sampling-level coefficients (and 95% CI) for each species

sampling_dotplot2<-function(best_mod,              # output of regimix function for final model
                            boot_obj,              # output of regiboot function for final model
                            legend_fact,           # levels of categorical sampling variable to plot. Coefficients relative to first level of factor. 
                                                   # So usually 2:n levels(sampling_variable)
                            CI=c(0.025, 0.975),    # confidence interval to plot
                            col="black",           # colour/s of dots and CI lines. Specified in same way as col is usually specified
                            lty=1)                 # lty= line type of CI lines. Specified in same way as lty is usually specified
{
  
  require(lattice)
  # gammas<- paste0("gamma", 1: (length(Species)*length(sampling_names)))
  gammas<-grepl("gamma",dimnames(boot_obj)[[2]])
  temp_dat<-boot_obj[,gammas]
  
  colnames(temp_dat)<-gsub("^A_","A-",colnames(temp_dat)) #not super generalized
  colnames(temp_dat)<-gsub("^I_","I-",colnames(temp_dat))
  temp<-data.frame(avs=as.numeric(unname(colMeans(temp_dat))),
                   t(apply(temp_dat, 2, quantile, probs=CI)),
                   #sampling_var=rep(sampling_names, each=length(length(best_mod$names$spp))),
                   sampling_var=sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 3),
                   #Species=factor(rep(best_mod$names$spp,length(sampling_names))))
                   Species=factor(sapply(strsplit(dimnames(temp_dat)[[2]],"_"), "[", 1)))
  
  names(temp)[2:3]<-c("lower", "upper")
  temp$Species<-gsub("."," ", temp$Species, fixed=TRUE) #get rid of '.' in species names
  
  temp$Species<-as.factor(temp$Species) #convert back to factor
  temp$Species <- factor(temp$Species, levels=rev(levels(temp$Species))) 
  
  trellis.par.set(superpose.symbol=list(pch=16,col=col, cex=1.2),
                  superpose.line=list(col="transparent"))
  dotplot(Species ~ avs, groups=sampling_var, data=temp, cols=col, lty=lty, low=temp$lower, high=temp$upper, subscript=TRUE, 
          auto.key=list(space="top", columns=2, cex=1.4, text=legend_fact),
          ylab=list(cex=1.4), xlab=list("Coefficient",cex=1.4),
          scales = list(tck = c(1, 0), x=list(cex=1.2), y=list(cex=1.2)),
          prepanel = function(x, y, ...) { list(xlim=range(temp$lower, temp$upper)) }, 
          panel=panel.superpose,
          panel.groups=function(x, y, subscripts, group.number, cols, low, high, ...)
          {
            if(group.number==1) jiggle <- 0.1 else jiggle <- -0.1
            panel.abline(v=0, lty=2)
            panel.abline(h=1:length(best_mod$names$spp), col.line="light grey", lty=1)
            panel.dotplot(x, y+jiggle, group.number, ...)
            panel.arrows(low[subscripts], y+jiggle, high[subscripts], y+jiggle, code=3, angle=90, 
                         length=0.05, col=cols[group.number], lty=lty[group.number])
            #panel.segments(temp$lower, y+jiggle, + temp$upper, y+jiggle, lty = lty, col =col, lwd=2, cex=1.2)
          })            
  
}



#####predict_maps2_SDF2----
# Plot RCP predictions
# S. Foster version, modified for vertical instead of horizontal RCP layout




predict_maps2_SDF2 <- function(predictions,     #output from predict.regimix
                               pred_space,      #dataframe containing coordinates for the prediction space
                               pred_crop,       #raster of extent of prediction space (used in rasterize function)
                               nRCP,            #the number of RCPs
                               my.ylim=NULL,    #y limit to plot
                               my.xlim=NULL,    #x limit to plot
                               my.asp=1)        #aspect for plotting
  {
  require(rasterVis)

  preds<-predictions
  av_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootPreds)), pred_crop)
  av_pred <- dropLayer( av_pred, 1)
  low_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootCIs[1:nrow(pred_space),1:nRCP,1])), pred_crop)
  low_pred <- dropLayer( low_pred, 1)
  upp_pred<-rasterize(SpatialPointsDataFrame(coords= pred_space[,1:2], data=as.data.frame(preds$bootCIs[1:nrow(pred_space),1:nRCP,2])), pred_crop)
  upp_pred <-  dropLayer( upp_pred, 1)
  
  colour <- c("#dddddd","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d")#,"#000000")
  breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  
  nPlots <- nRCP*3
  layMat <- matrix( NA, ncol=3+1, nrow=nRCP+2)
  layMat[1,1] <- 0
  layMat[1,-1] <- 1:3
  layMat[-1,1] <- 3+1:(nRCP+1)
  layMat[nrow( layMat),1] <- 0
  tmp <- matrix( 2 + (nrow( layMat)-1) + 1:(3*nRCP), ncol=3, byrow=TRUE)
  layMat[-(c(1,nrow( layMat))),-1] <- tmp
  layMat[nrow( layMat), -1] <- max( layMat, na.rm=TRUE) + 1
  #good grief!
  par( oma=c(1,1,1,1))
  layout( layMat, heights=c(0.15, rep(1, nRCP), 0.5), widths=c(0.5, rep( 1, 3)))
  #first plot (empty)
  
  #second plot (lower CI label)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Lower CI",cex=1.5)
  #third plot (pt pred)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Point Prediction",cex=1.5)
  #fourth plot (upper CI label)
  par( mar=c( 0,0,0,0))
  plot.new()
  text(0.5,0.5,"Upper CI",cex=1.5)
  #Next row labels
  for( ii in 1:nRCP){
    par( mar=c( 0,0,0,0))
    plot.new()
    text( 0.5,0.5, paste("RCP", ii, sep=" "), srt=90, cex=1.5)
  }
  #  plot.new()
  
  par( mar=c( 0.7,0.7,0,0)+0.1)
  for( ii in 1:nRCP){
    my.ylab <- c( paste( "RCP ", ii, sep=""), "", "")
    my.yaxt <- c( "s", "n", "n")
    my.xlab <- rep( "", 3)
    if( ii == 1)
      my.xaxt <- c( "n", "n", "n")
    else if( ii==nRCP)
      my.xaxt <- rep( "s", 3)
    else
      my.xaxt <- rep( "n", 3)
    tmpDat <- subset( low_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[1], yaxt=my.yaxt[1], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
    tmpDat <- subset( av_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[2], yaxt=my.yaxt[2], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
    tmpDat <- subset( upp_pred, ii)
    if( is.null( my.xlim))
      my.xlim <- range( coordinates( tmpDat)[,1])
    if( is.null( my.ylim))
      my.ylim <- range( coordinates( tmpDat)[,2])
    image( tmpDat, breaks=breaks, col=colour, ylab="", xlab="", main="", asp=my.asp, xaxt=my.xaxt[3], yaxt=my.yaxt[3], xlim=my.xlim, ylim=my.ylim,
           zlim=c(0,1))
    box()
  }
 
  plot.new()
  require( fields)#for image.plot
  image.plot( x=breaks, y=breaks, z=matrix( seq( from=0, to=1, length=length( breaks)^2), ncol=length( breaks)), horizontal=TRUE, legend.only=TRUE, legend.shrink=0.8, col=colour, smallplot=c(0.05,0.95,0.55,0.7), axis.args=list( at=seq( from=-0.055, to=1.055, length=6),#c(-0.05,0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05),# c(0,0.25,0.5,0.75,1),
                                                                                                                                                                                                                                    labels=seq( from=0, to=1, by=0.2),#c(0,0.25,0.5,0.75,1),
                                                                                                                                                                                                                                    cex.axis=1.1), legend.args=list(text='Probability of RCP membership', side=1, line=1.9, cex=0.8))
  
  return( NULL)
}
