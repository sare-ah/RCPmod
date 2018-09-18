#######################################################################################################################
# Identify best RCP model
# 
# Objective:  Loop through iterations of environmental variables, number of RCPs, sampling effect (yes/no), and rare species cutoff (none, 5%)
#             to find model with lowest BIC. Also finds "variable importance" by comparing models with single enviromental variables removed
#
# Author:     Sarah Davies and Katie Gale
#             Sarah.Davies@dfo-mpo.gc.ca/Katie.Gale@dfo-mpo.gc.ca
# Date:       August 2, 2018
######################################################################################################################

library(combinat)
library(RCPmod)

##############################################################
#Load species-enviro dataset and set up config info
##############################################################

setwd("F:/LoopModelSelection/2018-09-10")
rcp_data<-readRDS("rcp_data.rds")
head(rcp_data)

id_vars="SourceKey" # Unique identifier for each record
siteNo <- rcp_data[,id_vars] 
species <- names(rcp_data)[3:161] #159 species
env<-names(rcp_data)[162:ncol(rcp_data)] #11 enviro

sample_vars="Survey" # Field or fields that describe potential sampling effects
nstarts <- 10  # Editted this line from nstarts <- 1000
#max.nRCP <- 6   # Not used when manually testing number of RCPs
distribution <- "Bernoulli"  # Change to "Bernoulli" for P/A data or NegBin for abund
gen.start.val <- "random2"

# OPTIONAL
min.prevalence <- T # specify T/F to subset data based on species prevalence
species.n <- 40  # Minimum species prevalence  (20 is 1%)
subset.data <- T # Specify T/F to subset data or use all communities 
subset.size <- 600 # Specifcy random subset size --- originally set to 1000. 600 is 30%


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

length(species) #38 species (159 full; 133 with 1% cutoff; 119 with 1% and 1000 site cutoff; 70 with 2% and 600 site cutoff)
length(env) #11 enviro (11 full)
env<-env[!env%in%c("salt_min1", "dist_sand1","Slope1")]
length(site.names) #200 (1994 full; 1000 or 600 with cutoff)


##############################################################
## Create all possible iterations of environmental parameters
##############################################################

# Get list of all combinations of environmental variables
list1<-list()
for(i in 1:length(env)){
  list1[[i]]<-list(combn(env, i, fun=NULL, simplify=F)) }

# list2<-list()
# for (i in 2:length(env)){  #Exclude 1-variable solutions
#   list2[[i]]<-apply(do.call("rbind",list1[[i]][[1]]),1,paste,collapse="+") }

list2<-list()
#for (i in c(2:3, (length(env)-1):length(env))){  #Exclude 1, 4:9-variable solutions = 232 versions
  for (i in c(2:3, (length(env)-1):length(env))){
  list2[[i]]<-apply(do.call("rbind",list1[[i]][[1]]),1,paste,collapse="+") }

iterations<-unlist(list2)



allForms<-list()
for (n in 1:length(iterations)){
  allForms[[n]] <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",c(paste(iterations[n], collapse="+")))) }

length(allForms) #18 env is 262125 model forms. 9 is 511. 4 is 11. 11 is 2036. Removing 1 and 4:9 knocks it down to 232. This leaves 2 and 3 solutions, as well as 10 and 11 to get variable importance.

min.nRCP<-6
max.nRCP<-9

#AllFormula table
options<-data.frame(form=as.character(paste(allForms)))
for (x in 1:nrow(options)){ options$form2[x]<-strsplit(as.character(options$form[x]), "\\~\\s")[[1]][2]}

nRCP_list<-list()      
formula<-list()

for (nRCP in c(min.nRCP:max.nRCP)){
  for (form in c(1:length(allForms))){
    cat(c("\n\nRCP",nRCP,"of",max.nRCP, "/// formula",form, "of", length(allForms), "\n\n"))
    iteration<-regimix.multifit(form.RCP=allForms[form][[1]], data=rcp_data, nRCP=nRCP,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
        formula[[form]] <- sapply(iteration, function(y) y$BIC)
    names(formula)[[form]]<-paste0("formula",form)
  }
  
  nRCP_list[[nRCP]]<-formula
  names(nRCP_list)[[nRCP]]<-paste0("nRCP",nRCP)
}

# for (i in 1:length(nRCP_list)){
#   nRCP_list[[i]][sapply(nRCP_list[i], is.null)] <- NULL
# }


for (i in 1:length(nRCP_list)){
nRCP_list[[i]]<-sapply(nRCP_list[[i]], function(x) min(x))
}

nRCP_list[[1]]<-NULL

res<-data.frame(minBIC=unlist(nRCP_list))
res$iteration<-row.names(res)
res$nRCP=gsub("\\..*","",gsub("nRCP","",row.names(res)))
res$formula=gsub("nRCP[[:digit:]]*\\.","",gsub("formula","",row.names(res)))

plot(res$nRCP, res$minBIC, pch=res$formula, col=res$formula, xlab="n RCPs",  ylab="min BIC for formula")

png("./LoopModelSelection_Results_run2/Formula_BICs.png", res=200, height=10,width=10, units="cm")
plot(res$nRCP, res$minBIC, pch=res$formula, col=res$formula, xlab="n RCPs",  ylab="min BIC for formula")
dev.off()

res<-res[order(res$minBIC),]
head(res)

best <- res[which.min(res$minBIC),]

#This gives us the minimum BIC
bestForm<-allForms[as.numeric((best$formula))]
bestnRCP<-as.numeric(best$nRCP)

saveRDS(nRCP_list, "./LoopModelSelection_Results_run2/nRCP_list.rds")
saveRDS(res, "./LoopModelSelection_Results_run2/res.rds")

res2<-res
res2$terms<-as.character(res2$terms)
write.csv(res2,"./LoopModelSelection_Results_run2/results_600sites_2,3,7,8vars_6-9RCPs_2018.09.14.csv")

# End run on alien
#####

#nRCP_list<-readRDS("nRCP_list.rds")
#res<-readRDS("res.rds")


res$form<-rep(options$form2)
rownames(res)<-NULL
res<-res[,c(2:5,1)]
res<-res[order(res$minBIC),]

res$nTerms<-stringr::str_count(res$form, "\\+")+1
res$terms<-strsplit(res$form, "\\s\\+\\s")

#get contribution of each individual variable
drop1<-res[res$nTerms %in% c(max(res$nTerms), max(res$nTerms)-1),]
drop1<-drop1[order(drop1$form),]
drop1

table<-list()
for (n in c(min(as.numeric(res$nRCP)):max(as.numeric(res$nRCP)))){
subsetRCP<-drop1[drop1$nRCP==n,]
start<-subsetRCP[subsetRCP$nTerms==max(res$nTerms),]

subsetRCP$whichMissing[subsetRCP$nTerms==max(res$nTerms)]<-0
for (i in which(subsetRCP$nTerms!=max(res$nTerms))){
  subsetRCP$whichMissing[i]<-    start$terms[[1]][!(start$terms[[1]]%in% subsetRCP$terms[[i]]  )]
  }

subsetRCP$deltaBIC<-subsetRCP$minBIC[subsetRCP$nTerms==max(res$nTerms)]-subsetRCP$minBIC
table[[n]]<-data.frame(subsetRCP[,names(subsetRCP)%in% c("whichMissing","nRCP","deltaBIC")])
}

#Great - this gets us the "effect" of each variable. 
table<-do.call("rbind", table)

png("./LoopModelSelection_Results_run2//VariableContribution.png", res=200, height=25,width=50, units="cm", pointsize=10)
plot(as.factor(table$whichMissing[table$whichMissing!=0]), table$deltaBIC[table$whichMissing!=0], ylab=c("Change in BIC from full model,\naveraged across nRCP iterations"), col="lightgray", )
abline(h=0, col="red")
dev.off()

png("./LoopModelSelection_Results_run2//VariableContribution_byRCPRun.png", res=200, height=25,width=50, units="cm", pointsize=10)
plot(as.factor(table$whichMissing[table$whichMissing!=0&table$nRCP==3]), table$deltaBIC[table$whichMissing!=0&table$nRCP==3], border="purple", ylab=c("Change in BIC from full model,\naveraged across nRCP iterations"), ylim=c(min(table$deltaBIC), max(table$deltaBIC)))
plot(as.factor(table$whichMissing[table$whichMissing!=0&table$nRCP==4]), table$deltaBIC[table$whichMissing!=0&table$nRCP==4], border="blue", add=T)
plot(as.factor(table$whichMissing[table$whichMissing!=0&table$nRCP==5]), table$deltaBIC[table$whichMissing!=0&table$nRCP==5], border="darkgreen", add=T)
abline(h=0, col="red")
dev.off()

#Then run full model 




