#sample dataset
library(combinat)
library(RCPmod)
setwd("D:/Documents/Projects/WorldClass_DiveSurveys/RCP//Models/2018-08-13")
rcp_data<-readRDS("rcp_data_subsetted.rds")


id_vars="SourceKey" # Unique identifier for each record
sample_vars="Survey" # Field or fields that describe potential sampling effects
nstarts <- 10  # Editted this line from nstarts <- 1000
#max.nRCP <- 6   # Not used when manually testing number of RCPs
distribution <- "Bernoulli"  # Change to "Bernoulli" for P/A data or NegBin for abund
gen.start.val <- "random2"

species<-names(rcp_data)[3:40]
#env<-names(rcp_data)[c(41:58)]
env<-names(rcp_data)[c(41,43,45,47)]

look<-list()
for(i in 2:length(env)){
look[[i]]<-list(combn(env, i, fun=NULL, simplify=F))
}

list2<-list()
for (i in 2:length(env)){
list2[[i]]<-apply(do.call("rbind",look[[i]][[1]]),1,paste,collapse="+")}

iterations<-unlist(list2)
allForms<-list()

for (n in 1:length(iterations)){
  allForms[[n]] <- as.formula(paste("cbind(",paste(species, collapse=", "),")~",c(paste(iterations[n], collapse="+"))))
}

length(allForms) #18 env is 262143 model forms. 9 is 511. 4 is 15.

max.nRCP<-10

#AllFormula table
options<-data.frame(form=as.character(paste(allForms)))
for (x in 1:nrow(options)){ options$form2[x]<-strsplit(as.character(options$form[x]), "\\~\\s")[[1]][2]}


nRCP_list<-list()      
formula<-list()

for (nRCP in c(3:max.nRCP)){
  for (form in c(1:length(allForms))){
    print(c("RCP ",nRCP, ", formula ",form))
    iteration<-regimix.multifit(form.RCP=allForms[form][[1]], data=rcp_data, nRCP=nRCP,inits=gen.start.val, nstart=nstarts, dist=distribution, mc.cores=1)
        formula[[form]] <- sapply(iteration, function(y) y$BIC)
    names(formula)[[form]]<-paste0("formula",form)
  }
  
  nRCP_list[[nRCP]]<-formula
  names(nRCP_list)[[nRCP]]<-paste0("nRCP",nRCP)
}

for (i in 1:length(nRCP_list)){
  nRCP_list[[i]][sapply(nRCP_list[[i]], is.null)] <- NULL
}

for (i in 1:length(nRCP_list)){
nRCP_list[[i]]<-sapply(nRCP_list[[i]], function(x) min(x))
}

res<-data.frame(minBIC=unlist(nRCP_list))
res$iteration<-row.names(res)
res$nRCP=gsub("\\..*","",gsub("nRCP","",row.names(res)))
res$formula=gsub("nRCP[[:digit:]]*\\.","",gsub("formula","",row.names(res)))

plot(res$nRCP, res$minBIC, pch=res$formula, col=res$formula)

best <- res[which.min(res$minBIC),]

#This gives us the minimum BIC
bestForm<-allForms[as.numeric((best$formula))]
bestnRCP<-as.numeric(best$nRCP)

saveRDS(nRCP_list, "../Models/2018-08-13/nRCP_list.rds")
saveRDS(res, "../Models/2018-08-13/res.rds")

nRCP_list<-readRDS("nRCP_list.rds")
res<-readRDS("res.rds")

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

plot(as.factor(table$whichMissing[table$whichMissing!=0]), table$deltaBIC[table$whichMissing!=0], ylab=c("Change in BIC from full model, averaged across nRCP iterations"), col="lightgray")


#Then run full model 



# Are any RCPs consisting of a small number of sites?  (A posteriori) If so remove.
RCPsamp_minPosteriorSites <- cbind( n.sites, sapply( nRCPs_samp[-1], function(y) sapply( y, function(x) min( colSums( x$postProbs)))))
RCPsamp_ObviouslyBad <- RCPsamp_minPosteriorSites < 2
RCPsamp_BICs[RCPsamp_ObviouslyBad] <- NA


