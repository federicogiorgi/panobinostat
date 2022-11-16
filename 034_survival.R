source("../shared/functions/stepsurvival_cox.R")



cands<-c("BAZ1A","HCFC1","ZNF146")



root<-"kocak_NBL"
load(paste0("D:/Dropbox/projects/masterset/data/",root,"-survival.rda"))
load(paste0("D:/Dropbox/projects/masterset/data/",root,"-expmat.rda"))

# Prepare combined signature
metasig<-expmat[cands,]
metasig<-t(apply(metasig,1,rank))
metasig<-apply(metasig,2,mean)

# Prepare dataset
common<-intersect(names(metasig),rownames(survival))
metasig<-metasig[common]
survival<-survival[common,]

# Fit model
png(paste0("plots/034_survival.png"),width=4000,height=4000,res=600)
plotstepsurv(metasig,survival,ngroups=8,ylab="BAZ1A+HCHFC1+ZNF146 mean rank expression",title="Kocak NBL dataset",plot=TRUE,legend=FALSE,cox_survival=TRUE)
dev.off()

