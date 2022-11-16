cands<-c("BAZ1A","HCFC1","HDGF","MAZ","ZNF146","ZNF532")
more<-c(
    #"ABCC3",
    "EZH2",
    #"GAPDH",
    "MYC",
    "MYCN"
)
cands<-c(cands,more)

source("../shared/functions/geneids.R")


# TARGET dataset
load("../pancancer/NBL/NBL-rawcounts.rda")
expmat<-apply(rawcounts,2,function(x){1E6*x/sum(x)})
rownames(expmat)<-eg2sym(rownames(expmat))
sub<-expmat[cands,]
means1<-apply(sub,1,mean)

# Brain dataset
load("../masterset/data/gtex_Brain-rpms.rda")
sub<-rpms[cands,]
means2<-apply(sub,1,mean)

### TODO png and error bars


barplot(rbind(means2,means1),beside=TRUE,col=c("cornflowerblue","salmon"),
        ylab="Reads per Million"
        )
legend(
    "top",
    legend=c("Normal Brain","Neuroblastoma"),
    col=c("cornflowerblue","salmon"),
    pch=15
)



