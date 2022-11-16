
### Classic Marina using both Gonzalo's networks
library(viper)
options(java.parameters = "-Xmx10000m");library(xlsx)

# Create treatment and control expression matrices
load("data/seurat.rda")
dim(expmat) # 18557  2389
expmat<-expmat[rownames(regexpmat),colnames(regexpmat)]
dim(expmat) # 18481  2183
hashtags<-hashtags[colnames(expmat)]

# Regulons (generated and available from the TEAD4 paper http://cancerdiscovery.aacrjournals.org/content/early/2018/03/06/2159-8290.CD-16-0861)
load("fromGonzalo/l_viper_nets.rda") # the Target regulon
load("fromGonzalo/n_viper_nets.rda") # the NRC regulon

### MRA 1: expression matrix pnb vs. ctrl, no regression ----
if(FALSE){
    trt<-expmat[,names(hashtags)[hashtags=="pnb"]]
    ctr<-expmat[,names(hashtags)[hashtags=="ctr"]]
    fname<-"results/016_mra1_original.rda"
    if(!file.exists(fname)){
        # Calculate basic signature treatment vs. control
        sig<-rowTtest(trt,ctr)$statistic[,1]
        # Sanity check: remove eventually NaN and NA
        sig<-sig[!is.nan(sig)]
        sig<-sig[!is.na(sig)]
        # Calculate a robust signature using permutation of dataset
        dnull<-ttestNull(trt,ctr,per=100) # Only 100 permutations to reduce computation time, but it is recommended to perform at least 1000 permutations
        
        # Run master regulator analysis (target dataset)
        tmra<-msviper(sig,regulon=l_vnet,dnull)
        
        # Run master regulator analysis (NRC dataset)
        nmra<-msviper(sig,regulon=n_vnet,dnull)
        
        save(tmra,nmra,file=fname)
    }else{load(fname)}
    
    png("plots/016_mra1_target.png",family="xkcd",w=2000,h=4000,res=400)
    plot(tmra, cex=.7,mrs=50)
    dev.off()
    tnes<-tmra$es$nes
    
    png("plots/016_mra1_nrc.png",family="xkcd",w=2000,h=4000,res=400)
    plot(nmra, cex=.7,mrs=50)
    dev.off()
    nnes<-nmra$es$nes
    
    # Compare two mras (master regulator analyses)
    common<-intersect(names(tnes),names(nnes))
    tnes<-tnes[common]
    nnes<-nnes[common]
    
    tops<-names(sort(abs(tnes)+abs(nnes),dec=TRUE)[1:20])
    source("../shared/functions/qol.R")
    
    png("plots/016_mra1_comparison.png",family="xkcd",w=3000,h=3000,res=450)
    x<-tnes
    y<-nnes
    plot(x,y,pch=20,main="Master Regulator Analysis",xlab="TARGET regulon",ylab="NRC regulon",
         col="#00000011",xlim=1.4*c(min(x),max(x)),ylim=1.4*c(min(y),max(y))
    )
    scc<-signif(cor(x,y,method="s"),4)
    mtext(paste0("Panobinostat vs. Ctrl Signature, SCC=",scc))
    textplot3(x[tops],y[tops],words=tops,col="black")
    dev.off()
}

### MRA 2: expression matrix pnb vs. ctrl, with regression (the one used in TABELLONE) ----
if(TRUE){
    trt<-regexpmat[,names(hashtags)[hashtags=="pnb"]]
    ctr<-regexpmat[,names(hashtags)[hashtags=="ctr"]]
    fname<-"results/016_mra2_regressed.rda"
    if(!file.exists(fname)){
        # Calculate basic signature treatment vs. control
        sig<-rowTtest(trt,ctr)$statistic[,1]
        # Sanity check: remove eventually NaN and NA
        sig<-sig[!is.nan(sig)]
        sig<-sig[!is.na(sig)]
        # Calculate a robust signature using permutation of dataset
        dnull<-ttestNull(trt,ctr,per=100) # Only 100 permutations to reduce computation time, but it is recommended to perform at least 1000 permutations
        
        # Run master regulator analysis (target dataset)
        tmra<-msviper(sig,regulon=l_vnet,dnull)
        
        # Run master regulator analysis (NRC dataset)
        nmra<-msviper(sig,regulon=n_vnet,dnull)
        
        save(tmra,nmra,file=fname)
    }else{load(fname)}
    
    toadd<-c("BAZ1A","HCFC1","MAZ","MYCN","ZNF146")
    
    png("plots/016_mra2_target.png",w=4000,h=4000,res=600)
    toshow<-names(sort(abs(tmra$es$nes),dec=TRUE))[1:16]
    toshow<-sort(unique(c(toshow,toadd)))
    plot(tmra, cex=.7,mrs=toshow)
    dev.off()
    tnes<-tmra$es$nes
    
    png("plots/016_mra2_nrc.png",w=4000,h=4000,res=600)
    toshow<-names(sort(abs(nmra$es$nes),dec=TRUE))[1:16]
    toshow<-sort(unique(c(toshow,toadd)))
    plot(nmra, cex=.7,mrs=toshow)
    dev.off()
    nnes<-nmra$es$nes
    
    # Compare two mras (master regulator analyses)
    common<-intersect(names(tnes),names(nnes))
    tnes<-tnes[common]
    nnes<-nnes[common]
    
    tops<-names(sort(abs(tnes)+abs(nnes),dec=TRUE)[1:20])
    tops<-unique(c(tops,toadd))
    source("../shared/functions/qol.R")
    
    png("plots/016_mra2_comparison.png",w=3000,h=3000,res=500)
    x<-tnes
    y<-nnes
    plot(x,y,pch=20,main="Master Regulator Analysis",xlab="TARGET network",ylab="NRC network",
         col="#00000022",xlim=1.4*c(min(x),max(x)),ylim=1.4*c(min(y),max(y))
    )
    scc<-signif(cor(x,y,method="s"),4)
    mtext(paste0("Panobinostat vs. Control Signature, SCC=",scc))
    textplot3(x[tops],y[tops],words=tops,col="black")
    dev.off()
    

}


