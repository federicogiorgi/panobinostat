# cd /mnt/d/Dropbox/projects/panobinostat
# tar cvzf viper.tar.gz viper
# remove.packages("viper")
# install.packages("viper.tar.gz",repos=NULL)
# Switched aecdf1 for aecdf somehwere in the viper code
# Fixed error in aecdf1 by using the aecdf (aecdf1 is used for large datasets, but it contains a bug)
library(viper)

# Create treatment and control expression matrices
source("../shared/functions/qol.R")
load("data/seurat.rda")
dim(regexpmat) # 18481  2183
expmat<-regexpmat
rm(rawcounts,regexpmat)
hashtags<-hashtags[colnames(expmat)]


# Regulons (generated and available from the TEAD4 paper http://cancerdiscovery.aacrjournals.org/content/early/2018/03/06/2159-8290.CD-16-0861)
load("fromGonzalo/l_viper_nets.rda") # the Target regulon
load("fromGonzalo/n_viper_nets.rda") # the NRC regulon


# # Add white noise
# sd(expmat) # 0.9190246
# set.seed(1)
# expmat<-expmat+abs(rnorm(ncol(expmat)*nrow(expmat),sd=0.01))

# Generate VIPER signature
fname<-"data/vipersignature.rda"
if(!file.exists(fname)){
    signature<-viperSignature(expmat, expmat, method="mean")
    save(signature,file=fname)
} else {load(fname)}


### Generate VIPER mats
# Error in aecdf1(tmp[i, ], symmetric = TRUE, es[i, ]) : 
# object 'tl2' not found
fname<-"data/vipermats.rda"
if(!file.exists(fname)){
    tvipermat<-viper(signature,regulon=l_vnet, pleiotropy=FALSE, minsize=10) # TARGET viper
    nvipermat<-viper(signature,regulon=n_vnet, pleiotropy=FALSE, minsize=10) # NRC viper
    save(tvipermat,nvipermat,file=fname)
} else {load(fname)}

# Rtsne
library(Rtsne)
fname<-"results/017_rtsne.rda"
if(!file.exists(fname)){
    set.seed(1)
    tartsne<-Rtsne(t(tvipermat))
    nrctsne<-Rtsne(t(nvipermat))
    save(tartsne,nrctsne,file=fname)
}else{load(fname)}

# We start with TARGET
x<-setNames(tartsne$Y[,1],colnames(tvipermat))
y<-setNames(tartsne$Y[,2],colnames(tvipermat))


##################
# PLOTTING
mytfs<-c("HCFC1","BAZ1A","MAZ","ZNF146", # 4 tested TFs
         "MYCN","STAT3","TEAD2","SMARCA4","TOP2B" # Other interesting TFs
)


### Smarter way to pick it (higher difference, higher intra-group variance)
tctr<-tvipermat[,names(hashtags)[hashtags=="ctr"]]
tpnb<-tvipermat[,names(hashtags)[hashtags=="pnb"]]
nctr<-nvipermat[,names(hashtags)[hashtags=="ctr"]]
npnb<-nvipermat[,names(hashtags)[hashtags=="pnb"]]

var_tctr<-apply(tctr,1,var)
var_tpnb<-apply(tpnb,1,var)
var_tdiff<-sort(var_tpnb-var_tctr,dec=TRUE)
var_nctr<-apply(nctr,1,var)
var_npnb<-apply(npnb,1,var)
var_ndiff<-sort(var_npnb-var_nctr,dec=TRUE)

# TFs whose variance is higher in panobinostat than in control
tvar<-names(var_tdiff)[1:5]
nvar<-names(var_ndiff)[1:5]

png("plots/017_allTFs_variance.png",w=6000,h=3000,res=600)
par(mfrow=c(1,2))
plot(var_tctr,var_tpnb,main="TF Activity Variance",pch=20,col="#00000011",
     xlab="Variance in Control",ylab="Variance in Panobinostat",xlim=c(-0.2,2)
)
mtext("TARGET network",cex=0.8)
textplot3(var_tctr[c(mytfs,tvar)],var_tpnb[c(mytfs,tvar)],words=c(mytfs,tvar),col="navy")
plot(var_nctr,var_npnb,main="TF Activity Variance",pch=20,col="#00000011",
     xlab="Variance in Control",ylab="Variance in Panobinostat",xlim=c(-0.2,2)
)
mtext("NRC network",cex=0.8)
textplot3(var_nctr[c(mytfs,nvar)],var_npnb[c(mytfs,nvar)],words=c(mytfs,nvar),col="navy")
dev.off()


# Simpler plot
png("plots/017_allTFs_variance_target.png",w=3000,h=3000,res=600)
plot(var_tctr,var_tpnb,main="TF Activity Variance",pch=20,col="#00000011",
     xlab="Variance in Control",ylab="Variance in Panobinostat",xlim=c(-0.2,1.2),ylim=c(-0.2,1.2)
)
mytfs<-c("HCFC1","BAZ1A","MAZ","ZNF146","MYCN")
textplot3(var_tctr[c(mytfs,tvar)],var_tpnb[c(mytfs,tvar)],words=c(mytfs,tvar),col="navy",font=2,cex=0.8)
lm1<-lm(var_tpnb~var_tctr)
abline(lm1$coef,lty=2,col="red3") # Fit
abline(0,1) # Bisector
legend("topleft",lty=c(1,2),col=c("black","red3"),legend=c("Bisector","Line of Best Fit"))
dev.off()


#### Pick TF
moretfs<-c(mytfs,tvar)
for(tf in mytfs){
    
    # TF nes to be plotted
    nes<-tvipermat[tf,]
    
    ## Colors
    # Turn this into a function (input: nes, output: readcols)
    colfunc<-colorRampPalette(c("#00008000","#00008099","#000080CC","#EEEEEECC","#810000CC","#81000099","#81000000"))
    z<-nes
    nbreaks<-1000
    extreme = round(max(abs(z)))
    breaks<-seq(-extreme, extreme, length = nbreaks)
    z<-z-mean(z)
    ncol <- length(breaks) - 1
    col<-colfunc(ncol)
    CUT <- cut(z, breaks = breaks)
    colorlevels <- col[match(CUT, levels(CUT))]
    names(colorlevels) <- names(z)
    readcols<-colorlevels
    
    # Symbols
    pch<-hashtags
    pch[hashtags=="ctr"]<-20
    pch[hashtags=="pnb"]<-15
    pch<-as.numeric(pch)
    
    # Plot
    png(paste0("plots/017_viper_target_",tf,".png"),w=4200,h=3000,res=500,family="xkcd")
    layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main=paste0(tf," activity in dataset"),type="n")
    points(x,y,col=readcols,pch=pch,cex=0.8)
    legend("topright",pch=c(20,15),legend=c("ctr","pnb"))
    plot(c(0,2),c(0,1),type='n',axes=F,xlab='',ylab="Normalized Enrichment Score",main=paste0(tf," Activity"))
    legend_image<-as.raster(rev(matrix(col, ncol=1)))
    rasterImage(legend_image,0,0,1,1)
    text(x=1.5,y=seq(0,1,l=length(pretty(breaks))),labels=pretty(breaks))
    dev.off()
}


