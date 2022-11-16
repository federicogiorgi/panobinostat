
### Pathway enrichment results panobinostat vs. control
load("results/014_gsea_giorgi.rda")

## Msigdb
load("results/015_mlist.rda")

## Pathways to show
paths<-c("FRIDMAN_SENESCENCE_UP","PEART_HDAC_PROLIFERATION_CLUSTER_DN","HELLER_HDAC_TARGETS_UP","GOBP_NEURON_DIFFERENTIATION")
mlist[paths]

## DE results
load("results/013_res_pnb_vs_ctr.rda")

## Signature
signature<-setNames(-log10(res$p)*sign(res$wstat),rownames(res))
signature<-signature[!is.na(signature)]


## Set up GSEA visualization
library(corto)
# path<-"FRIDMAN_SENESCENCE_UP"

for(path in paths){
  pval<-gsea[gsea$pathway==path,"pval"] # We take both p-value and NES from the fgsea calculation
  padj<-gsea[gsea$pathway==path,"padj"]
  nes<-gsea[gsea$pathway==path,"NES"]
  set<-mlist[[path]]
  obj<-gsea(signature,set,method="permutation")
  png(paste0("plots/035_gsea_",path,".png"),w=2000,h=2000,res=450)
  plot_gsea(obj,ext_nes=nes,ext_pvalue=padj,title=path)
  dev.off()
}


### Violin plot
library(vioplot)

## Load normalized RPM (reads per million mapped reads)
what<-load("data/seurat.rda")
what<-load("data/rpms.rda")

pnb<-apply(rpms[,names(hashtags)[hashtags=="pnb"]],1,mean)
ctr<-apply(rpms[,names(hashtags)[hashtags=="ctr"]],1,mean)



png("plots/035_vioplot.png",w=4000,h=1000,res=450)
par(mfrow=c(1,4))
i<-1
for(path in paths){
  genes<-intersect(mlist[[path]],rownames(rpms))
  vioplot(log10(pnb[genes]+0.0001),log10(ctr[genes]+0.0001),col=i+1,ylab="Log10 RPM",names=c("pnb","ctr"),main=path,cex.main=0.8)
  i<-i+1
}
dev.off()


