###### Single-cell GSEA
source("../shared/functions/qol.R")
library(fgsea)

load("data/seurat.rda")



## Load the msigdb database 7.0.1
fname<-"results/015_mlist.rda"
if(!file.exists(fname)){
  library(msigdbr)
  msigdbr_show_species()
  mdf<-msigdbr(species="Homo sapiens") # Retrieve all human gene sets
  length(mdf$gs_name) # Number of associations (2639029)
  length(unique(mdf$gs_name)) # Number of pathways (22596)
  head(mdf)
  mlist<-mdf %>% split(x=.$gene_symbol,f=.$gs_name)
  save(mlist,file=fname)
}else{load(fname)}
grep("(REPAIR|INTERFERON)",names(mlist),value="TRUE")
grep("NEURON_DIFFERENTIATION",names(mlist),value="TRUE")
plength<-sapply(mlist,length)
max(plength) # 2775, the biggest pathway



#### Calculate single-cell gsea (better to use area)
source("../shared/functions/area.R")

fname<-"results/015_ssgsea_regexpmat.rda"
if(!file.exists(fname)){
    scalemat<-t(scale(t(regexpmat)))
    ssgsea<-area(signatures=scalemat,groups=mlist,minsize=4)
    save(ssgsea,file=fname)    
}else{load(fname)}


##############################
#### Overlay pathways on cells
library(Seurat)
load("results/011_seuset_regressed.rda")
# Get Highly variable genes for clustering
seuset<-FindVariableFeatures(seuset,mean.function=ExpMean,dispersion.function=LogVMR,nfeatures=5000)
hvargenes<-seuset[["RNA"]]@var.features
# TSNE based on these genes
tsnemat<-regexpmat[hvargenes,]
fname<-"results/011b_rtsne.rda"
if(!file.exists(fname)){
    set.seed(1)
    ttt<-Rtsne(t(tsnemat))
    save(ttt,file=fname)
}else{load(fname)}
x<-setNames(ttt$Y[,1],colnames(tsnemat))
y<-setNames(ttt$Y[,2],colnames(tsnemat))

dim(ssgsea) # 21941  2183
length(x) # 2183

### Plot TSNE with RNA counts, hashtag and cell cycle
# Pick a path and color it
paths<-grep("(APOPTOSIS|SENESCENCE|MESENCHYMAL|EMT|DIFFERENTIATION|PROLIFERATION)",names(mlist),value="TRUE")
paths<-c("WONG_EMBRYONIC_STEM_CELL_CORE","HELLER_HDAC_TARGETS_UP",paths)

# Shape according to origin
pchs<-hashtags[names(x)]
pchs[pchs=="ctr"]<-16
pchs[pchs=="pnb"]<-15
pchs<-as.numeric(pchs)

dir.create("plots/015_ssgsea_paths")
for(path in paths){
    if(path %in%rownames(ssgsea)){ # Some are too small and were not calculated
        # Format path name
        nicepath<-gsub("_"," ",path)
        nicepath<-tools::toTitleCase(tolower(nicepath))
        nes<-ssgsea[path,colnames(tsnemat)]
        colfunc<-colorRampPalette(c("#00009999","gray95","#99000099"))
        readcols<-setNames(colfunc(100)[as.numeric(cut(nes,breaks=100))],names(nes))
        # Plot
        png(paste0("plots/015_ssgsea_paths/015_rtsne_",path,".png"),w=4000,h=3000,res=500)
        layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
        plot(x,y,xlab="TSNE1",ylab="TSNE2",main=nicepath,type="n")
        points(x,y,col=readcols,pch=pchs)
        legend("topright",pch=c(16,15),legend=c("Ctrl","Pnb"),col="grey")
        plot(c(0,2),c(0,1),type='n',axes=F,xlab='',ylab="Normalized Enrichment Score",main="Pathway Score")
        legend_image<-as.raster(rev(matrix(colfunc(100), ncol=1)))
        rasterImage(legend_image,0,0,1,1)
        #text(x=1.5,y=seq(0,1,l=5),labels=paste0(signif(quantile(nes),2)))
        fourpoints<-signif(c(-max(abs(nes)),-max(abs(nes))/2,0,max(abs(nes))/2,max(abs(nes))),2)
        text(x=1.5,y=seq(0,1,l=5),labels=fourpoints)
        dev.off()
    }
}




###### Use heatmaps to vary
# Prepare the matrix
load("results/014_topPathways.rda")
nicetop<-gsub("_"," ",top)
nicetop<-tools::toTitleCase(tolower(nicetop))
plotmat<-ssgsea[top,]
rownames(plotmat)<-nicetop

# # Scale it row-wise
# plotmat<-t(scale(t(plotmat)))
# Set a maximum
plotmat[plotmat>10]<-10
plotmat[plotmat<(-10)]<-(-10)

# Prepare the annotation
hashtags<-hashtags[colnames(plotmat)]
#plotmat<-plotmat[,order(hashtags)] # Order (optional)
colannot<-setNames(rep("gnegne",ncol(plotmat)),colnames(plotmat))
colannot[names(hashtags)[hashtags=="pnb"]]<-"black"
colannot[names(hashtags)[hashtags=="ctr"]]<-"white"
colannot<-colannot[colnames(plotmat)]
colannot<-t(t(colannot))
colnames(colannot)<-"Panobinostat (in black)"
rownames(plotmat)<-gsub("Mrna","mRNA",rownames(plotmat))

# Plot it
source("../shared/functions/heatmaps.R")
colfun<-colorRampPalette(c("navy","navy","white","red3","red3"))

png("plots/015_heatmap_topGSEAs.png",w=6000,h=3000,res=600)
heatmap.3(plotmat,ColSideColors=colannot,KeyValueName="Scaled NES",mar=c(0,16),
          col=colfun,cexRow= 0.6,Colv=FALSE,dendrogram="row",scale="row")
dev.off()




# Show the pathways most varying within PNB group


# Show the pathways most varying within PNB group and not varying in CTR group







