### Nice PCAs (check holding.png)
### TSNEs with different perplexities, as in Holding et al., 2019 (BBA review): 10, 25, 50
### With RNA counts, treatment, cell cycle phase
### Then regress out cell cycle phases and RNA counts

source("../shared/functions/qol.R")
source("../shared/functions/geneids.R")
library(Seurat)
library(Rtsne)


### The gene counts were calculated by CellRanger vs the hg38 annotation in Kelly cells
# Hashtag-1 is DMSO
# Hashtag-2 is Panobinostat
# Treatment (Hashtag-2) was Panobinostat 10 nmol/L for 8 hours
load("data/rawcounts.rda")
nrnas<-apply(rawcounts,2,sum)

mean(apply(rawcounts,2,mean))

### Hashtags
if(TRUE){
    # Stoeckius hashtags oligos
    x<-log10(rawcounts["Hashtag-1",]+0.1)
    y<-log10(rawcounts["Hashtag-2",]+0.1)
    png("plots/011_hashtags.png",w=3000,h=3000,res=500)
    plot(x,y,xlab="Log10 Counts of Hashtag-1 (Control)",ylab="Log10 Counts of Hashtag-2 (Panobinostat)",pch=2,main="Hashtag distribution")
    dev.off()
    
    # Group assignment
    x<-rawcounts["Hashtag-1",]
    y<-rawcounts["Hashtag-2",]
    ll<-100
    ratio<-10
    hashtags<-setNames(rep(NA,ncol(rawcounts)),colnames(rawcounts))
    hashtags[(y<=ll&(x/y>=ratio))|(y==0&x>0)]<-"ctr"
    hashtags[(x<=ll&(y/x>=ratio))|(x==0&y>0)]<-"pnb"
    hashtags[is.na(hashtags)]<-"multiplet"
    pchs<-htcols<-hashtags
    htcols[hashtags=="ctr"]<-"blue"
    htcols[hashtags=="pnb"]<-"red3"
    htcols[hashtags=="multiplet"]<-"purple"
    pchs[hashtags=="ctr"]<-3
    pchs[hashtags=="pnb"]<-5
    pchs[hashtags=="multiplet"]<-9
    pchs<-as.numeric(pchs)
    x<-log10(x+0.1)
    y<-log10(y+0.1)
    
    png("plots/011_hashtags_assignment.png",w=3000,h=3000,res=500)
    plot(x,y,xlab="Log10 Counts of Hashtag-1 (Control)",ylab="Log10 Counts of Hashtag-2 (Panobinostat)",frame.plot=FALSE,
         pch=pchs,col=htcols,main="Hashtag distribution in Dataset")
    legend("bottomleft",pch=c(3,5,9),col=c("blue","red3","purple"),bg="ivory",
           legend=c(
               paste0("ctr=",table(hashtags)["ctr"]),
               paste0("pnb=",table(hashtags)["pnb"]),
               paste0("multiplet=",table(hashtags)["multiplet"])
           )
    )
    dev.off()
}


###### Seurat pipeline ?seurat
if(TRUE){
    ## We will keep all genes expressed in >=3 cells and all cells with at least 1000 detected genes:
    dim(rawcounts) # 33460 genes, 2678 samples
    seuset<-CreateSeuratObject(counts=rawcounts,project="panobinostat",min.cells=3,min.features=1000)
    seuset # 18557 genes, 2389 samples
    
    ## Normalize with Seurat
    # By default, Seurat employs a global-scaling normalization method LogNormalize that
    # normalizes the gene expression measurements for each cell by the total expression,
    # multiplies this by a scale factor (10,000 by default), and log-transforms the result:
    seuset<-NormalizeData(seuset,normalization.method="LogNormalize",scale.factor=10000)
    
    ## Mean variability plot showing most expressed and variable genes
    expmat<-as.matrix(seuset[["RNA"]]@data)
    dim(expmat) # 18564 genes, 2659 samples
    avgexp<-apply(expmat,1,ExpMean)
    dispersion<-apply(expmat,1,LogVMR) # variance of mean ratio
    png(paste0("plots/011_meanvarplot.png"),w=4000,h=3000,res=500)
    plot(avgexp,dispersion,xlab="Average LogNorm Gene Expression",ylab="Dispersion",pch=20,main="Gene Expression Distribution",
         col="salmon",xlim=c(0,max(avgexp)*1.1))
    toshow<-c("MYCN",names(sort(rank(dispersion)+rank(avgexp),dec=TRUE))[1:30])
    textplot3(avgexp[toshow],dispersion[toshow],words=toshow,font=2,cex=0.9,line.col="#00000033",padding="  ")
    grid()
    dev.off()
}


#### TSNE 1: simple normalization
if(TRUE){
    # Get Highly variable genes for clustering
    seuset<-FindVariableFeatures(seuset,mean.function=ExpMean,dispersion.function=LogVMR,nfeatures=5000)
    hvargenes<-seuset[["RNA"]]@var.features
    # TSNE based on these genes
    tsnemat<-expmat[hvargenes,]
    fname<-"results/011a_rtsne.rda"
    if(!file.exists(fname)){
        set.seed(1)
        ttt<-Rtsne(t(tsnemat))
        save(ttt,file=fname)
    }else{load(fname)}
    x<-ttt$Y[,1]
    y<-ttt$Y[,2]
    
    ### Plot TSNE with RNA counts, hashtag and cell cycle
    ### nrnas as color
    here_nrnas<-nrnas[colnames(tsnemat)]
    colfunc<-colorRampPalette(c("blue","green","red"))
    readcols<-setNames(colfunc(100)[as.numeric(cut(here_nrnas,breaks=100))],names(here_nrnas))
    # readcols<-readcols[colnames(tsnemat)]
    
    png("plots/011a_rtsneA_nrnas.png",w=3500,h=3000,res=600)
    layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
    par(mar=c(4,4,3,0))
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="RNA counts",col=readcols,pch=20)
    par(mar=c(0,0,0,2))
    plot(c(0,2),c(0,1),type='n',axes=F,xlab='',ylab='')
    legend_image<-as.raster(rev(matrix(colfunc(100), ncol=1)))
    rasterImage(legend_image,0,0,1,1)
    text(x=1.5,y=seq(0,1,l=5),labels=paste0(round(quantile(here_nrnas/1E3)),"k"))
    dev.off()
    
    ### Hashtags as colors
    htcols<-hashtags
    htcols[hashtags=="ctr"]<-"blue"
    htcols[hashtags=="pnb"]<-"red3"
    htcols[hashtags=="multiplet"]<-"purple"
    htcols<-htcols[colnames(tsnemat)]
    png("plots/011a_rtsneB_hashtags.png",w=3000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Hashtags",col=htcols,pch=20)
    legend("topright",pch=20,col=c("blue","red3","purple"),bg="ivory",pt.cex=2,cex=0.7,
           legend=c(
               paste0("ctr=",table(hashtags[colnames(tsnemat)])["ctr"]),
               paste0("pnb=",table(hashtags[colnames(tsnemat)])["pnb"]),
               paste0("multiplet=",table(hashtags[colnames(tsnemat)])["multiplet"])
           )
    )
    dev.off()
    
    
    ### Cell cycle phases as colors
    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    s.genes<-cc.genes$s.genes
    g2m.genes<-cc.genes$g2m.genes
    ## Assign Cell-Cycle Scores
    seuset<-CellCycleScoring(seuset,s.features=s.genes,g2m.features=g2m.genes,set.ident=TRUE)
    # View cell cycle scores and phase assignments
    head(seuset@meta.data)
    # Use Seurat assignment
    phases<-cyclecols<-setNames(as.character(seuset@meta.data$Phase),rownames(seuset@meta.data))
    cyclecols[phases=="G1"]<-"salmon"
    cyclecols[phases=="G2M"]<-"seagreen"
    cyclecols[phases=="S"]<-"cornflowerblue"
    cyclecols<-cyclecols[colnames(tsnemat)]

    png("plots/011a_rtsneC_CellCycle.png",w=3000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Cell Cycle",col=cyclecols,pch=20)
    legend("topright",pch=19,col=c("salmon","seagreen","cornflowerblue"),legend=c("G1","G2M","S"),bg="ivory",title="Cell Cycle Phase",cex=0.7)
    dev.off()
}

### Some numbers from the initial assignment
if(TRUE){
    library(gridExtra)
    library(grid)
    
    hashtags<-hashtags[names(phases)]
    distro<-table(hashtags,phases)
    png("plots/011_cellCycle_cellNumbers.png",w=3000,h=2000,res=500)
    grid.table(distro)
    dev.off()
    #           phases
    # hashtags     G1 G2M   S
    # ctr       344 340 671
    # multiplet  14  30  52
    # pnb       447 151 340
    tmp<-apply(distro,1,sum)
    distroperc<-as.matrix(distro)
    for(i in 1:nrow(distroperc)){
        distroperc[i,]<-paste0(signif(100*as.numeric(distroperc[i,])/tmp[i],4),"%")
    }
    png("plots/011_cellCycle_cellPercentages.png",w=3000,h=2000,res=500)
    grid.table(distroperc)
    dev.off()
    #               phases
    # hashtags    G1     G2M    S     
    # ctr       25.39% 25.09% 49.52%
    # multiplet 14.58% 31.25% 54.17%
    # pnb       47.65% 16.1%  36.25%
}


### Regress out cell cycle and RNA counts
if(TRUE){
    fname<-"results/011_seuset_regressed.rda"
    if(!file.exists(fname)){
        # Remove multiplet cells and cells with low nrnas
        cells<-names(nrnas)[nrnas>=10000]
        cells<-intersect(names(hashtags)[hashtags%in%c("pnb","ctr")],cells)
        pnb<-intersect(cells,names(hashtags)[hashtags=="pnb"])
        ctr<-intersect(cells,names(hashtags)[hashtags=="ctr"])
        length(pnb) # 905
        length(ctr) # 1278
        subcounts<-rawcounts[,c(pnb,ctr)]
        
        # Repeat previous steps just to be sure
        seuset<-CreateSeuratObject(counts=subcounts,project="panobinostat",min.cells=3,min.features=1000)
        seuset<-NormalizeData(seuset,normalization.method="LogNormalize",scale.factor=10000)
        seuset<-FindVariableFeatures(seuset,mean.function=ExpMean,dispersion.function=LogVMR,nfeatures=5000)
        seuset<-CellCycleScoring(seuset,s.features=s.genes,g2m.features=g2m.genes,set.ident=TRUE)
        seuset<-ScaleData(object=seuset,vars.to.regress=c("nCount_RNA","S.Score","G2M.Score"),features=rownames(seuset))
        save(seuset,file=fname)
    } else {load(fname)}
    regexpmat<-as.matrix(seuset[["RNA"]]@scale.data)
    dim(regexpmat) # 18481 2183
}


#### TSNE 2: post regression
if(TRUE){
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
    x<-ttt$Y[,1]
    y<-ttt$Y[,2]
    
    ### Plot TSNE with RNA counts, hashtag and cell cycle
    ### nrnas as color
    nrnas<-here_nrnas[colnames(tsnemat)]
    colfunc<-colorRampPalette(c("blue","green","red"))
    readcols<-setNames(colfunc(100)[as.numeric(cut(here_nrnas,breaks=100))],names(here_nrnas))
    # readcols<-readcols[colnames(tsnemat)]
    
    png("plots/011b_rtsneA_nrnas.png",w=4000,h=3000,res=500)
    layout(matrix(1:2,ncol=2), width = c(3,1),height = c(1,1))
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="RNA counts in Cell Clustering",col=readcols,pch=20)
    plot(c(0,2),c(0,1),type='n',axes=F,xlab='',ylab='',main='RNA counts')
    legend_image<-as.raster(rev(matrix(colfunc(100), ncol=1)))
    rasterImage(legend_image,0,0,1,1)
    text(x=1.5,y=seq(0,1,l=5),labels=paste0(round(quantile(here_nrnas/1E3)),"k"))
    dev.off()
    
    ### Hashtags as colors
    htcols<-hashtags
    htcols[hashtags=="ctr"]<-"blue"
    htcols[hashtags=="pnb"]<-"red3"
    htcols[hashtags=="multiplet"]<-"purple"
    htcols<-htcols[colnames(tsnemat)]
    png("plots/011b_rtsneB_hashtags.png",w=3000,h=3000,res=500)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Hashtags in Cell Clustering",col=htcols,pch=20)
    legend("topleft",pch=20,col=c("blue","red3","purple"),bg="ivory",pt.cex=2,cex=0.7,
           legend=c(
               paste0("ctr=",table(hashtags[colnames(tsnemat)])["ctr"]),
               paste0("pnb=",table(hashtags[colnames(tsnemat)])["pnb"]),
               paste0("multiplet=",table(hashtags[colnames(tsnemat)])["multiplet"])
           )
    )
    dev.off()
    
    
    ### Cell cycle phases as colors
    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
    # segregate this list into markers of G2/M phase and markers of S phase
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    ## Assign Cell-Cycle Scores
    seuset<-CellCycleScoring(seuset,s.features=s.genes,g2m.features=g2m.genes,set.ident=TRUE)
    # View cell cycle scores and phase assignments
    head(seuset@meta.data)
    # Use Seurat assignment
    phases<-cyclecols<-setNames(as.character(seuset@meta.data$Phase),rownames(seuset@meta.data))
    cyclecols[phases=="G1"]<-"salmon"
    cyclecols[phases=="G2M"]<-"seagreen"
    cyclecols[phases=="S"]<-"cornflowerblue"
    cyclecols<-cyclecols[colnames(tsnemat)]

    png("plots/011b_rtsneC_CellCycle.png",w=3000,h=3000,res=500)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Cell Cycle Clustering",col=cyclecols,pch=20)
    legend("topleft",pch=19,col=c("salmon","seagreen","cornflowerblue"),legend=c("G1","G2M","S"),bg="ivory",title="Cell Cycle Phase",cex=0.7)
    dev.off()
}






### Based on the simple normalized data we will calculate global differential expression
save(hashtags,phases,nrnas,rawcounts,expmat,regexpmat,file="data/seurat.rda")







