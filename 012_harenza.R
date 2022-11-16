### This dataset from the Maris lab provides raw counts for most datasets
# https://www.nature.com/articles/sdata201733

# setwd("D:/Dropbox/projects/panobinostat")
# annotation<-read.delim("data/harenza/nbm2017annot.txt",sep="\t",as.is=TRUE)
# harenzas<-annotation[annotation$dataset=="harenza",]
# oldpath<-"H:/data/nbl/nbm2017/bams/"
# newpath<-"D:/Dropbox/projects/bigdata/nbmdata/harenza/"
# 
# for(i in 1:nrow(harenzas)){
#     srid<-harenzas$run[i]
#     cell<-harenzas$cell[i]
#     message("Processing ",cell)
#     oldbam<-paste0(oldpath,srid,".sorted.bam")
#     oldbai<-paste0(oldpath,srid,".sorted.bam.bai")
#     newbam<-paste0(newpath,"harenza_",cell,"_",srid,".sorted.bam")
#     newbai<-paste0(newpath,"harenza_",cell,"_",srid,".sorted.bam.bai")
#     file.copy(oldbam,newbam)
#     file.copy(oldbai,newbai)
# }


# Run subread
export PATH="/mnt/d/Dropbox/programs/subread-1.6.5-Linux-x86_64/bin/:$PATH"
cd /mnt/d/Dropbox/projects/bigdata/nbmdata/harenza
gtf=/mnt/d/Dropbox/projects/shared/human/hg38/annotations/gencode.v29.annotation.gtf
featureCounts -T 4 -t exon -g gene_id -a $gtf -o harenza.counts.txt *sorted.bam
cp harenza.counts.txt /mnt/d/Dropbox/projects/panobinostat/data/harenza

# Level : meta-feature level
# Paired-end : no
# Strand specific : no
# Multimapping reads : not counted
# Multi-overlapping reads : not counted
# Read orientations : fr

# Load annotation file /mnt/d/Dropbox/projects/shared/human/hg38/annotat...
# Features : 1262773
# Meta-features : 58721
# Chromosomes/contigs : 25




### Back to R -  Format Harenza dataset
fname<-"data/harenza/rawcounts.rda"
if(!file.exists(fname)){
    rawcounts<-read.delim("data/harenza/harenza.counts.txt",as.is=TRUE,skip=1,row.names=1)
    rawcounts<-as.matrix(rawcounts[,6:ncol(rawcounts)])
    # Convert to gene symbols
    source("../shared/functions/geneids.R")
    source("../shared/functions/expmerger.R")
    ensg<-rownames(rawcounts)
    ensg<-gsub("\\..+","",ensg)
    entrezs<-ens2eg(ensg)
    symbols<-eg2sym(entrezs)
    symbols[is.na(symbols)]<-ensg[is.na(symbols)]
    names(symbols)<-ensg
    length(unique(ensg)) # 58676
    length(unique(symbols)) # 58581
    # Sum reads mapping to the same symbol
    newcounts<-squish(rawcounts,symbols,method="sum",verbose=TRUE)
    rawcounts<-newcounts
    colnames(rawcounts)<-gsub("\\.sorted\\.bam","",colnames(rawcounts))
    colnames(rawcounts)<-gsub("harenza_","",colnames(rawcounts))
    colnames(rawcounts)<-gsub("_.+","",colnames(rawcounts))
    colnames(rawcounts)<-gsub("\\.","-",colnames(rawcounts))
    # RPM
    rpms<-apply(rawcounts,2,function(x){1E6*x/sum(x)})
    save(rpms,file="data/harenza/rpms.rda")
    # VST
    source("../shared/functions/vst.R")
    expmat<-vst(rawcounts)
    save(expmat,file="data/harenza/expmat.rda")
    save(rawcounts,file=fname)
} else {load(fname)}
harenza<-rawcounts
rm(rawcounts)


## Load the Giorgi dataset
load("data/seurat.rda") # hashtags,phases,nrnas,rawcounts,expmat,regexpmat



## Compare RPMs in Harenza Kelly untreated vs Giorgi Kelly untreated
# Overall gene expression is conserved between single cell and bulk rnaseq
if(TRUE){
    fname<-"data/rpms.rda"
    if(!file.exists(fname)){
        rpms<-apply(rawcounts,2,function(x){1E6*x/sum(x)})
        save(rpms,file=fname)
    } else {load(fname)}
    grpm<-rpms[,names(hashtags)[hashtags=="ctr"]]
    grpm<-apply(grpm,1,mean)
    hrpm<-1E6*harenza[,"KELLY"]/sum(harenza[,"KELLY"])
    common<-intersect(names(grpm),names(hrpm))
    grpm<-grpm[common]
    hrpm<-hrpm[common]
    
    source("../shared/functions/qol.R")
    x<-log10(grpm+0.01)
    y<-log10(hrpm+0.01)
    n<-15
    up<-names(sort((rank(x)+rank(y))/2,dec=TRUE))[1:n]
    extra<-c("MYCN","ACTB")
    toshow<-union(up,extra)
    lowexp<-c("MYC","MYCL")
    
    png("plots/012_RPM_Harenza_vs_Giorgi.png",w=3000,h=3000,res=550)
    plot(x,y,pch=20,col="#00000011",
         xlab="Log10 RPM (Single Cell RNA-Seq)",
         ylab="Log10 RPM (Bulk RNA-Seq, Harenza et al. Dataset)",
         main="Overall Gene Expression in Untreated Kelly Cells",
         xlim=c(-3,5),ylim=c(-2,4.5)
    )
    scc<-cor.test(x,y,method="s")
    mtext(paste0("SCC=",signif(scc$estimate,3)," (p=",signif(scc$p.value,3),")"))
    lm1<-lm(y~x)
    abline(lm1$coef,lwd=2,col="#FF0000AA",lty=1)
    textplot3(x[toshow],y[toshow],words=toshow,font=2,cex=0.8)
    textplot3(x[lowexp],y[lowexp],words=lowexp,font=2,cex=0.8,col="darkgrey")
    dev.off()
    
    ### What is the table of correlations with Harenza?
    grpm<-rpms[,names(hashtags)[hashtags=="ctr"]]
    grpm<-apply(grpm,1,mean)
    output<-matrix(NA,nrow=ncol(harenza),ncol=4)
    colnames(output)<-c("PCC","PCC p-value","SCC","SCC p-value")
    rownames(output)<-colnames(harenza)
    for (i in 1:ncol(harenza)){
        cellname<-colnames(harenza)[i]
        message("Doing ",cellname)
        hrpm<-1E6*harenza[,cellname]/sum(harenza[,cellname])
        common<-intersect(names(grpm),names(hrpm))
        grpm<-grpm[common]
        hrpm<-hrpm[common]
        x<-log10(grpm+0.01)
        y<-log10(hrpm+0.01)
        
        ii<-1
        for(cortype in c("p","s")){
            cortest<-cor.test(x,y,method=cortype)
            coeff<-signif(cortest$estimate,4)
            p<-signif(cortest$p.value,3)
            output[i,ii]<-coeff
            ii<-ii+1
            output[i,ii]<-p
            ii<-ii+1
        }
    }
    library(gridExtra)
    library(grid)
    output<-output[order(-output[,1]),]
    png("plots/012_cortests_harenza_vs_ctr.png",w=2000,h=5000,res=350)
    grid.table(output)
    dev.off()
    
    ### What is the table of correlations with Harenza? IN PANOBINOSTAT
    grpm<-rpms[,names(hashtags)[hashtags=="pnb"]]
    grpm<-apply(grpm,1,mean)
    output<-matrix(NA,nrow=ncol(harenza),ncol=4)
    colnames(output)<-c("PCC","PCC p-value","SCC","SCC p-value")
    rownames(output)<-colnames(harenza)
    for (i in 1:ncol(harenza)){
        cellname<-colnames(harenza)[i]
        message("Doing ",cellname)
        hrpm<-1E6*harenza[,cellname]/sum(harenza[,cellname])
        common<-intersect(names(grpm),names(hrpm))
        grpm<-grpm[common]
        hrpm<-hrpm[common]
        x<-log10(grpm+0.01)
        y<-log10(hrpm+0.01)
        
        ii<-1
        for(cortype in c("p","s")){
            cortest<-cor.test(x,y,method=cortype)
            coeff<-signif(cortest$estimate,4)
            p<-signif(cortest$p.value,3)
            output[i,ii]<-coeff
            ii<-ii+1
            output[i,ii]<-p
            ii<-ii+1
        }
    }
    library(gridExtra)
    library(grid)
    output<-output[order(-output[,1]),]
    png("plots/012_cortests_harenza_vs_pnb.png",w=2000,h=5000,res=350)
    grid.table(output,theme=ttheme_default(base_colour="red3"))
    dev.off()
    
}



rm(list=ls())




############################### Overlay Harenza centroids over Giorgi dataset
#windowsFonts("xkcd" = windowsFont("xkcd-Regular"))
# Load Giorgi annotation
load("data/seurat.rda") # hashtags,phases,nrnas,rawcounts,expmat,regexpmat
# In regexpmat, we kept all genes expressed in >=3 cells and all cells with at least 1000 detected genes
# In regexpmat, cell cycle and nr. transcripts were regressed out
dim(regexpmat) # 18481  2183
dim(expmat) # 18557  2389
expmat<-expmat[rownames(regexpmat),colnames(regexpmat)]
dim(expmat) # 18481  2183
giorgi<-expmat


# Load RPMs for both datasets
load("data/harenza/expmat.rda")
harenza<-expmat
harenza<-harenza[,setdiff(colnames(harenza),"COG-N-440")] # Outlier
rm(expmat)

# Make them comparable (same genes)
common<-intersect(rownames(harenza),rownames(giorgi))
harenza<-harenza[common,]
giorgi<-giorgi[common,]



### Plot without normalization/regression, Bulk, Pnb and Control
if(TRUE){
    # Merge them
    expmat<-cbind(giorgi,harenza)
    # Rtsne
    topvars<-names(sort(apply(expmat,1,var),decreasing=TRUE))[1:5000]
    tsnemat<-expmat[topvars,]
    library(Rtsne)
    fname<-"results/012_rtsne.rda"
    if(!file.exists(fname)){
        set.seed(1)
        ttt<-Rtsne(t(tsnemat))
        save(ttt,file=fname)
    }else{load(fname)}
    x<-ttt$Y[,1]
    y<-ttt$Y[,2]
    # Use Seurat assignment for cell cycles
    cyclecols<-setNames(rep("magenta",ncol(tsnemat)),colnames(tsnemat))
    cyclecols[names(phases)[phases=="G1"]]<-"#fa807299"
    cyclecols[names(phases)[phases=="G2M"]]<-"#2e8b5799"
    cyclecols[names(phases)[phases=="S"]]<-"#6495ed99"
    cyclecols[names(cyclecols)%in%colnames(harenza)]<-"#00000099"
    # Sample group
    group<-setNames(rep(17,ncol(tsnemat)),colnames(tsnemat))
    group[names(hashtags)[hashtags=="pnb"]]<-18
    group[names(hashtags)[hashtags=="ctr"]]<-20
    # And plot
    png("plots/012b_tsne_with_harenza.png",w=4000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Single Cell vs. Bulk RNA-Seq",col=cyclecols,pch=group)
    legend("bottomright",pch=15,col=c("#fa807299","#2e8b5799","#6495ed99"),
           legend=c("G1","G2M","S"),bg="ivory",title="Cell Cycle Phase",cex=0.7,pt.cex=1.2)
    legend("bottomleft",pch=c(17,18,20),col="#00000099",
           legend=c("Bulk","SC pnb","SC ctrl"),bg="ivory",title="Sample Group",cex=0.7,pt.cex=1.2)
    dev.off()
}

### Plot without normalization/regression, Bulk and Control (no Pnb)
if(TRUE){
    # Merge them
    onlyctr<-intersect(colnames(giorgi),names(hashtags)[hashtags=="ctr"])
    expmat<-cbind(giorgi[,onlyctr],harenza)
    # Rtsne
    topvars<-names(sort(apply(expmat,1,var),decreasing=TRUE))[1:5000]
    tsnemat<-expmat[topvars,]
    library(Rtsne)
    fname<-"results/012_rtsne2.rda"
    if(!file.exists(fname)){
        set.seed(1)
        ttt<-Rtsne(t(tsnemat))
        save(ttt,file=fname)
    }else{load(fname)}
    x<-setNames(ttt$Y[,1],colnames(tsnemat))
    y<-setNames(ttt$Y[,2],colnames(tsnemat))
    # Use Seurat assignment for cell cycles
    cyclecols<-setNames(rep("magenta",ncol(tsnemat)),colnames(tsnemat))
    cyclecols[names(phases)[phases=="G1"]]<-"#fa807299"
    cyclecols[names(phases)[phases=="G2M"]]<-"#2e8b5799"
    cyclecols[names(phases)[phases=="S"]]<-"#6495ed99"
    cyclecols[names(cyclecols)%in%colnames(harenza)]<-"#00000099"
    # Sample group
    group<-setNames(rep(17,ncol(tsnemat)),colnames(tsnemat))
    group[names(hashtags)[hashtags=="ctr"]]<-20
    # And plot
    png("plots/012c_tsne_with_harenza_nopnb.png",w=4000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Single Cell vs. Bulk RNA-Seq",col=cyclecols,pch=group)
    legend("topleft",pch=15,col=c("#fa807299","#2e8b5799","#6495ed99"),
           legend=c("G1","G2M","S"),bg="ivory",title="Cell Cycle Phase",cex=0.7,pt.cex=1.2)
    legend("left",pch=c(17,20),col="#00000099",
           legend=c("Bulk","Single Cell"),bg="ivory",title="Sample Group",cex=0.7,pt.cex=1.2)
    dev.off()
}



###### Remove the bulk/sc differences. Only in control, of course!
onlyctr<-intersect(colnames(giorgi),names(hashtags)[hashtags=="ctr"])
expmat<-cbind(giorgi[,onlyctr],harenza)
origin<-setNames(c(rep("sc",length(onlyctr)),rep("bulk",ncol(harenza))),colnames(expmat))
# Regression, gene by gene
fname<-"results/012_regression.rda"
if(!file.exists(fname)){
    regression<-apply(expmat,1,function(vector){
        mylm<-lm(vector~origin)
        return(mylm$residuals)
    })
    save(regression,file=fname)    
}else{load(fname)}
### Plot with normalization/regression
if(TRUE){
    expmat<-t(regression)
    # Rtsne
    topvars<-names(sort(apply(expmat,1,var),decreasing=TRUE))[1:5000]
    tsnemat<-expmat[topvars,]
    library(Rtsne)
    fname<-"results/012_rtsne3.rda"
    if(!file.exists(fname)){
        set.seed(1)
        ttt<-Rtsne(t(tsnemat))
        save(ttt,file=fname)
    }else{load(fname)}
    x<-setNames(ttt$Y[,1],colnames(tsnemat))
    y<-setNames(ttt$Y[,2],colnames(tsnemat))
    # Use Seurat assignment for cell cycles
    cyclecols<-setNames(rep("magenta",ncol(tsnemat)),colnames(tsnemat))
    cyclecols[names(phases)[phases=="G1"]]<-"#fa807299"
    cyclecols[names(phases)[phases=="G2M"]]<-"#2e8b5799"
    cyclecols[names(phases)[phases=="S"]]<-"#6495ed99"
    cyclecols[names(cyclecols)%in%colnames(harenza)]<-"#000000AA"
    # And plot
    source("../shared/functions/qol.R")
    png("plots/012d_tsne_with_harenza_regressed.png",w=4000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Single Cell vs. Bulk RNA-Seq",col=cyclecols,pch=20,
         xlim=1.2*c(min(x),max(x)),ylim=1.2*c(min(y),max(y))
    )
    legend("topleft",pch=20,col=c("#fa807299","#2e8b5799","#6495ed99","#00000099"),
           legend=c("Single Cell G1","Single Cell G2M","Single Cell S","Bulk"),bg="ivory",title="Legend",cex=0.7,pt.cex=1.2)
    dev.off()
    
    
    cyclecols[names(cyclecols)%in%colnames(harenza)]<-"#00000033"
    png("plots/012d_tsne_with_harenza_regressed_withlabels.png",w=4000,h=3000,res=600)
    plot(x,y,xlab="TSNE1",ylab="TSNE2",main="Single Cell vs. Bulk RNA-Seq",col=cyclecols,pch=20,
         xlim=1.2*c(min(x),max(x)),ylim=1.2*c(min(y),max(y))
    )
    legend("topleft",pch=20,col=c("#fa807299","#2e8b5799","#6495ed99","#00000099"),
           legend=c("Single Cell G1","Single Cell G2M","Single Cell S","Bulk"),bg="ivory",title="Legend",cex=0.7,pt.cex=1.2)
    # toshow<-c("KELLY","SK-N-BE-2--C")
    # textplot3(x[toshow],y[toshow],words=c("KELLY","BE2C"))
    toshow<-colnames(harenza)
    textplot3(x[toshow],y[toshow],words=toshow,rstep=1.0,cex=0.5,font=2,padding="  ")
    dev.off()
    
    
}


rm(list=ls())


############################### What if the SUM of single cell equals bulk rnaseq?
windowsFonts("xkcd" = windowsFont("xkcd-Regular"))

load("data/harenza/rawcounts.rda")
harenza<-rawcounts
wass<-load("data/seurat.rda")

### Sum of control cells
ctr<-rawcounts[,intersect(colnames(rawcounts),names(hashtags)[hashtags=="ctr"])]
pnb<-rawcounts[,intersect(colnames(rawcounts),names(hashtags)[hashtags=="pnb"])]
ctr<-apply(ctr,1,sum)
pnb<-apply(pnb,1,sum)
common<-intersect(names(ctr),rownames(harenza))
newmat<-cbind(harenza[common,],ctr[common],pnb[common])
colnames(newmat)[(ncol(newmat)-1):ncol(newmat)]<-c("Kelly_SC_ctr","Kelly_SC_pnb")
source("../shared/functions/vst.R")
expmat<-vst(newmat)

### Rtsne
# Rtsne
topvars<-names(sort(apply(expmat,1,var),decreasing=TRUE))[1:5000]
tsnemat<-expmat[topvars,]
library(Rtsne)
set.seed(1)
ttt<-Rtsne(t(tsnemat),perplexity=10)
x<-setNames(ttt$Y[,1],colnames(tsnemat))
y<-setNames(ttt$Y[,2],colnames(tsnemat))

source("../shared/functions/qol.R")
png("plots/012e_tsne_bulk.png",w=3000,h=3000,res=600)
plot(x,y,pch=20,xlab="TSNE1",ylab="TSNE2",main="Neuroblastoma Cell Lines",
     xlim=1.4*c(min(x),max(x)),ylim=1.2*c(min(y),max(y)),
     col=c(rep("#00000099",length(x)-2),"cornflowerblue","red3")
)
set.seed(3)
textplot3(x,y,words=names(x),cex=0.8,
          font=c(rep(1,length(x)-2),2,2),
          col=c(rep("black",length(x)-2),"cornflowerblue","red3"),padding="  "
)
dev.off()







