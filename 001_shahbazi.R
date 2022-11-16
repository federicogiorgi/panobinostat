### Dataset by Shahbazi (Panobinostat on Neuroblstoma) GSE68690
# GSM1679076	vehicle control - replicate 1
# GSM1679078	10 nM panobinostat - replicate 1
# GSM1679080	vehicle control - replicate 2
# GSM1679082	10 nM panobinostat - replicate 2
# GSM1679084	vehicle control - replicate 3
# GSM1679086	10 nM panobinostat - replicate 3
# Affymetrix HuGene-2_0-st
# SK-N-BE(2) cells

# Create a CustomCDF annotation from CustomCDF v23.0
if(FALSE){
  library(pdInfoBuilder)
  library(foreach)
  z<-cdf2table("data/GSE68690/hugene20st_Hs_ENTREZG.cdf")
  seed<-new("GenericPDInfoPkgSeed",table=z,author="me",
            email="me@mine.org",species="Homo sapiens",pkgName="hugene20"
  )
  makePdInfoPackage(seed)
}
install.packages("hugene20/",repos=NULL,type="source")

# Load micraorrays
library(oligo)
cels<-list.celfiles("data/GSE68690/",full.names=TRUE,listGzipped=TRUE)
library(hugene20)
raw<-read.celfiles(cels,pkgname="hugene20")
# The raw container is a GeneFeatureSet

# Preprocessing
expset<-rma(raw)
colnames(expset)<-c("Ctrl1","Pano1","Ctrl2","Pano2","Ctrl3","Pano3")
expmat<-exprs(expset)

# Sample names simplified
png("plots/001_shahbazi_boxplot.png",w=1500,h=1500,p=45)
par(las=2,mar=c(5,4,3,1))
boxplot(expset,main="Shahbazi Panobinostat Dataset",ylab="RMA-normalized expression")
dev.off()

### PCA
ppp<-prcomp(t(expmat))
PC1<-ppp$x[,1]
PC2<-ppp$x[,2]
samples<-colnames(expmat)
library(ggplot2)
png("plots/001_shahbazi_pca.png",w=1200,h=1000,res=300)
q<-qplot(PC1,PC2)+geom_point(aes(colour=samples),size=5)+theme_bw()+
  ggtitle("Shahbazi Panobinostat Dataset")
print(q)
#library(directlabels)
#direct.label(q)
dev.off()


### Convert probesets to entrez ids
rownames(expmat)<-gsub("_at","",rownames(expmat))
expmat[1:5,]
save(expmat,file="results/001_shahbazi_expmat.rda")
dim(expmat) # 29635  6


### And finally, limma
library(limma)
design<-cbind(WT=1,panoVSctrl=grepl("Pano",colnames(expmat)))
rownames(design)<-colnames(expmat)
fit<-lmFit(expmat,design)
fit<-eBayes(fit)
res<-topTable(fit,coef="panoVSctrl",n=Inf)


## Annotate entrez as gene symbols
source("../shared/functions/geneids.R")
symbols<-eg2sym(rownames(res))
resannot<-cbind(symbols,res)
save(resannot,file="results/001_shahbazi_resannot.rda")
write.csv(resannot,file="results/001_shahbazi_resannot.csv")


