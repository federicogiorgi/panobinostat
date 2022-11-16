### Load a list of tfs (ordered as in 030_tabellone_forGiorgio) to indicate RPM and %cells
tfs<-read.delim("results/030_tflist.txt",as.is=TRUE,header=FALSE)[,1]

# Load single cell data
load("data/seurat.rda")
rpms<-apply(rawcounts,2,function(x){1E6*x/sum(x)})
pnb<-names(hashtags[hashtags=="pnb"])
ctr<-names(hashtags[hashtags=="ctr"])
output<-cbind(
    apply(rpms[tfs,ctr],1,mean),
    apply(rpms[tfs,pnb],1,mean)
)

# % Cells
perc<-100*apply(rawcounts[tfs,],1,function(x){sum(x>0)/length(x)})
output<-cbind(output,perc)

rm(rawcounts)

# Load # Harenza
load("data/harenza/rawcounts.rda")
harenza<-rawcounts[tfs,"KELLY"]
harenza<-(1E6*harenza)/sum(harenza)
output<-cbind(output,harenza)

# Print
colnames(output)<-c(
    "RPM_scCtrl","RPM_scPnb","% scCells","RPM_bulk_KELLY"
)
write.table(output,file="results/031_extraColumns.txt",quote=FALSE,sep="\t")
