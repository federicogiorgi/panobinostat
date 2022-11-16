### We operate on the dataset scaled by cell cycle score and mRNA counts (regexpmat)

# VIPER activity (two networks: TARGET and NRC)
load("data/vipermats.rda")

# Candidate TFs
cands<-c("BAZ1A","HCFC1","MAZ","MYCN","ZNF146")

# Classic Rtsne to visualize progress
fname<-"data/regexpmat_plus_tsne.rda"
if(!file.exists(fname)){
    load("data/seurat.rda")
    rm(expmat)
    vargenes<-apply(regexpmat,1,var)
    topvargenes<-names(sort(vargenes,dec=TRUE))[1:5000]
    tsnemat<-regexpmat[topvargenes,]
    library(Rtsne)
    set.seed(2)
    ttt<-Rtsne(t(tsnemat))
    hashtags<-hashtags[colnames(regexpmat)]
    nrnas<-nrnas[colnames(regexpmat)]
    phases<-phases[colnames(regexpmat)]
    rawcounts<-rawcounts[,colnames(regexpmat)]
    rpms<-apply(rawcounts,2,function(x){1E6*x/sum(x)})
    rownames(ttt$Y)<-colnames(regexpmat)
    save(ttt,regexpmat,rawcounts,rpms,hashtags,nrnas,phases,file=fname)
} else {load(fname)}


# Plot Rtsne
x<-ttt$Y[,2]
y<-ttt$Y[,1]
pchs<-hashtags
pchs[hashtags=="ctr"]<-21
pchs[hashtags=="pnb"]<-22
pchs<-setNames(as.numeric(pchs),names(hashtags))
plot(x,y,pch=pchs,col="gray40")

# Centroid coordinates in TSNE space
centroid_ctr<-apply(ttt$Y[hashtags=="ctr",],2,mean)
points(centroid_ctr[1],centroid_ctr[2],pch=19,cex=3)
centroid_pnb<-apply(ttt$Y[hashtags=="pnb",],2,mean)
points(centroid_pnb[1],centroid_pnb[2],pch=15,cex=3)

# Furthest point from centroids (start and end)
distcalc<-rbind(ttt$Y,centroid_ctr,centroid_pnb)
distmat<-as.matrix(dist(distcalc,method="euclidean"))
start<-names(which.max(distmat["centroid_pnb",]))
end<-names(which.max(distmat[start,]))

# Plot start and end
plot(x,y,pch=pchs,col="gray40")
points(x[start],y[start],pch="S",cex=2)
points(x[end],y[end],pch="E",cex=2)

# Loop to include further cells, increasingly distant from start
distmat<-as.matrix(dist(ttt$Y,method="euclidean"))
faraway<-sort(distmat[start,])
tot<-length(faraway) # 2183

# Which MRs to show
load("data/tabellone.rda")
tabellone<-tabellone[order(as.numeric(tabellone[,"MeanRank_Evidences"])),]
toshow<-cands
#scalemat<-t(scale(t(regexpmat)))
scalemat<-t(scale(t(tvipermat)))

# Frames and steps
nrsteps<-100
nrcells<-20
#steps<-cut(faraway,nrsteps)
steps<-split(faraway,ceiling(seq_along(faraway)/nrcells))

### Make a movie
png("plots/032/input%03d.png",width=2000,height=1500,res=250)
minmax<-c(-2,2)
pb<-txtProgressBar(0,length(unique(steps)),style=3)
#for(i in 1:length(unique(steps))){
for(i in 1:length(steps)){
    setTxtProgressBar(pb,i)
    # Wavefront
    par(mar=c(5,5,3,1))
    lmatrix<-rbind(c(1,1),c(2,2))
    layout(lmatrix)
    plot(x,y,pch=pchs,col="gray40",xlab="TSNE1",ylab="TSNE2",main="Panobinostat effect on Kelly cells")
    text(x[start],y[start],labels="Ctrl",cex=2)
    text(x[end],y[end],labels="Pnb",cex=2)
    cells<-names(steps[[i]])
    
    # Add neighbors
    if(i>1&i<length(steps)){
        aftercells<-names(steps[[i+1]])
        beforecells<-names(steps[[i-1]])
        cells<-c(aftercells,cells,beforecells)
    }
    
    herepchs<-pchs[cells]
    herepchs[herepchs==21]<-19
    herepchs[herepchs==22]<-15
    points(x[cells],y[cells],pch=herepchs)
    legend("topright",legend=c("ctr","pnb"),pch=c(21,22))
    
    # Barplot
    toplot<-scalemat[toshow,cells]
    toplot<-apply(toplot,1,mean)
    col<-ifelse(toplot<0,"cornflowerblue","salmon")
    par(las=1,mar=c(4,6,0,2))
    barplot(toplot,horiz=TRUE,col=col,xlim=minmax,xlab="TF Score",cex.names=1.5)
    par(las=0)
}
dev.off()
png_files<-dir("plots/032",full.names=TRUE)

library(av)
av_encode_video(png_files,"plots/032_runningPlot.mp4", framerate = 24)

