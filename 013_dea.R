
source("../shared/functions/qol.R")

load("results/001_shahbazi_resannot.rda")
shahbazi<-resannot
rm(resannot)

load("data/seurat.rda") # hashtags,phases,nrnas,rawcounts,expmat,regexpmat
# Focus on normalized matrix
dim(regexpmat) # 18481  2183
dim(expmat) # 18557  2389

# Remove multiplets (already removed if regexpmat)
expmat<-expmat[,names(hashtags)[hashtags!="multiplet"]]
dim(expmat) # 18557  2293

# Remove low coverage cells (already removed if regexpmat)
keep<-intersect(colnames(expmat),names(nrnas)[nrnas>=10000])
expmat<-expmat[,keep]
dim(expmat) # 18557  2183

###### Define contrast (Wilcoxon tests)
### Function to calculate Differential expression with wilcoxon tests
wexp<-function(matx,maty){
    columns<-c("log2fc","wstat","p","fdr")
    output<-matrix(NA,nrow=nrow(matx),ncol=length(columns))
    colnames(output)<-columns
    rownames(output)<-rownames(matx)
    pb<-txtProgressBar(0,nrow(matx),style=3)
    for(i in 1:nrow(matx)){
        x<-matx[i,]
        y<-maty[i,]
        l2fc<-log2(mean(x)/mean(y))
        wt<-wilcox.test(x,y)
        p<-wt$p.value
        if(p==0){p<-1e-301}
        stat<--log10(p)*sign(l2fc)
        output[i,]<-c(l2fc,stat,p,NA)
        setTxtProgressBar(pb,i)
    }
    output[,"fdr"]<-p.adjust(output[,"p"],method="BH")
    return(as.data.frame(output,stringsAsFactors=FALSE))
}
### pnb vs ctr
fname<-"results/013_res_pnb_vs_ctr.rda"
if(!file.exists(fname)){
    pnb<-intersect(names(hashtags)[hashtags=="pnb"],colnames(expmat))
    pnb<-expmat[,pnb]
    ctr<-intersect(names(hashtags)[hashtags=="ctr"],colnames(expmat))
    ctr<-expmat[,ctr]
    res<-wexp(pnb,ctr)
    save(res,file=fname)
}else{load(fname)}

# Classic volcano plot
top<-setNames(res$wstat,rownames(res))
up<-names(sort(top,dec=TRUE))[1:10]
dn<-names(sort(top,dec=FALSE))[1:10]
dn<-c(dn,"MYCN")
xxx<-setNames(res$log2fc,rownames(res))
yyy<-setNames(-log10(res$fdr),rownames(res))

whichup<-names(xxx)[xxx>1&yyy>=2]
whichdn<-names(xxx)[xxx<(-1)&yyy>=2]

png("plots/013_volcano.png",w=3000,h=3000,res=600)
plot(xxx,yyy,xlab="Log2FC",ylab="-log10(pFDR)",main="Panobinostat vs. Control",pch=20,col="#00000011",ylim=c(0,350),xlim=c(-7,7))
abline(h=0)
textplot3(xxx[up],yyy[up],words=up,cex=0.8,col=2,font=2)
textplot3(xxx[dn],yyy[dn],words=dn,cex=0.8,col=4,font=2)
points(xxx[whichup],yyy[whichup],col="#99000033",pch=20)
points(xxx[whichdn],yyy[whichdn],col="#00009933",pch=20)
mtext(paste0("Genes upregulated: ",length(whichup),", downregulated: ",length(whichdn)," (pFDR<=0.01, |log2FC|>1)"),cex=0.7,font=2)
dev.off()

### First test: compare with CellRanger
cellranger<-read.csv("cellranger/Pnb-VS-CTRL_Sig_Genes.csv",as.is=TRUE)
cr<-setNames(cellranger$Pnb.Log2.Fold.Change,cellranger$FeatureName)
my<-setNames(res$log2fc,rownames(res))
common<-intersect(names(cr),names(my))
cr<-cr[common]
my<-my[common]
png("plots/013_cellranger_vs_seurat.png",w=3000,h=3000,res=500)
plot(cr,my,xlab="CellRanger",ylab="Seurat Wilcoxon",main="Log2FC Panobinostat dataset")
dev.off()


### Second test: compare with Shahbazi
giorgi<-setNames(-log10(res$p)*sign(res$wstat),rownames(res))
shah<-setNames(-log10(shahbazi$P.Value)*sign(shahbazi$t),shahbazi$symbols)
common<-intersect(names(giorgi),names(shah))
giorgi<-giorgi[common]
shah<-shah[common]

png("plots/013_giorgi_vs_shahbazi.png",w=3000,h=3000,res=600)
par(mar=c(5,5,3,1))
plot(giorgi,shah,
     xlab=expression('Current Dataset: -log'[10]*'(p) * sign of log'[FC]),
     ylab=expression('Shahbazi Dataset: -log'[10]*'(p) * sign of log'[FC]),
     main="DE Statistics, Panobinostat vs. Control",pch=20,col="#00000088",
     xlim=c(1.1*c(min(giorgi),max(giorgi))))
scc<-cor.test(giorgi,shah,method="s")
mtext(paste0("SCC=",signif(scc$estimate,3)," (p=",signif(scc$p.value,3),")"))
lm1<-lm(shah~giorgi)
abline(lm1$coef,lwd=2,col="#FF0000AA",lty=1)
abline(h=c(2,-2),v=c(2,-2),lty=2,col="#00000033")
dev.off()


## Show some genes
source("../shared/functions/qol.R")
# library(extrafont)
#font_import()
# fonts()
# loadfonts(device="win")
# https://www.reddit.com/r/xkcd/comments/2i1i37/download_and_use_randalls_handwriting_as_a_font/
# saved in $DROPBOX/projects/shared/xkcd-Regular.otf
# windowsFonts("xkcd" = windowsFont("xkcd-Regular"))

n<-15
up<-names(sort((rank(giorgi)+rank(shah))/2,dec=TRUE))[1:n]
dn<-names(sort((rank(giorgi)+rank(shah))/2,dec=FALSE))[1:n]
up2<-names(sort(rank(giorgi),dec=TRUE))[1:n]
dn2<-names(sort(rank(giorgi),dec=FALSE))[1:n]
up<-union(up,up2)
dn<-union(dn,dn2)
add<-c("HAND2-AS1")
toshow<-union(add,union(up,dn))

png("plots/013_giorgi_vs_shahbazi_withgenes.png",w=3000,h=3000,res=600)
par(mar=c(5,5,3,1))
plot(giorgi,shah,
     xlab=expression('Current Dataset: -log'[10]*'(p) * sign of log'[FC]),
     ylab=expression('Shahbazi Dataset: -log'[10]*'(p) * sign of log'[FC]),
     main="DE Statistics, Panobinostat vs. Control",pch=20,col="#00000044",
     xlim=c(1.1*c(min(giorgi),max(giorgi))))
scc<-cor.test(giorgi,shah,method="s")
mtext(paste0("SCC=",signif(scc$estimate,3)," (p=",signif(scc$p.value,3),")"))
lm1<-lm(shah~giorgi)
abline(lm1$coef,lwd=2,col="#FF0000AA",lty=1)
abline(h=c(2,-2),v=c(2,-2),lty=2,col="#00000033")
set.seed(1)
textplot3(giorgi[toshow],shah[toshow],words=toshow,cex=0.5,font=2,line.col="#00000033",rstep=0.1,padding=" ")
dev.off()



### Annotate results and save
source("../shared/functions/geneids.R")
options(java.parameters = "-Xmx8G")
library(xlsx)
entrezs<-sym2eg(rownames(res))
annotations<-annotategene(entrezs)
resannot<-cbind(rownames(res),entrezs,as.matrix(annotations),as.matrix(res))
colnames(resannot)[1:3]<-c("SYMBOL","ENTREZID","NAME")
resannot<-as.data.frame(resannot,stringsAsFactors=FALSE)
resannot$log2fc<-as.numeric(resannot$log2fc)
resannot$wstat<-as.numeric(resannot$wstat)
resannot$p<-as.numeric(resannot$p)
resannot$fdr<-as.numeric(resannot$fdr)
resannot<-resannot[order(resannot$p),]
write.xlsx2(resannot,file="results/013_resannot.xlsx",row.names=FALSE,sheetName="PNB_vs_CTRL")
save(resannot,file="results/013_resannot.rda")


