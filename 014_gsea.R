options(java.parameters = "-Xmx10000m");library(xlsx)


source("../shared/functions/qol.R")
## Fast GSEA library
library(fgsea)

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


###### Gene set enrichment analysis
if(FALSE){
    ### Load PNB vs CTR Giorgi signature
    load("results/013_res_pnb_vs_ctr.rda")
    load("results/001_shahbazi_resannot.rda")
    shahbazi<-resannot
    rm(resannot)
    
    ## Define signature (Giorgi)
    signature<-setNames(-log10(res$p)*sign(res$wstat),rownames(res))
    fname<-"results/014_gsea_giorgi.rda"
    if(!file.exists(fname)){
        gsea<-fgsea(pathways=mlist,stats=signature,nperm=1E6,minSize=4,maxSize=Inf,nproc=7)
        save(gsea,signature,file=fname)
    }else{load(fname)}
    shah<-gsea
    
    ## Define signature (Shahbazi)
    signature<-setNames(-log10(shahbazi$P.Value)*sign(shahbazi$t),shahbazi$symbols)
    fname<-"results/014_gsea_shah.rda"
    if(!file.exists(fname)){
        gsea<-fgsea(pathways=mlist,stats=signature,nperm=1E6,minSize=4,maxSize=Inf,nproc=7)
        save(gsea,file=fname)
    }else{load(fname)}
}


### Top pathways
load("results/014_gsea_giorgi.rda") # gsea, signature
gsea<-gsea[order(gsea$pval),]
gsea[1:20,1:7]
smallgsea<-gsea[gsea$size<=100,]
smallgsea[1:20,1:7]

write.xlsx2(gsea,file="results/014_gsea_giorgi.xlsx",sheetName="PNB_vs_CTRL",row.names=FALSE)


## Plot the top ones as a GSEA plot

# Top gseas (all sizes)
topup<-setNames(gsea$pval[gsea$NES>0],gsea$pathway[gsea$NES>0])
topup<-names(sort(topup)[1:10])
topdn<-setNames(gsea$pval[gsea$NES>0],gsea$pathway[gsea$NES<0])
topdn<-names(sort(topdn)[1:10])
top<-c(topup,topdn)

source("../shared/functions/plotEnrichment.R")
if(FALSE){ # I turned off plotting
    for(path in c(top)){
        nicepath<-gsub("_"," ",path)
        nicepath<-tools::toTitleCase(tolower(nicepath))
        
        png(paste0("plots/014_gsea/gsea_",path,".png"),w=3000,h=2000,res=600)
        gp<-plotEnrichment2(mlist[[path]],signature,ticks=0.4) +
            labs(title=nicepath,subtitle=paste0("p=",signif(gsea$padj[gsea$pathway==path]))) +
            xlab("Panobinostat Signature")
        print(gp)
        dev.off()
    }
}

# Top gseas (smaller pathways)
stopup<-setNames(smallgsea$pval[smallgsea$NES>0],smallgsea$pathway[smallgsea$NES>0])
stopup<-names(sort(stopup)[1:10])
stopdn<-setNames(smallgsea$pval[smallgsea$NES>0],smallgsea$pathway[smallgsea$NES<0])
stopdn<-names(sort(stopdn)[1:10])
stop<-c(stopup,stopdn)

source("../shared/functions/plotEnrichment.R")
if(FALSE){
    dir.create("plots/014_gsea")
    for(path in stop){
        nicepath<-gsub("_"," ",path)
        nicepath<-tools::toTitleCase(tolower(nicepath))
        
        png(paste0("plots/014_gsea/smallgsea_",path,".png"),w=3000,h=2000,res=600)
        gp<-plotEnrichment2(mlist[[path]],signature,ticks=0.4,colticks = "#00000099") +
            labs(title=nicepath,subtitle=paste0("p=",signif(smallgsea$padj[smallgsea$pathway==path]))) +
            xlab("Panobinostat Signature")
        print(gp)
        dev.off()
    }
}

# Save pathway names for plotting in 015 heatmaps
top<-union(top,stop)
save(top,file="results/014_topPathways.rda")


