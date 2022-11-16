### The analysis should predict transcription factors that are
# Downregulated by panobinostat (expression in Shahbazi)
# Downregulated by panobinostat (activity, TARGET inference)
# Downregulated by panobinostat (activity, NRC inference)
# Negatively associated to survival (four NBL datasets)
# Negatively associated to survival (pan-cancer)
# Upregulated in cancer (orphan TFs)
# Targetable in Drugbank

source("../shared/functions/qol.R")
options(java.parameters = "-Xmx10000m");library(xlsx)


### Gene expression analysis
# Tao
load("results/001_shahbazi_resannot.rda")
tao<-setNames(p2z(resannot$P.Value)*sign(resannot$logFC),resannot$symbols)
rm(resannot)

# Us (without regression)
load("results/013_res_pnb_vs_ctr.rda")
exp<-setNames(p2z(res$p)*sign(res$log2fc),rownames(res))
rm(res)

### Master regulator analysis
# # On original data
# load("results/016_mra1_original.rda")

# On cell cycle regressed data
load("results/016_mra2_regressed.rda")
tmra<-tmra$es$nes
nmra<-nmra$es$nes


### Survival analysis
load("../masterset/pansurvival/pansurvival.rda")
surv<-t(apply(pansurvival,1,function(x){
    x<-z2p(abs(x))*sign(x)
    return(x)
}))
rm(pansurvival)

# Integrate survival
survnbl<-apply(surv[,grep("NBL",colnames(surv))],1,function(x){stouffer(-x)})
survtcga<-apply(surv[,grep("tcga",colnames(surv))],1,function(x){stouffer(-x)})


### Let's limit this to TFs
tfs<-intersect(names(tmra),names(nmra))
tfs<-intersect(tfs,rownames(surv))
length(tfs) # 1331
tmra<-tmra[tfs]
nmra<-nmra[tfs]
tao<-tao[tfs]
exp<-exp[tfs]
survnbl<-survnbl[tfs]
survtcga<-survtcga[tfs]

# Exp vs. Activity plot
toshow<-c("HCFC1","BAZ1A","MAZ","ZNF146")
png("plots/030_exp_vs_mra.png",w=3000,h=3000,res=600)
plot(exp,tmra,pch=20,col="#00000033",xlab="Differential Expression (NES)",
     ylab="Differential Activity (NES)",main="PNB vs. Control Single Cell Dataset")
textplot3(exp[toshow],tmra[toshow],words=toshow,font=2,col=4)
dev.off()


### Tabellone
tabellone<-cbind(tao,exp,tmra,nmra,survnbl,survtcga)
allranks<-apply(tabellone,2,rank)
meanrank<-apply(allranks,1,mean)
tabellone<-cbind(tabellone,meanrank)
rownames(tabellone)<-tfs

tabellone<-tabellone[!is.na(tabellone[,"exp"]),]
tabellone<-tabellone[order(as.numeric(tabellone[,"tmra"])),]

# Drugbank
load("data/drugbank.rda")
drugs<-setNames(rep(NA,nrow(tabellone)),rownames(tabellone))
for(t in rownames(tabellone)){
    if(t%in%names(drugbank)){
        tdrugs<-paste0(drugbank[[t]],collapse=" ; ")
        drugs[t]<-tdrugs
    } else {
        drugs[t]<-NA
    }
}
tabellone<-cbind(tabellone, drugs)


# Final format
tabellone<-as.data.frame(tabellone)
for(i in 1:7){
    tabellone[,i]<-as.numeric(tabellone[,i])
}
colnames(tabellone)<-c("Tao_Expression_NES","Milazzo_Expression_NES","Milazzo_TARGET_TFactivity_NES",
                       "Milazzo_NRC_TFactivity_NES","Survival_NBL_integrated_NES",
                       "Survival_TCGA_integrated_NES","MeanRank_Evidences","DrugBank"
                       )
tabellone<-tabellone[order(tabellone$MeanRank_Evidences),]

save(tabellone,file="data/tabellone.rda")
write.xlsx2(tabellone,file="results/030_tabellone.xlsx")







