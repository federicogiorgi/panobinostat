# https://www.drugbank.ca/releases/latest#protein-identifiers
# All targets, downloaded on Sep-22 2019 (DrugBank Release Version 5.1.4)
targets<-read.delim("data/drugbank/targets/all.csv",as.is=TRUE,sep=",",header=TRUE)
drugs<-read.delim("data/drugbank/drugs/drugbank vocabulary.csv",as.is=TRUE,sep=",",header=TRUE)
drugs<-setNames(drugs$Common.name,drugs$DrugBank.ID)



### Create a relationship object
dim(targets) # 4991 13
targets<-targets[targets$Species=="Humans",]
dim(targets) # 3052 13
targets<-targets[,c("Name","Gene.Name","GenAtlas.ID","Drug.IDs")]
dim(targets) # 3052 4

# Some genes are listed more than once (e.g. HTR2A)
drugbank<-list()
for(t in unique(targets$Gene.Name)){
    rows<-targets[targets$Gene.Name==t,]
    tdrugs<-unique(unlist(strsplit(rows$Drug.IDs,"; ")))
    tdrugs<-drugs[tdrugs]
    drugbank[[t]]<-tdrugs
}
save(drugbank,file="data/drugbank.rda")



