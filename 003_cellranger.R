
### Understanding cellranger output:
### https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview

# Installing Cellranger 3.1.0 (July 24, 2019)
cd /mnt/d/Dropbox/programs
wget -O cellranger-3.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz?Expires=1566964056&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU2Njk2NDA1Nn19fV19&Signature=j7yh2-GDGi2sdek92D9H722K9pB3RK~bwnT6sZr9XFm7xstoBlCAzCnDa6xL0K7AD2d3iLCdFNvL5-gRTPxSQLE7d6Wz2l1qYuY9kgGk9PmLOAkblPVbx2hL3tQrOFWBSGJBh~8CGOm3YYUbNfycF9euD4xXivH2jwUM7c57Swt5X03qyLk~lx4bD~HJRmnHRRedVOlQSMFj6Co4y8GWubO8EKt-4D1AfbWJelXxNkW8uH1FTFDkVChFqPJKRRURx7RI~PmhJd2ykdmIaSe5xEdttbhNtn3~yrcWCLLDGnLaLzt5P79mcB1MQR4Ilkx2oI0uRZfdttr9CjPjPSleJg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar xvzf cellranger-3.1.0.tar.gz

# Installing references for Cellranger 3.1.0 (July 24, 2019)
cd /mnt/d/Dropbox/programs/
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-and-mm10-3.1.0.tar.gz

# Each matrix is stored in the Market Exchange Format (MEX) for sparse matrices.
# It also contains gzipped TSV files with feature and barcode sequences
# corresponding to row and column indices respectively. For example,
# the matrices output may look like:
# filtered_feature_bc_matrix
# ├── barcodes.tsv.gz
# ├── features.tsv.gz
# └── matrix.mtx.gz

### Convert Cell ranger files to gene counts
cd /mnt/g/delivery_20190826/outs
alias cellranger="/mnt/d/Dropbox/programs/cellranger-3.1.0/cellranger"

# # convert from MEX
# cellranger mat2csv filtered_feature_bc_matrix panobinostat.csv
# or, convert from H5
cellranger mat2csv filtered_feature_bc_matrix.h5 panobinostat.csv
# this matrix has 33540 x 2678 (89820120 total) elements, 87.283473% of which are zero.
gzip panobinostat.csv
cp panobinostat.csv.gz /mnt/d/Dropbox/projects/panobinostat/data


### Loading Matrices into R
# The R package Matrix supports loading MEX format data, and can be easily
# used to load the sparse feature-barcode matrix, as shown in the example
# code below.
library(Matrix)
matrix_dir<-"G:/delivery_20190826/outs/filtered_feature_bc_matrix/"
barcode.path<-paste0(matrix_dir,"barcodes.tsv.gz")
features.path<-paste0(matrix_dir,"features.tsv.gz")
matrix.path<-paste0(matrix_dir,"matrix.mtx.gz")
rawcounts_sparsed<-readMM(file=matrix.path)
feature.names<-read.delim(features.path,header=FALSE,stringsAsFactors=FALSE)
barcode.names<-read.delim(barcode.path,header=FALSE,stringsAsFactors=FALSE)
colnames(rawcounts_sparsed)<-barcode.names$V1
rownames(rawcounts_sparsed)<-feature.names$V1
save(rawcounts_sparsed,file="data/rawcounts_sparsed.rda")
dim(rawcounts_sparsed) # 33540  2678

# As an alternative, we can load the non-sparse CSV from the cellranger mat2csv conversion
rawcounts<-as.matrix(read.csv("data/panobinostat.csv.gz",as.is=TRUE,row.names=1))
dim(rawcounts) # 33540  2678
save(rawcounts,file="data/rawcounts_ensg.rda")

# Convert to gene symbols
source("../shared/functions/geneids.R")
source("../shared/functions/expmerger.R")
ensg<-rownames(rawcounts)
entrezs<-ens2eg(ensg)
symbols<-eg2sym(entrezs)
symbols[is.na(symbols)]<-ensg[is.na(symbols)]
names(symbols)<-ensg
length(unique(ensg)) # 33540
length(unique(symbols)) # 33460

# Sum reads mapping to the same symbol
newcounts<-squish(rawcounts,symbols,method="sum",verbose=TRUE)
rawcounts<-newcounts
save(rawcounts,file="data/rawcounts.rda")













