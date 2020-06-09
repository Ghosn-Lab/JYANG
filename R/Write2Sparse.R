library(Matrix)
library(Seurat)

setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq")
PBMC3<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
setwd('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature')
DSB1_matrix<-read.csv("DSB_adt+vdj_PBMC3.csv",row.names = 1)
colnames(PBMC3$`Gene Expression`)<-sub("-1$","_3",colnames(PBMC3$`Gene Expression`))

PBMC3_5000<-CreateSeuratObject(counts = PBMC3$`Gene Expression`)
# a<-matrix(1,1,length(PBMC3_5000$orig.ident),dimnames = list("PBMC3",names(PBMC3_5000$orig.ident)))
# DSB1_matrix<-rbind(DSB1_matrix,a)
PBMC3_5000[['Protein']] = CreateAssayObject(counts = DSB1_matrix)
PBMC3_5000 <- NormalizeData(PBMC3_5000, normalization.method = "LogNormalize", scale.factor = 10000)
PBMC3_5000<-FindVariableFeatures(PBMC3_5000,selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(PBMC3_5000)
PBMC3_var<-PBMC3_5000@assays$RNA@var.features
PBMC3_5000 <- ScaleData(PBMC3_5000, features = all.genes)




PBMC3_5000<-SCTransform(PBMC3_5000,verbose = T)
PBMC3_5000<-RunPCA(PBMC3_5000,npcs = 50,verbose = F)

raw<-PBMC3_5000@assays$RNA@counts
Log<-as(PBMC3_5000@assays$RNA@scale.data,"sparseMatrix")
ADT<-PBMC3_5000@assays$Protein@counts
SCT<-as(PBMC3_5000@assays$SCT@scale.data,"sparseMatrix")
SCT_PCA<-as(t(PBMC3_5000@reductions$pca@cell.embeddings),"sparseMatrix")

ADT_raw<-rbind(raw,ADT)
ADT1_Log<-rbind(Log,ADT)
ADT_raw_pca<-rbind(ADT_raw,SCT_PCA)




setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/ADT+GEX")
###Sparse matrix
writeMM(ADT_raw,file ='PBMC3_raw/matrix.mtx')
writeMM(ADT_raw_pca,'PBMC3_SCT/matrix.mtx')
writeMM(ADT_Log,'PBMC3_log/matrix.mtx')

### feature names
feature_raw<-data.frame(a=rownames(raw),b=rownames(raw),c="Gene Expression")
feature_ADT<-data.frame(a=rownames(ADT),b=rownames(ADT),c="Antibody Capture")
feature_SCT_PCA<-data.frame(a=rownames(SCT_PCA),b=rownames(SCT_PCA),c="Dimensionality Reduction")
Feature_raw<-rbind(feature_raw,feature_ADT)
Feature_raw_pca<-rbind(Feature_raw,feature_SCT_PCA)

write.table(Feature_raw,file = "PBMC3_raw/features.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(Feature_raw,file = "PBMC3_log/features.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(Feature_raw_pca,file = "PBMC3_SCT/features.tsv",sep = "\t",col.names = F,row.names = F,quote = F)

### barcodes
UV<-c(PBMC3_var,PBMC3_var,PBMC3_var)
Var<-unique(UV)

barcodes<-data.frame(a=colnames(raw))
write.table(barcodes,file = "PBMC3_raw/barcodes.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(barcodes,file = "PBMC3_SCT/barcodes.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
write.table(barcodes,file = "PBMC3_log/barcodes.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
###########
###########
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq")
PBMC3<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
setwd('C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/new_feature')
DSB1_matrix<-read.csv("DSB_adt+vdj_PBMC3.csv",row.names = 1)
colnames(PBMC3$`Gene Expression`)<-paste(colnames(PBMC3$`Gene Expression`),"3",sep = "_")

PBMC3_5000<-CreateSeuratObject(counts = PBMC3$`Gene Expression`)
PBMC3_5000[['Protein']] = CreateAssayObject(counts = DSB1_matrix)

ADT<-cbind(ADT1,ADT2)
ADT<-cbind(ADT,ADT3)

Log<-cbind(Log1[Var,],Log2[Var,])
Log<-cbind(Log,Log3[Var,])

Big<-cbind(ADT_raw2,ADT_raw3)
Big<-cbind(Big,ADT_raw)

feature_gene<-data.frame(a=rownames(raw),b=rownames(raw),c="Gene Expression")
feature_ADT<-data.frame(a=rownames(ADT),b=rownames(ADT),c="Antibody Capture")
Feature_log<-rbind(feature_gene,feature_ADT)
write.table(Feature_log,file = "PBMC323_log/features.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
writeMM(Big,file ='PBMC323_log/matrix.mtx')
barcodes<-data.frame(a=colnames(Big))
write.table(barcodes,file = "PBMC323_log/barcodes.tsv",sep = "\t",col.names = F,row.names = F,quote = F)
