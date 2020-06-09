library(Seurat)
library(qusage)
library(ggplot2)
library(Matrix)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/filtered/")



########## Log normalization
filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM1",gene.column = 2, unique.features = TRUE)
colnames(filtered$`Gene Expression`)<-sub("-1$","_1",colnames(filtered$`Gene Expression`))
C7<-CreateSeuratObject(counts = filtered$`Gene Expression`)
rm(filtered)

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM2",gene.column = 2, unique.features = TRUE)
colnames(filtered$`Gene Expression`)<-sub("-1$","_2",colnames(filtered$`Gene Expression`))
C5<-CreateSeuratObject(counts = filtered$`Gene Expression`)
rm(filtered)
 
filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM3",gene.column = 2, unique.features = TRUE)
colnames(filtered$`Gene Expression`)<-sub("-1$","_3",colnames(filtered$`Gene Expression`))
D5<-CreateSeuratObject(counts = filtered$`Gene Expression`)
rm(filtered)

C7<-AddMetaData(C7,'BM1',col.name = 'batch')
C5<-AddMetaData(C5,'BM2',col.name = 'batch')
D5<-AddMetaData(D5,'BM3',col.name = 'batch')

C7<-SubsetData(C7,subset.name = "nFeature_RNA",low.threshold = 40)
C7<-SubsetData(C7,subset.name = "nCount_RNA",low.threshold = 100)
C5<-SubsetData(C5,subset.name = "nFeature_RNA",low.threshold = 40)
C5<-SubsetData(C5,subset.name = "nCount_RNA",low.threshold = 100)
D5<-SubsetData(D5,subset.name = "nFeature_RNA",low.threshold = 40)
D5<-SubsetData(D5,subset.name = "nCount_RNA",low.threshold = 100)


C7 <- NormalizeData(C7, normalization.method = "LogNormalize", scale.factor = 10000)
C7 <- FindVariableFeatures(C7, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(C7)
C7 <- ScaleData(C7, features = all.genes)
cv<-C7@assays$RNA@var.features

C5 <- NormalizeData(C5, normalization.method = "LogNormalize", scale.factor = 10000)
C5 <- FindVariableFeatures(C5, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(C5)
C5 <- ScaleData(C5, features = all.genes)
c5v<-C5@assays$RNA@var.features

D5 <- NormalizeData(D5, normalization.method = "LogNormalize", scale.factor = 10000)
D5 <- FindVariableFeatures(D5, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(D5)
D5 <- ScaleData(D5, features = all.genes)
dv<-D5@assays$RNA@var.features

Uv<-unique(c(cv,c5v,dv))

C7<-SCTransform(C7,verbose = F)
C5<-SCTransform(C5,verbose = F)
D5<-SCTransform(D5,verbose = F)

Uv_SCT<-unique(c(C7@assays$SCT@var.features,C5@assays$SCT@var.features,D5@assays$SCT@var.features))

###var scale Merge 
###
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
###var genes prior merge
protease@assays$RNA@var.features<-Uv
###var genes after merge
# protease <- FindVariableFeatures(protease, selection.method = "vst", nfeatures = 3000)

Mscale<-cbind(C7@assays$RNA@scale.data[protease@assays$RNA@var.features,],C5@assays$RNA@scale.data[protease@assays$RNA@var.features,])
Mscale<-cbind(Mscale,D5@assays$RNA@scale.data[protease@assays$RNA@var.features,])

protease@assays$RNA@scale.data<-Mscale


protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)

DimPlot(protease,group.by = "batch")

### Merge var and scale
###
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease <- FindVariableFeatures(protease, selection.method = "vst", nfeatures = 3000)
var.genes <- protease@assays$RNA@var.features
protease <- ScaleData(protease, features = var.genes)
protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)

DimPlot(protease,group.by = "batch")

### var Merge scale
###
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease@assays$RNA@var.features<-Uv
var.genes <- protease@assays$RNA@var.features
protease <- ScaleData(protease, features = var.genes)
protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)

DimPlot(protease,group.by = "batch")

### Merge SCT

protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease<-SCTransform(protease,verbose = F)
protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)
DimPlot(protease,group.by = "batch")

### SCT Merge(default merge)
# 
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease@assays$SCT@var.features<-rownames(protease@assays$SCT@scale.data)
protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)
DimPlot(protease,group.by = "batch")

### SCT Merge(full merge)
# 
C7<-SCTransform(C7,verbose = F,return.only.var.genes = F)
C5<-SCTransform(C5,verbose = F,return.only.var.genes = F)
D5<-SCTransform(D5,verbose = F,return.only.var.genes = F)

Uv_SCT<-unique(c(C7@assays$SCT@var.features,C5@assays$SCT@var.features,D5@assays$SCT@var.features))
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease@assays$SCT@var.features<-Uv_SCT
protease@assays$SCT@var.features<-intersect(Uv_SCT,rownames(protease@assays$SCT@scale.data))

Mscale<-merge(C7@assays$SCT@scale.data,C5@assays$SCT@scale.data,all=T,by="row.names")
rownames(Mscale)<-Mscale$Row.names
Mscale$Row.names<-NULL
Mscale<-merge(Mscale,D5@assays$SCT@scale.data,all = T,by="row.names")
rownames(Mscale)<-Mscale$Row.names
Mscale$Row.names<-NULL
Mscale<-as.matrix(Mscale[Uv_SCT,])
Mscale[is.na(Mscale)]<-0
protease@assays$SCT@scale.data<-Mscale

protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)
DimPlot(protease,group.by = "batch")


######################