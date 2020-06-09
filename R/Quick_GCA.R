library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)

### filtered matrix
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/filtered/")
filtered<-Read10X("GCA1_prot_01788_5GEX_C7",gene.column = 2, unique.features = TRUE)
C7<-CreateSeuratObject(counts = filtered)
rm(filtered)
C7[["percent.mt"]] <- PercentageFeatureSet(C7, pattern = "^MT-")



filtered<-Read10X("GCA2_prot_21394_5GEX_C5",gene.column = 2, unique.features = TRUE)
C5<-CreateSeuratObject(counts = filtered)
rm(filtered)
C5[["percent.mt"]] <- PercentageFeatureSet(C5, pattern = "^MT-")

filtered<-Read10X("GCA3_prot_21460_5GEX_D5",gene.column = 2, unique.features = TRUE)
D5<-CreateSeuratObject(counts = filtered)
rm(filtered)
D5[["percent.mt"]] <- PercentageFeatureSet(D5, pattern = "^MT-")

filtered<-Read10X("GCA4_CD45neg_21461_5GEX_G5",gene.column = 2, unique.features = TRUE)
G5<-CreateSeuratObject(counts = filtered)
rm(filtered)
G5[["percent.mt"]] <- PercentageFeatureSet(G5, pattern = "^MT-")

filtered<-Read10X("GCA4_CD45pos_21461_5GEX_F5",gene.column = 2, unique.features = TRUE)
F5<-CreateSeuratObject(counts = filtered)
rm(filtered)
F5[["percent.mt"]] <- PercentageFeatureSet(F5, pattern = "^MT-")

filtered<-Read10X("GCA4_total_21461_5GEX_E5",gene.column = 2, unique.features = TRUE)
E5<-CreateSeuratObject(counts = filtered)
rm(filtered)
E5[["percent.mt"]] <- PercentageFeatureSet(E5, pattern = "^MT-")

filtered<-Read10X("GCA5_CD45neg_21478_5GEX_E7",gene.column = 2, unique.features = TRUE)
E7<-CreateSeuratObject(counts = filtered)
rm(filtered)
E7[["percent.mt"]] <- PercentageFeatureSet(E7, pattern = "^MT-")

filtered<-Read10X("GCA5_CD45pos_21478_5GEX_A6",gene.column = 2, unique.features = TRUE)
A6<-CreateSeuratObject(counts = filtered)
rm(filtered)
A6[["percent.mt"]] <- PercentageFeatureSet(A6, pattern = "^MT-")

filtered<-Read10X("GCA5_prot_21478_5GEX_H5",gene.column = 2, unique.features = TRUE)
H5<-CreateSeuratObject(counts = filtered)
rm(filtered)
H5[["percent.mt"]] <- PercentageFeatureSet(H5, pattern = "^MT-")

filtered<-Read10X("GCA5_total_21478_5GEX_D7",gene.column = 2, unique.features = TRUE)
D7<-CreateSeuratObject(counts = filtered)
rm(filtered)
D7[["percent.mt"]] <- PercentageFeatureSet(D7, pattern = "^MT-")

filtered<-Read10X("GCA6_CD45neg_20006_5GEX_D6",gene.column = 2, unique.features = TRUE)
D6<-CreateSeuratObject(counts = filtered)
rm(filtered)
D6[["percent.mt"]] <- PercentageFeatureSet(D6, pattern = "^MT-")

filtered<-Read10X("GCA6_CD45pos_20006_5GEX_C6",gene.column = 2, unique.features = TRUE)
C6<-CreateSeuratObject(counts = filtered)
rm(filtered)
C6[["percent.mt"]] <- PercentageFeatureSet(C6, pattern = "^MT-")

filtered<-Read10X("GCA6_total_20006_5GEX_B6",gene.column = 2, unique.features = TRUE)
B6<-CreateSeuratObject(counts = filtered)
rm(filtered)
B6[["percent.mt"]] <- PercentageFeatureSet(B6, pattern = "^MT-")

Idents(A6)<-"A6"
Idents(B6)<-"B6"
Idents(C5)<-"C5"
Idents(C6)<-"C6"
Idents(C7)<-"C7"
Idents(D5)<-"D5"
Idents(D6)<-"D6"
Idents(D7)<-"D7"
Idents(E5)<-"E5"
Idents(E7)<-"E7"
Idents(F5)<-"F5"
Idents(G5)<-"G5"
Idents(H5)<-"H5"

Merged<-merge(A6, y = c(B6,C5,C6,C7,D5,D6,D7,E5,E7,F5,G5,H5), project = "Merged", merge.data = TRUE)


plot<-VlnPlot(Merged, features = c("nCount_RNA"),pt.size = 0)+stat_summary(fun.y = median, geom='point', size = 5, colour = "darkred", shape = 95)
VlnPlot(Merged, features = c("percent.mt"),pt.size = 0,y.max = 50)+stat_summary(fun.y = median, geom='point', size = 5, colour = "darkred", shape = 95)

setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/updates/2020-04-21/GCA")
png('A6.png')
FeatureScatter(A6, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

C7<-AddMetaData(C7,'GCA1',col.name = 'batch')
C5<-AddMetaData(C5,'GCA2',col.name = 'batch')
D5<-AddMetaData(D5,'GCA3',col.name = 'batch')

C7<-SCTransform(C7,verbose = F)
C5<-SCTransform(C5,verbose = F)
D5<-SCTransform(D5,verbose = F)


Merged<-merge(C7, y = c(C5,D5), project = "Merged", merge.data = TRUE)
Merged<-SCTransform(Merged,verbose = F)
Merged<-RunPCA(Merged,verbose = F)
Merged<-RunUMAP(Merged,verbose = F,dims = 1:30)
Merged<-FindNeighbors(Merged,dims = 1:30,verbose = F)
Merged<-FindClusters(Merged,verbose = F,resolution = 0.8)

SCT_merged<-merge(C7, y = c(C5,D5), project = "Merged", merge.data = TRUE)
SCT_merged<-RunPCA(SCT_merged,verbose = F)
SCT_merged<-RunUMAP(SCT_merged,verbose = F,dims = 1:30)
SCT_merged<-FindNeighbors(SCT_merged,dims = 1:30,verbose = F)
SCT_merged<-FindClusters(SCT_merged,verbose = F,resolution = 0.8)
