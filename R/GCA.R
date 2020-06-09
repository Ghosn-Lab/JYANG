library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)
### raw matrix
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/")
raw<-Read10X("GCA4_CD45neg_21461_5GEX_G5",gene.column = 2, unique.features = TRUE)
raw_matrix<-CreateSeuratObject(counts = raw)
rm(raw)
### filtered matrix
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/filtered/")
filtered<-Read10X("GCA4_CD45neg_21461_5GEX_G5",gene.column = 2, unique.features = TRUE)
filtered_matrix<-CreateSeuratObject(counts = filtered)
rm(filtered)

### select "real" cells from raw
raw_matrix<-SubsetData(raw_matrix,subset.name = "nFeature_RNA",low.threshold = 40)
raw_matrix<-SubsetData(raw_matrix,subset.name = "nCount_RNA",low.threshold = 100)


raw_matrix<-AddMetaData(raw_matrix,'raw',col.name = 'batch')
filtered_matrix<-AddMetaData(filtered_matrix,'filtered',col.name = 'batch')

raw_matrix<-SCTransform(raw_matrix,verbose = F)
filtered_matrix<-SCTransform(filtered_matrix,verbose = F)


raw_matrix<-RunPCA(raw_matrix,verbose = FALSE)
raw_matrix<-RunUMAP(raw_matrix,dims = 1:30,verbose = FALSE)
raw_matrix<-FindNeighbors(raw_matrix,dims = 1:30,verbose = FALSE)
raw_matrix<-FindClusters(raw_matrix,verbose = FALSE,resolution =0.8)

filtered_matrix<-RunPCA(filtered_matrix,verbose = FALSE)
filtered_matrix<-RunUMAP(filtered_matrix,dims = 1:30,verbose = FALSE)
filtered_matrix<-FindNeighbors(filtered_matrix,dims = 1:30,verbose = FALSE)
filtered_matrix<-FindClusters(filtered_matrix,verbose = FALSE,resolution =0.8)


extra_cells<-setdiff(names(raw_matrix$orig.ident),names(filtered_matrix$orig.ident))
a<-raw_matrix$orig.ident
levels(a)<-c("filtered","extra")
a[extra_cells]<-"extra"
raw_matrix$orig.ident<-a

DimPlot(raw_matrix,label = T,cells.highlight = extra_cells)
DimPlot(raw_matrix,split.by = "orig.ident",group.by = "seurat_clusters",label = T)

setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/")
save(filtered_matrix, raw_matrix, file="GCA6_total_20006_5GEX_B6.RData")


