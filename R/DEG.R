setwd("C:/Users/yjk12/Dropbox/Batch_correction/5GEX_PBMC1_36yoF_F1/")
library(magrittr)
library(dplyr)
library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)


Blood3<-Read10X(data.dir = 'Batch2', gene.column = 2, unique.features = TRUE)

#seurat_object<-CreateSeuratObject(counts = Blood1)


F1B2<-CreateSeuratObject(counts = Blood3$`Gene Expression`)
F1B2[['Protein']] = CreateAssayObject(counts = Blood3$`Antibody Capture`)
F1B2 <- NormalizeData(F1B2, assay = "Protein", normalization.method = "CLR")
F1B2 <- ScaleData(F1B2, assay = "Protein")

F1B2<-SCTransform(F1B2,verbose = T)

#F1B1<-RunPCA(F1B1,verbose = F)
F1B2<-RunPCA(F1B2,verbose = F)
F1B2<-RunUMAP(F1B2,dims = 1:30,verbose = FALSE)
F1B2<-FindNeighbors(F1B2,dims = 1:30,verbose = FALSE)
F1B2<-FindClusters(F1B2,verbose = FALSE,resolution =0.8)
DimPlot(F1B2,label = TRUE)   #, cells.highlight = cell


FeaturePlot(F1B2, features = c( "CD8a-TotalSeqC","CD4-TotalSeqC"), ncol = 2,label = T)

### Rank Sum test
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/JC/DEG")
DefaultAssay(F1B2)<-"SCT"
marker9.wilcox<-FindMarkers(F1B2,ident.1 = 9,only.pos = T)
marker9.wilcox<-marker9.wilcox[which(marker9.wilcox['p_val_adj']<0.05),]
DotPlot(F1B2,features = rownames(marker9.wilcox),col.min = 0,dot.scale = 4)+RotatedAxis()

marker9_CD4<-FindMarkers(F1B2,ident.1 = 9,ident.2 = c(0,1,14,16,19),only.pos = T)
marker9_CD4<-marker9_CD4[which(marker9_CD4['p_val_adj']<0.05),]
DotPlot(F1B2,features = rownames(marker9_CD4),col.min = 0,dot.scale = 4)+RotatedAxis()
write.table(marker9.wilcox,file = "")

### MAST
DefaultAssay(F1B2)<-"RNA"
F1B2 <- NormalizeData(F1B2, assay = "RNA")
#F1B2 <- FindVariableFeatures(F1B2, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(F1B2)
#F1B2 <- ScaleData(F1B2, assay = "RNA",features = all.genes)
marker9.MAST<-FindMarkers(F1B2,ident.1 = 9,only.pos = T,test.use = "MAST",assay = "RNA")
marker9.MAST<-marker9.MAST[which(marker9.MAST['p_val_adj']<0.05),]
DotPlot(F1B2,features = rownames(marker9.MAST),col.min = 0,dot.scale = 4)+RotatedAxis()

marker9_CD4<-FindMarkers(F1B2,ident.1 = 9,ident.2 = c(0,1,14,16,19),only.pos = T,test.use = "MAST",assay = "RNA")
marker9_CD4<-marker9_CD4[which(marker9_CD4['p_val_adj']<0.05),]
DotPlot(F1B2,features = rownames(marker9_CD4),col.min = 0,dot.scale = 4)+RotatedAxis()



### DESeq2
DefaultAssay(F1B2)<-"RNA"
F1B2 <- NormalizeData(F1B2, assay = "RNA")
#F1B2 <- FindVariableFeatures(F1B2, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(F1B2)
#F1B2 <- ScaleData(F1B2, assay = "RNA",features = all.genes)
marker9.DESeq2<-FindMarkers(F1B2,ident.1 = 9,only.pos = T,test.use = "DESeq2",assay = "RNA",max.cells.per.ident = 100)
marker9.DESeq2<-marker9.DESeq2[which(marker9.DESeq2['p_val_adj']<0.05),]
DotPlot(F1B2,features = rownames(marker9.DESeq2),col.min = 0,dot.scale = 4)+RotatedAxis()