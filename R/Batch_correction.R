setwd("C:/Users/yjk12/Dropbox/Batch_correction/5GEX_PBMC1_36yoF_F1/")
library(magrittr)
library(dplyr)
library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)

Blood2<-Read10X(data.dir = 'Batch1', gene.column = 2, unique.features = TRUE)
Blood3<-Read10X(data.dir = 'Batch2', gene.column = 2, unique.features = TRUE)

#seurat_object<-CreateSeuratObject(counts = Blood1)
F1B1<-CreateSeuratObject(counts = Blood2$`Gene Expression`)
F1B1[['Protein']] = CreateAssayObject(counts = Blood2$`Antibody Capture`)

F1B2<-CreateSeuratObject(counts = Blood3$`Gene Expression`)
F1B2[['Protein']] = CreateAssayObject(counts = Blood3$`Antibody Capture`)

F1B1<-SCTransform(F1B1, verbose = T)
F1B2<-SCTransform(F1B2,verbose = T)
F1B1<-AddMetaData(F1B1,'B1',col.name = 'batch')
F1B2<-AddMetaData(F1B2,'B2',col.name = 'batch')
var.gene<-union(F1B1@assays$SCT@var.features,F1B2@assays$SCT@var.features)

Merged <- merge(F1B1, y = c(F1B2), project = "Merged", merge.data = TRUE)
Merged@assays$SCT@var.features<-var.gene
#F1B1<-RunPCA(F1B1,verbose = F)
#F1B2<-RunPCA(F1B2,verbose = F)
Merged<-RunPCA(Merged,verbose = F,assay = "SCT")
pancreas.list <- list('B1'=F1B1,'B2'=F1B2)
wilcox.p <- function(n){
  #print(VlnPlot(F1, features = paste0('PC_', n), group.by = "orig.ident"))
  t <- wilcox.test(Merged@reductions$pca@cell.embeddings[which(Merged$batch=="B1"), paste0('PC_', n)], 
                   Merged@reductions$pca@cell.embeddings[which(Merged$batch=="B2"), paste0('PC_', n)])
  return(-log(t$p.value))
}
x <- sapply(1:50, FUN = function(n) wilcox.p(n))
plot(x[1:50])



VlnPlot(Merged, "PC_1", group.by = "batch", cols = viridis::viridis(3))



pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)



### Identify Anchors
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
