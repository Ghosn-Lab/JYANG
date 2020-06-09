setwd('F:/Ghosn_lab/AWS/force_cell_compare')
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
options(future.globals.maxSize = 4000 * 1024^2)
### Read data
sample1<-Read10X(data.dir = 'F:/Ghosn_lab/AWS/force_cell_compare/F1/original', gene.column = 2, unique.features = TRUE)
original<-CreateSeuratObject(counts = sample1$`Gene Expression`)
original[['Protein']] = CreateAssayObject(counts = sample1$`Antibody Capture`)
original<-AddMetaData(original,'original',col.name = 'batch')



sample2<-Read10X(data.dir = 'F:/Ghosn_lab/AWS/force_cell_compare/F1/4000cells', gene.column = 2, unique.features = TRUE)
modify<-CreateSeuratObject(counts = sample2$`Gene Expression`)
modify[['Protein']] = CreateAssayObject(counts = sample2$`Antibody Capture`)
modify<-AddMetaData(modify,'4000',col.name = 'batch')




pancreas.list <- list('original'=original,'modify'=modify)
### Normalization
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)



### Identify Anchors
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

### Dimensionality reduction

pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, combine = FALSE,group.by = 'batch',split.by = 'batch')
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)

### Clustering
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30, verbose = FALSE)
pancreas.integrated <- FindClusters(pancreas.integrated, verbose = FALSE)
DimPlot(pancreas.integrated, label = TRUE,split.by = 'batch') + NoLegend()
D<-setdiff(colnames(sample2$`Gene Expression`),colnames(sample1$`Gene Expression`))
for (i in 1:length(D))
  D[i]=paste(D[i],'_2',sep = '')

DimPlot(pancreas.integrated, label = TRUE,cells.highlight = D) + NoLegend()# Localize extra cells

### Markers
markers <- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay='SCT')
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) %>%write.table(file = 'F1_diff.csv' ,sep = ',',col.names=TRUE,row.names = FALSE)

marker.adt<-FindAllMarkers(pancreas.integrated,only.pos = TRUE,assay = 'Protein')
marker.adt %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC) %>%write.table(file = 'F1_diff_adt.csv' ,sep = ',',col.names=TRUE,row.names = FALSE)
head(markers[markers['cluster']==0,])      ##specify a cluster

##Feature plot
rownames(pancreas.integrated@assays$Protein)
pancreas.integrated <- NormalizeData(pancreas.integrated, assay = "Protein", normalization.method = "CLR")
pancreas.integrated <- ScaleData(pancreas.integrated, assay = "Protein")
#B
FeaturePlot(pancreas.integrated, features = c( "CD19-TotalSeqC","CD20-TotalSeqC"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
VlnPlot(pancreas.integrated, features = c( "CD19-TotalSeqC","CD20-TotalSeqC"),group.by = "seurat_clusters",pt.size=0,split.by = "batch")

#T
FeaturePlot(pancreas.integrated, features = c( "CD4-TotalSeqC","CD8a-TotalSeqC","CD3-TotalSeqC","CD127-TotalSeqC"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
#Monocyte
FeaturePlot(pancreas.integrated, features = c( "CD16-TotalSeqC","CD45-TotalSeqC","CD56-TotalSeqC","CD19-TotalSeqC"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)


### Dot plot
DotPlot(pancreas.integrated, features = c( "CD19-TotalSeqC","CD20-TotalSeqC","CD34-TotalSeqC","CD27-TotalSeqC","CD4-TotalSeqC","CD8a-TotalSeqC",
                                           "CD3-TotalSeqC","CD127-TotalSeqC","CD16-TotalSeqC","CD45-TotalSeqC","CD56-TotalSeqC","CD25-TotalSeqC","CD14-TotalSeqC"),
                                          group.by = "seurat_clusters",assay = 'Protein',split.by = "batch",col.min = 0,c( "red", "blue"),dot.scale = 4) + RotatedAxis()


DoHeatmap(pancreas.integrated,features = c( "CD19-TotalSeqC","CD20-TotalSeqC","CD34-TotalSeqC","CD27-TotalSeqC","CD4-TotalSeqC","CD8a-TotalSeqC",
                                            "CD3-TotalSeqC","CD127-TotalSeqC","CD16-TotalSeqC","CD45-TotalSeqC","CD56-TotalSeqC","CD25-TotalSeqC","CD14-TotalSeqC") ,
                                            assay = 'Protein',size = 3,  angle = 90) + NoLegend()


VlnPlot(pancreas.integrated, features = c( "CD16-TotalSeqC","CD25-TotalSeqC"),group.by = "seurat_clusters",pt.size=0,split.by = "batch")

cluster20.markers <- FindMarkers(pancreas.integrated, ident.1 = 20, min.pct = 0.25,only.pos = TRUE)
hi<-rownames(head(cluster20.markers, n = 20))

DotPlot(pancreas.integrated, features = hi,group.by = "seurat_clusters",assay = 'RNA',col.min = 0,dot.scale = 4,split.by = 'batch',c( "red", "blue")) + RotatedAxis()
