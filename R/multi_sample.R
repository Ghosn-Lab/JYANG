setwd("F:/Ghosn_lab/ADT/Fetal Samples")
library(Seurat)
library(ggplot2)
library(SeuratData)

LH100<-Read10X(data.dir = 'LH28188_100dM', gene.column = 2, unique.features = TRUE)
SH100<-Read10X(data.dir = 'SH28188_100dM', gene.column = 2, unique.features = TRUE)

LH<-CreateSeuratObject(counts = LH100$`Gene Expression`)
LH[['Protein']] = CreateAssayObject(counts = LH100$`Antibody Capture`)

SH<-CreateSeuratObject(counts = SH100$`Gene Expression`)
SH[['Protein']] = CreateAssayObject(counts = SH100$`Antibody Capture`)

pancreas.list <- c(SH=SH,LH=LH)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


pancreas.list$SH$orig.ident<-'SH'
pancreas.list$LH$orig.ident<-'LH'

reference.list <- pancreas.list[c('SH','LH')]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
DimPlot(pancreas.integrated, reduction = "umap",group.by = 'orig.ident')



### Clustering
pancreas.integrated<-FindNeighbors(pancreas.integrated,dims = 1:30,verbose = FALSE)
pancreas.integrated<-FindClusters(pancreas.integrated,verbose = FALSE)
DimPlot(pancreas.integrated,label = TRUE)+NoLegend()






#########################################
#########################################


##SC Transform
pancreas.list <- c(SH=SH,LH=LH)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}
pancreas.list$SH$orig.ident<-'SH'
pancreas.list$LH$orig.ident<-'LH'

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, group.by = 'orig.ident', combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, 
                                                                                                              byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)


### Clustering
pancreas.integrated<-FindNeighbors(pancreas.integrated,dims = 1:30,verbose = FALSE)
pancreas.integrated<-FindClusters(pancreas.integrated,verbose = FALSE)
DimPlot(pancreas.integrated,label = TRUE)+NoLegend()