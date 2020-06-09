setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Devon/Adding sequence/")
setwd("F:/Ghosn_lab/Devon_test/Devon_R/Devon_R/Devon/")
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Troy/")
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/Cite_seq/")
library(Seurat)
library(ggplot2)

#Blood1<-Read10X(data.dir = 'BM2', gene.column = 2, unique.features = TRUE)
Blood2<-Read10X(data.dir = 'LungTotalnew_01', gene.column = 2, unique.features = TRUE)
Blood3<-Read10X(data.dir = 'Total', gene.column = 2, unique.features = TRUE)

#seurat_object<-CreateSeuratObject(counts = Blood1)
seurat_object<-CreateSeuratObject(counts = Blood2$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = Blood2$`Antibody Capture`)

origin<-CreateSeuratObject(counts = Blood3$`Gene Expression`)
origin[['Protein']] = CreateAssayObject(counts = Blood3$`Antibody Capture`)

##Normalize RNA data ---SC transform
#seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
#seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_object <- SCTransform(seurat_object, verbose = FALSE)
origin<-SCTransform(origin)
#Log normalization
seurat_object <- NormalizeData(seurat_object,verbose = T)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
##Normalize ADT data ---CLR transform
seurat_object <- NormalizeData(seurat_object, assay = "Protein", normalization.method = "CLR")
seurat_object <- ScaleData(seurat_object, assay = "Protein")
origin <- NormalizeData(origin, assay = "Protein", normalization.method = "CLR")
origin <- ScaleData(origin, assay = "Protein")

##Dimensionality reduction AND Clustering using RNA

seurat_object<-RunPCA(seurat_object,verbose = FALSE)
seurat_object<-RunUMAP(seurat_object,dims = 1:30,verbose = FALSE)

origin<-RunPCA(origin,verbose = F)
origin<-RunUMAP(origin,dims = 1:30,verbose = F)
seurat_object <- RunTSNE(seurat_object, dims = 1:25, method = "FIt-SNE")

seurat_object<-FindNeighbors(seurat_object,dims = 1:30,verbose = FALSE)
seurat_object<-FindClusters(seurat_object,verbose = FALSE,resolution =0.8)
DimPlot(seurat_object,label = TRUE)   #, cells.highlight = cell
origin<-FindNeighbors(origin,dims = 1:30,verbose = F)
origin<-FindClusters(origin,verbose = F,resolution = 0.8)
DimPlot(origin,label = T)

cell=c("AACTGGTCAGTCAGAG","ACGGAGACAGAAGCAC","CACATAGAGGACCACA","CAGTAACCAGGAACGT","CGGTTAAAGAGGACGG","CTTTGCGCACAACTGT","GCTCCTACAGCCAATT","GTTAAGCGTTGTACAC","TATGCCCCAGACAAAT","TTTACTGCATCCCACT")

###Find markders that defefine clusters
seurat.rna.markers <- FindAllMarkers(seurat_object, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)

new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs",'a','b','c')
names(new.cluster.ids) <- levels(seurat_object)

##Feature plot
FeaturePlot(seurat_object, features = c( "PR8-NA-AF38912-1","tdTomato",'PR8-rc',"PR8-M-AF389121-1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
RidgePlot(seurat_object, features = c('CD8a-TotalSeqC','CD34-TotalSeqC'), ncol = 2)
FeatureScatter(seurat_object, feature1 = "CD8a-TotalSeqC", feature2 = "CD8A")

##DE of ADT
seurat.small <- subset(seurat_object, downsample = 300)
seurat.markers <- FindAllMarkers(seurat.small, assay = "Protein", only.pos = TRUE)
DoHeatmap(seurat.small, features = unique(seurat.markers$gene), assay = "Protein", angle = 90) + NoLegend()


DotPlot(seurat_object, features = c( "tdTomato.dna-rc",'PR8-rc',"PR8-M-AF389121-1","PR8-NA-AF38912-1","PR8-NP-AF38911-1",
                                     "PR8-NS-AF38912-1","PR8-PA-AF38911-1","PR8-PB1-AF3891-1","PR8-PB2-AF3891-1"),
                                    group.by = "seurat_clusters",col.min = 0,dot.scale = 4) + RotatedAxis()

adt<-rownames(seurat_object@assays$Protein)
DotPlot(seurat_object, features = adt,group.by = "seurat_clusters",assay = 'Protein',col.min = 0,dot.scale = 4) + RotatedAxis()


cluster5.markers <- FindMarkers(seurat_object, ident.1 = 8, min.pct = 0.25,only.pos = TRUE)
hi<-rownames(head(cluster5.markers, n = 20))
head(cluster5.markers, n = 20)
DotPlot(seurat_object, features = hi,group.by = "seurat_clusters",assay = 'RNA',col.min = 0,dot.scale = 4) + RotatedAxis()
#######################################
#######################################
##Dimensionality reduction AND Clustering using ADT
DefaultAssay(seurat_object) <- "Protein"
seurat_object <- RunPCA(seurat_object, features = rownames(seurat_object), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
DimPlot(seurat_object, reduction = "pca_adt")


FeaturePlot(seurat_object, features = c( "YFP"), ncol = 2)
FeaturePlot(seurat_object, features = c( "Itgax","Cybb","Sirpa","Cx3cr1"), min.cutoff = "q05", max.cutoff = "q95", ncol = 2)
