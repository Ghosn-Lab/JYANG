setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Devon/Adding sequence/")
library(magrittr)
library(dplyr)
library(Seurat)
library(ggplot2)

#Blood1<-Read10X(data.dir = 'BM2', gene.column = 2, unique.features = TRUE)
Blood2<-Read10X(data.dir = 'LungTotalnew_01', gene.column = 2, unique.features = TRUE)


setwd("C:/Users/yjk12/Box/Devon's cell-cell")
Blood3<-Read10X(data.dir = 'LungTotalnew_01', gene.column = 2, unique.features = TRUE)


#seurat_object<-CreateSeuratObject(counts = Blood1)
seurat_object<-CreateSeuratObject(counts = Blood2$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = Blood2$`Antibody Capture`)
seurat_object<-AddMetaData(seurat_object,'extra',col.name = 'batch')
origin<-CreateSeuratObject(counts = Blood3$`Gene Expression`)
origin[['Protein']] = CreateAssayObject(counts = Blood3$`Antibody Capture`)
origin<-AddMetaData(origin,'origin',col.name = 'batch')
Merged <- merge(seurat_object, y = c(origin), project = "Merged", merge.data = TRUE)

#seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
#seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_object <- SCTransform(seurat_object, verbose = FALSE)
#origin <- PercentageFeatureSet(origin, pattern = "^MT-", col.name = "percent.mt")
#origin <- SCTransform(origin, vars.to.regress = "percent.mt", verbose = FALSE)
origin<-SCTransform(origin,verbose = F)
Merged<-SCTransform(Merged,verbose = F)
##### Different Highly variable genes
extra<-setdiff(seurat_object@assays$SCT@var.features,origin@assays$SCT@var.features)
old<-setdiff(origin@assays$SCT@var.features,seurat_object@assays$SCT@var.features)
#####

seurat_object <- NormalizeData(seurat_object, assay = "Protein", normalization.method = "CLR")
seurat_object <- ScaleData(seurat_object, assay = "Protein")
origin <- NormalizeData(origin, assay = "Protein", normalization.method = "CLR")
origin <- ScaleData(origin, assay = "Protein")
Merged<-NormalizeData(Merged,assay="Protein",normalization.method = "CLR")
Merged<-NormalizeData(Merged,assay="Protein")

seurat_object<-RunPCA(seurat_object,verbose = FALSE)
seurat_object<-RunUMAP(seurat_object,dims = 1:30,verbose = FALSE)

origin<-RunPCA(origin,verbose = F)
origin<-RunUMAP(origin,dims = 1:30,verbose = F)

Merged<-RunPCA(Merged,verbose = F)
Merged<-RunUMAP(Merged,verbose = F,dims = 1:30)

seurat_object<-FindNeighbors(seurat_object,dims = 1:30,verbose = FALSE)
seurat_object<-FindClusters(seurat_object,verbose = FALSE,resolution =0.8)
DimPlot(seurat_object,label = TRUE)   #, cells.highlight = cell

origin<-FindNeighbors(origin,dims = 1:30,verbose = F)
origin<-FindClusters(origin,verbose = F,resolution = 0.8)
DimPlot(origin,label = T)

Merged<-FindNeighbors(Merged,dims = 1:30,verbose = F)
Merged<-FindClusters(Merged,verbose = F,resolution = 0.8)
DimPlot(Merged,label = T,split.by = "batch",group.by = "batch")



extra.marker<-FindAllMarkers(seurat_object,only.pos = TRUE,assay = "SCT")
origin.marker<-FindAllMarkers(origin,only.pos = TRUE,assay = "SCT")

extra.marker %>% group_by(cluster) %>% top_n(n = 20) %>%write.table(file = 'DEG_total.csv' ,sep = ',',col.names=TRUE,row.names = FALSE)




######Dot plot
FeaturePlot(seurat_object, features = c( "PR8-NA-AF38912-1","tdTomato",'PR8-rc',"Tcf7"), max.cutoff = "q95", ncol = 2,label = F)

DotPlot(seurat_object, features = c( "tdTomato",'PR8-rc',"PR8-M-AF389121-1","PR8-NA-AF38912-1","PR8-NP-AF38911-1",
                                     "PR8-NS-AF38912-1","PR8-PA-AF38911-1","PR8-PB1-AF3891-1","PR8-PB2-AF3891-1"),
        group.by = "seurat_clusters",col.min = 0,dot.scale = 4) + RotatedAxis()


ADT<-rownames(seurat_object@assays$Protein)
DotPlot(seurat_object, features = ADT, group.by = "seurat_clusters",assay = 'Protein',col.min = 0,dot.scale = 4) + RotatedAxis()



DotPlot(seurat_object,features = head(extra.marker[which(extra.marker$cluster==6),"gene"],n=20),group.by = "seurat_clusters", assay = "SCT", col.min = 0,dot.scale = 4)+RotatedAxis()





### Highly expressed gene
######################################################
######################################################
cluster5.markers <- FindMarkers(origin, ident.1 = 5, min.pct = 0.1,only.pos = TRUE)
ori_hi_5<-rownames(head(cluster5.markers, n = 20))
head(cluster5.markers, n = 20)

cluster6.markers <- FindMarkers(origin, ident.1 = 6, min.pct = 0.1,only.pos = TRUE)
ori_hi_6<-rownames(head(cluster6.markers, n = 20))
head(cluster6.markers, n = 20)

cluster0.markers <- FindMarkers(origin, ident.1 = 0, min.pct = 0.1,only.pos = TRUE)
ori_hi_0<-rownames(head(cluster0.markers, n = 20))
head(cluster0.markers, n = 20)

### Extra cluster
cluster1.markers <- FindMarkers(origin, ident.1 = 1,ident.2 = 5 ,min.pct = 0.1,only.pos = TRUE)
ori_hi_1<-rownames(head(cluster1.markers, n = 20))
head(cluster1.markers, n = 20)

add_cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, min.pct = 0.1,only.pos = TRUE)
add_hi_0<-rownames(head(add_cluster0.markers, n = 20))
head(add_cluster0.markers, n = 20)

add_cluster1.markers <- FindMarkers(seurat_object, ident.1 = 1, min.pct = 0.1,only.pos = TRUE)
add_hi_1<-rownames(head(add_cluster1.markers, n = 20))
head(add_cluster1.markers, n = 20)

add_cluster5.markers <- FindMarkers(seurat_object, ident.1 = 5, min.pct = 0.1,only.pos = TRUE)
add_hi_5<-rownames(head(add_cluster5.markers, n = 20))
head(add_cluster5.markers, n = 20)

###################################################
###################################################
DefaultAssay(origin)<-"SCT"
oldmarker5_1<-FindMarkers(origin,ident.1 = 5,ident.2 = 1)
oldmarker5_1<-oldmarker5_1[which(oldmarker5_1$p_val_adj<0.05),]
FeaturePlot(Merged,features = c("Trbc2","Ighm","Trac","Trbc1"),label = T)
FeaturePlot(Merged, features = c( "PR8-NP-AF38911-1","PR8-M-AF389121-1",'PR8-rc',"PR8-NS-AF38912-1"), max.cutoff = "q95", ncol = 2,label = F)

marker0_1<-FindMarkers(Merged,ident.1 = 0,ident.2 = 1)
marker5_4<-FindMarkers(Merged,ident.1 = 5,ident.2 = 4)
marker5_4<-marker5_4[which(marker5_4$p_val_adj<0.05),]

a<-colnames(Merged@assays$SCT@counts[,which(Merged$batch=="extra")])
b<-colnames(Merged@assays$SCT@counts[,which(Merged$batch=="origin")])
Idents(Merged,cells=a)<-"extra"
Idents(Merged,cells=b)<-"origin"
c<-FindMarkers(Merged,ident.1 = "extra",ident.2 = "origin")

new.cluster.ids <- c("Naive CD4 T", "IFN-induced Th1", "Macrophages", "CD8 T", "Activated CD4 T", "Exhausted T","Lag3+ CD4 T",
                      "alvM¦µ", "B","Neutrophils","T reg","NKT","IAV_induced Proliferating T")
names(new.cluster.ids) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
