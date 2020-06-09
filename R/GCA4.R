library(Seurat)
library(qusage)
library(ggplot2)
library(Matrix)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/filtered/")



########## Log normalization
filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA4_CD45neg_21461_5GEX_G5",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"1",sep = "_")
C7<-CreateSeuratObject(counts = filtered)
rm(filtered)

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA5_CD45neg_21478_5GEX_E7",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"2",sep = "_")
C5<-CreateSeuratObject(counts = filtered)
rm(filtered)

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA6_CD45neg_20006_5GEX_D6",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"3",sep = "_")
D5<-CreateSeuratObject(counts = filtered)
rm(filtered)

C7<-AddMetaData(C7,'GCA4_total',col.name = 'batch')
C5<-AddMetaData(C5,'GCA4_CD45pos',col.name = 'batch')
D5<-AddMetaData(D5,'GCA4_CD45neg',col.name = 'batch')

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

###var scale Merge 
###
protease<-merge(C7, y = c(C5,D5), project = "protease", merge.data = TRUE)
protease<-SCTransform(protease,verbose = F)
###var genes prior merge
protease@assays$RNA@var.features<-Uv
Mscale<-cbind(C7@assays$RNA@scale.data[protease@assays$RNA@var.features,],C5@assays$RNA@scale.data[protease@assays$RNA@var.features,])
Mscale<-cbind(Mscale,D5@assays$RNA@scale.data[protease@assays$RNA@var.features,])

protease@assays$RNA@scale.data<-Mscale
rm(Mscale)

protease<-RunPCA(protease,verbose = F)
protease<-RunUMAP(protease,verbose = F,dims = 1:30)
protease<-FindNeighbors(protease,dims = 1:30,verbose = F)
protease<-FindClusters(protease,verbose = F,resolution = 0.8)

DimPlot(protease,group.by = "batch",cols = c("green","yellow","red"))



C7<-SCTransform(C7,verbose = F)
C7<-RunPCA(C7,verbose = F)
C7<-RunUMAP(C7,verbose = F,dims = 1:30)
C7<-FindNeighbors(C7,dims = 1:30,verbose = F)
C7<-FindClusters(C7,verbose = F,resolution = 0.8)

DimPlot(C7,label = T)

C5<-SCTransform(C5,verbose = F)
C5<-RunPCA(C5,verbose = F)
C5<-RunUMAP(C5,verbose = F,dims = 1:30)
C5<-FindNeighbors(C5,dims = 1:30,verbose = F)
C5<-FindClusters(C5,verbose = F,resolution = 0.8)

DimPlot(C5,label = T)

D5<-SCTransform(D5,verbose = F)
D5<-RunPCA(D5,verbose = F)
D5<-RunUMAP(D5,verbose = F,dims = 1:30)
D5<-FindNeighbors(D5,dims = 1:30,verbose = F)
D5<-FindClusters(D5,verbose = F,resolution = 0.8)

DimPlot(D5,label = T)


c17<-FindMarkers(D5,ident.1 = 8,only.pos = T)
DotPlot(D5,features = head(rownames(c17),n=20),col.min = 0,dot.scale = 4)+RotatedAxis()

c7_7<-names(C7$seurat_clusters[C7$seurat_clusters==7])
c7_2<-names(C7$seurat_clusters[C7$seurat_clusters==2])
c7_6<-names(C7$seurat_clusters[C7$seurat_clusters==6])
c7_5<-names(C7$seurat_clusters[C7$seurat_clusters==5])
c7_DCs<-names(C7$seurat_clusters[C7$seurat_clusters==8])


c7_B<-c(c7_2,c7_6)
c7_T<-c(c7_5)
Non_c7<-colnames(C7)
Non_c7<-setdiff(Non_c7,c(c7_B,c7_5,c7_DCs,c7_7))



c5_3<-names(C5$seurat_clusters[C5$seurat_clusters==3])
c5_11<-names(C5$seurat_clusters[C5$seurat_clusters==11])
c5_DCs<-names(C5$seurat_clusters[C5$seurat_clusters==6])

c5_T<-c(c5_3,c5_11)
Non_c5<-colnames(C5)
Non_c5<-setdiff(Non_c5,c(c5_T,c5_DCs))


d5_5<-names(D5$seurat_clusters[D5$seurat_clusters==5])
d5_11<-names(D5$seurat_clusters[D5$seurat_clusters==11])
d5_4<-names(D5$seurat_clusters[D5$seurat_clusters==4])
d5_7<-names(D5$seurat_clusters[D5$seurat_clusters==7])

d5_T<-c(d5_5,d5_11)

Non_d5<-colnames(D5)
Non_d5<-setdiff(Non_d5,c(d5_T))


p_5<-names(D5$seurat_clusters[D5$seurat_clusters==5])
p_11<-names(D5$seurat_clusters[D5$seurat_clusters==11])
p_4<-names(D5$seurat_clusters[D5$seurat_clusters==4])
p_7<-names(D5$seurat_clusters[D5$seurat_clusters==7])

DimPlot(protease,cells.highlight = list("GCA1_B"=c7_B),cols.highlight = c("indianred1"))
DimPlot(protease,cells.highlight = list("GCA1_T"=c7_T,"GCA2_T"=c5_T,"GCA3_T"=d5_T),cols.highlight = c("dodgerblue1","palegreen3","indianred1"))
DimPlot(protease,cells.highlight = list("GCA1_DCs"=c7_DCs,"GCA2_DCs"=c5_DCs,"GCA3_DCs"=d5_DCs),cols.highlight = c("dodgerblue1","palegreen3","indianred1"))
DimPlot(protease,cells.highlight = list("GCA1_Non_immune"=Non_c7,"GCA2_Non_immune"=Non_c5,"GCA3_Non_immune"=Non_d5),cols.highlight = c("dodgerblue1","palegreen3","indianred1"))

DotPlot(protease,features = head(rownames(a),n=20),col.min = 0,dot.scale = 4)+RotatedAxis()
protease[["percent.mt"]] <- PercentageFeatureSet(protease, pattern = "^MT-")
