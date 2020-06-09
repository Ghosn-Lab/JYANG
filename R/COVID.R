library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)
library(dsb)
### filtered matrix
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/COVID-19/mito/Sputum4102/")

filtered<-Read10X("full_mito",gene.column = 2, unique.features = TRUE)
# colnames(filtered$`Gene Expression`)<-sub("-1$","_1",colnames(filtered$`Gene Expression`))
# colnames(filtered$`Antibody Capture`)<-sub("-1$","_1",colnames(filtered$`Antibody Capture`))
E7<-CreateSeuratObject(counts = filtered$`Gene Expression`)
E7[['Protein']] = CreateAssayObject(counts = filtered$`Antibody Capture`)
rm(filtered)

a<-Read10X("../filtered/Sputum_4102",gene.column = 2, unique.features = TRUE)
extra_cell<-setdiff(colnames(E7),colnames(a$`Gene Expression`))

pbmc.htos<-Read10X("Sputum_4102",gene.column = 2, unique.features = TRUE)
hto<-CreateSeuratObject(counts = pbmc.htos$`Gene Expression`)
hto[['Protein']]<-CreateAssayObject(counts = pbmc.htos$`Antibody Capture`)
E7<-SubsetData(hto,subset.name = "nFeature_RNA",low.threshold = 40)
E7<-SubsetData(E7,subset.name = "nCount_RNA",low.threshold = 100)
rm(pbmc.htos)
real<-SubsetData(hto,subset.name = "nFeature_RNA",low.threshold = 40)
real<-SubsetData(real,subset.name = "nCount_RNA",low.threshold = 100)
real<-as.matrix(real@assays$Protein@counts)
background<-SubsetData(hto,subset.name = "nFeature_RNA",high.threshold = 40)
background<-as.matrix(background@assays$Protein@counts)

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = real, empty_drop_matrix = background)
E7[['Protein']] = CreateAssayObject(counts = normalized_matrix)
E7[['Protein']] = CreateAssayObject(counts = real)
real<-real[which(rowSums(real)>0),]

filtered<-Read10X("Sputum_4106",gene.column = 2, unique.features = TRUE)
# colnames(filtered$`Gene Expression`)<-sub("-1$","_2",colnames(filtered$`Gene Expression`))
# colnames(filtered$`Antibody Capture`)<-sub("-1$","_2",colnames(filtered$`Antibody Capture`))
A6<-CreateSeuratObject(counts = filtered$`Gene Expression`)
A6[['Protein']] = CreateAssayObject(counts = filtered$`Antibody Capture`)
rm(filtered)

b<-Read10X("../filtered/Sputum_4106",gene.column = 2, unique.features = TRUE)

pbmc.htos<-Read10X("Sputum_4106",gene.column = 2, unique.features = TRUE)
hto<-CreateSeuratObject(counts = pbmc.htos$`Gene Expression`)
hto[['Protein']]<-CreateAssayObject(counts = pbmc.htos$`Antibody Capture`)
A6<-SubsetData(hto,subset.name = "nFeature_RNA",low.threshold = 40)
A6<-SubsetData(A6,subset.name = "nCount_RNA",low.threshold = 100)
rm(pbmc.htos)
real<-SubsetData(hto,subset.name = "nFeature_RNA",low.threshold = 40)
real<-SubsetData(real,subset.name = "nCount_RNA",low.threshold = 100)
real<-as.matrix(real@assays$Protein@counts)
background<-SubsetData(hto,subset.name = "nFeature_RNA",high.threshold = 40)
background<-as.matrix(background@assays$Protein@counts)

normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = real, empty_drop_matrix = background)
A6[['Protein']] = CreateAssayObject(counts = normalized_matrix)
A6[['Protein']] = CreateAssayObject(counts = real)
real<-real[which(rowSums(real)>0),]
background<-background[which(rowSums(background)>0),]

filtered<-Read10X("Blood_4102",gene.column = 2, unique.features = TRUE)
# colnames(filtered$`Gene Expression`)<-sub("-1$","_3",colnames(filtered$`Gene Expression`))
# colnames(filtered$`Antibody Capture`)<-sub("-1$","_3",colnames(filtered$`Antibody Capture`))
D7<-CreateSeuratObject(counts = filtered$`Gene Expression`)
rm(filtered)

filtered<-Read10X("Blood_4106",gene.column = 2, unique.features = TRUE)
# colnames(filtered$`Gene Expression`)<-sub("-1$","_4",colnames(filtered$`Gene Expression`))
# colnames(filtered$`Antibody Capture`)<-sub("-1$","_4",colnames(filtered$`Antibody Capture`))
H5<-CreateSeuratObject(counts = filtered$`Gene Expression`)
rm(filtered)

E7<-SubsetData(E7,subset.name = "nFeature_RNA",low.threshold = 40)
E7<-SubsetData(E7,subset.name = "nCount_RNA",low.threshold = 100)
A6<-SubsetData(A6,subset.name = "nFeature_RNA",low.threshold = 40)
A6<-SubsetData(A6,subset.name = "nCount_RNA",low.threshold = 100)
D7<-SubsetData(D7,subset.name = "nFeature_RNA",low.threshold = 40)
D7<-SubsetData(D7,subset.name = "nCount_RNA",low.threshold = 100)
H5<-SubsetData(H5,subset.name = "nFeature_RNA",low.threshold = 40)
H5<-SubsetData(H5,subset.name = "nCount_RNA",low.threshold = 100)

E7<-AddMetaData(E7,'Sputum_4102',col.name = 'batch')
A6<-AddMetaData(A6,'Sputum_4106+',col.name = 'batch')
H5<-AddMetaData(H5,'Blood_4102',col.name = 'batch')
D7<-AddMetaData(D7,'Blood_4102',col.name = 'batch')


Collagenase<-merge(D7, y = c(E7,A6), project = "Collagenase", merge.data = TRUE)

Merged<-merge(H5, y = c(D7,E7,A6), project = "Merged", merge.data = TRUE)
a<-Merged$orig.ident
levels(a)<-c("Protease","Collagenase")
Pro<-rownames(Merged@meta.data[which(Merged@meta.data["batch"]=="Protease"),])
Colla<-rownames(Merged@meta.data[which(Merged@meta.data["batch"]!="Protease"),])
a[Pro]<-"Protease"
a[Colla]<-"Collagenase"
Merged$orig.ident<-a

Collagenase<-SCTransform(Collagenase,verbose = F)
Collagenase<-RunPCA(Collagenase,verbose = F)
Collagenase<-RunUMAP(Collagenase,verbose = F,dims = 1:30)
Collagenase<-FindNeighbors(Collagenase,dims = 1:30,verbose = F)
Collagenase<-FindClusters(Collagenase,verbose = F,resolution = 0.8)

Merged<-SCTransform(Merged,verbose = F,return.only.var.genes = F)
Merged<-RunPCA(Merged,verbose = F)
Merged<-RunUMAP(Merged,verbose = F,dims = 1:30)
Merged<-FindNeighbors(Merged,dims = 1:30,verbose = F)
Merged<-FindClusters(Merged,verbose = F,resolution = 0.8)




D7<-SCTransform(D7,verbose = F)
ev<-E7@assays$SCT@var.features
av<-A6@assays$SCT@var.features
dv<-D7@assays$SCT@var.features
Cov<-unique(c(ev,av,dv))

####
E7[["percent.mt"]] <- PercentageFeatureSet(E7, pattern = "^MT-")
E7<-subset(E7, subset = nFeature_RNA > 40 & nCount_RNA >100)
E7_filtered<-subset(E7, subset = percent.mt < 20)
E7 <- NormalizeData(E7, assay = "Protein", normalization.method = "CLR")
E7 <- ScaleData(E7, assay = "Protein")
E7<-SCTransform(E7,verbose = F)
E7<-RunPCA(E7,verbose = F)
E7<-RunUMAP(E7,verbose = F,dims = 1:30)
E7<-FindNeighbors(E7,dims = 1:30,verbose = F)
E7<-FindClusters(E7,verbose = F,resolution = 0.8)
DimPlot(E7,label = T)
DefaultAssay(E7)<-"RNA"
E7<-NormalizeData(E7,normalization.method = "LogNormalize", scale.factor = 10000)
FeaturePlot(E7,features = c("CD3-TotalA","CD19-TotalA","CD20-TotalA","percent.mt"),label = T)
FeaturePlot(E7,features = c("CD3D","CD19","MS4A1","CD4"),label = T)
FeaturePlot(E7,features = c("ACE2","TMPRSS2","TMPRSS4","CTSL"),label = T) ### bronchial epithelial
FeaturePlot(E7,features = c("CTSL","DPP4","ANPEP","ACE"),label = T) ###


####

A6[["percent.mt"]] <- PercentageFeatureSet(A6, pattern = "^MT-")
A6<-subset(A6, subset = nFeature_RNA > 40 & nCount_RNA >100)
A6_filtered<-subset(A6, subset = percent.mt < 20)
A6 <- NormalizeData(A6, assay = "Protein", normalization.method = "CLR")
A6 <- ScaleData(A6, assay = "Protein")
A6<-SCTransform(A6,verbose = F)
A6<-RunPCA(A6,verbose = F)
A6<-RunUMAP(A6,verbose = F,dims = 1:30)
A6<-FindNeighbors(A6,dims = 1:30,verbose = F)
A6<-FindClusters(A6,verbose = F,resolution = 0.8)
DimPlot(A6,label = T)

DefaultAssay(A6)<-"RNA"
A6<-NormalizeData(A6,normalization.method = "LogNormalize", scale.factor = 10000)
FeaturePlot(A6,features = c("CD3-TotalA","CD19-TotalA","CD20-TotalA","percent.mt"),label = T)
FeaturePlot(A6,features = c("CD3D","CD19","MS4A1","CD4"),label = T)
FeaturePlot(A6,features = c("ACE2","TMPRSS2","TMPRSS4","CTSL"),label = T) ### bronchial epithelial
FeaturePlot(A6,features = c("CTSL","DPP4","ANPEP","ACE"),label = T) ###


a<-as.matrix(tail(A6@assays$RNA@counts,n=13))
rowSums(a)


