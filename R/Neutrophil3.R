library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/")
options(future.globals.maxSize = 4000 * 1024^2)

gene.list<-read.gmt("Geneset-Library-Granulocytes.gmt")
PBMC1_origin<-Read10X(data.dir = 'PBMC4', gene.column = 2, unique.features = TRUE)
PBMC1<-Read10X(data.dir = 'BM2_4000', gene.column = 2, unique.features = TRUE)
PBMC2_origin<-Read10X(data.dir = 'PBMC2', gene.column = 2, unique.features = TRUE)
PBMC2<-Read10X(data.dir = 'PBMC4', gene.column = 2, unique.features = TRUE)
PBMC3_origin<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
PBMC3d<-Read10X(data.dir = 'PBMC3_7000', gene.column = 2, unique.features = TRUE)
### DSB data
DSB1_matrix<-read.csv("DSB_adt+vdj_PBMC1_filtered.csv",row.names = 1)
DSB2_matrix<-read.csv("DSB_adt+vdj_PBMC2_filtered.csv",row.names = 1)
DSB3_matrix<-read.csv("DSB_adt+vdj_PBMC3_filtered.csv",row.names = 1)
gene<-gene.list$Up_Neutrophils[gene.list$Up_Neutrophils%in%rownames(PBMC1$`Gene Expression`)]
colnames(PBMC1$`Gene Expression`)<-paste(colnames(PBMC1$`Gene Expression`),"1",sep = "_")
colnames(PBMC2$`Gene Expression`)<-paste(colnames(PBMC2$`Gene Expression`),"2",sep = "_")
colnames(PBMC3d$`Gene Expression`)<-paste(colnames(PBMC3d$`Gene Expression`),"3",sep = "_")
PBMC1_5000<-CreateSeuratObject(counts = PBMC1$`Gene Expression`)
PBMC1_5000[['Protein']] = CreateAssayObject(counts = normalized_matrix)
PBMC2_3500<-CreateSeuratObject(counts = PBMC2$`Gene Expression`)
PBMC2_3500[['Protein']] = CreateAssayObject(counts = DSB2_matrix)
PBMC3<-CreateSeuratObject(counts = PBMC3d$`Gene Expression`)
PBMC3[['Protein']] = CreateAssayObject(counts = DSB3_matrix)

### Raw ADT
PBMC1_5000<-CreateSeuratObject(counts = PBMC1$`Gene Expression`)
PBMC1_5000[['Protein']] = CreateAssayObject(counts = PBMC1$`Antibody Capture`)
PBMC2_3500<-CreateSeuratObject(counts = PBMC2$`Gene Expression`)
PBMC2_3500[['Protein']] = CreateAssayObject(counts = PBMC2$`Antibody Capture`)
PBMC3<-CreateSeuratObject(counts = PBMC3d$`Gene Expression`)
PBMC3[['Protein']] = CreateAssayObject(counts = PBMC3d$`Antibody Capture`)
###
PBMC1_5000<-AddMetaData(PBMC1_5000,'PBMC1_5000',col.name = 'batch')
PBMC2_3500<-AddMetaData(PBMC2_3500,'PBMC2_3500',col.name = 'batch')
PBMC3<-AddMetaData(PBMC3,'PBMC3',col.name = 'batch')

### PBMC1_5000
###############

PBMC1_5000 <- NormalizeData(PBMC1_5000, assay = "Protein", normalization.method = "CLR")
PBMC1_5000 <- ScaleData(PBMC1_5000, assay = "Protein")

PBMC1_5000<-SCTransform(PBMC1_5000,verbose = T)
PBMC1_5000<-RunPCA(PBMC1_5000,verbose = F)
PBMC1_5000<-RunUMAP(PBMC1_5000,dims = 1:30)
PBMC1_5000<-FindNeighbors(PBMC1_5000,dims = 1:30,verbose = F)
PBMC1_5000<-FindClusters(PBMC1_5000,verbose = F,resolution = 0.8)
### Gate neutrophil
plot<-FeatureScatter(PBMC1_5000,feature1 = "CD16-TotalSeqC",feature2 = "HLA-DR-TotalSeqC",slot = "counts")
neutrophil<-CellSelector(plot,PBMC1_5000,ident = "neutrophil")
plot2<-DimPlot(PBMC1_5000,label = T)
neutrophil<-CellSelector(plot2,PBMC1_5000,ident = "neutrophil")
plot3<-FeaturePlot(PBMC1_5000,label = T,features = "MPO")
neutrophil<-CellSelector(plot3,PBMC1_5000,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(PBMC1_5000@meta.data[which(Idents(PBMC1_5000)=="8"),])
Idents(PBMC1_5000)<-"Non_neutro"
PBMC1_5000<-SetIdent(PBMC1_5000,cells=c(neutrophil),value = "neutrophil")

DimPlot(PBMC1_5000,cells.highlight = neutrophil,label = T)

#FeaturePlot(PBMC1_5000,features = gene.list$Up_Neutrophils)
#DoHeatmap(PBMC1_5000,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(PBMC1_5000, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### new gating

plot<-FeatureScatter(PBMC1_5000,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
base1<-CellSelector(plot,PBMC1_5000,ident = "base")
potential1<-rownames(base1@meta.data[which(Idents(base1)=="base"),])

plot2<-FeaturePlot(base1,features = "MPO",cells = potential1)
neutrophil1<-CellSelector(plot2,PBMC1_5000,ident = "neutrophil")
cell_name1<-rownames(neutrophil1@meta.data[which(Idents(neutrophil1)=="neutrophil"),])
neutrophil1<-SetIdent(base1,cells=c(cell_name1),value = "neutrophil")
FeatureScatter(neutrophil1,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data",cells = potential1)





### DEG
mc1.markers11<-FindMarkers(PBMC1_5000,ident.1 = 11,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers11),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers19<-FindMarkers(PBMC1_5000,ident.1 = 19,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers19),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers4<-FindMarkers(PBMC1_5000,ident.1 = 4,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers4),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers2<-FindMarkers(PBMC1_5000,ident.1 = 2,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers2),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers18<-FindMarkers(PBMC1_5000,ident.1 = 18,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers18),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers20<-FindMarkers(PBMC1_5000,ident.1 = 20,only.pos = T)
DotPlot(PBMC1_5000,features = head(rownames(mc1.markers20),15),col.min = 0,dot.scale = 4)+RotatedAxis()

###sta
Data_neutrophil<-data.frame(umi=PBMC1_5000$nCount_RNA[cell_name1], row.names=names(PBMC1_5000$nCount_RNA[cell_name1]),group="neutrohil")
mu <- ddply(Data_neutrophil, "group", summarise, grp.mean=mean(umi))
p<-ggplot(Data_neutrophil,aes(x=umi,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")


### Sta
umisum<-colSums(as.matrix(PBMC1$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/PBMC1_5000$nCount_RNA[neutrophil]
control<-PBMC1$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/PBMC1_5000$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")

mono1<-names(PBMC1_5000$seurat_clusters[PBMC1_5000$seurat_clusters==2])
mono1<-c(mono1,names(PBMC1_5000$seurat_clusters[PBMC1_5000$seurat_clusters==21]))
mono1<-c(mono1,names(PBMC1_5000$seurat_clusters[PBMC1_5000$seurat_clusters==15]))

mono1.2<-sample5_5000$seurat_clusters[sample5_5000$seurat_clusters==14]

### PBMC2_3500
#################
PBMC2_3500 <- NormalizeData(PBMC2_3500, assay = "Protein", normalization.method = "CLR")
PBMC2_3500 <- ScaleData(PBMC2_3500, assay = "Protein")


PBMC2_3500<-SubsetData(PBMC2_3500,subset.name = "nFeature_RNA",low.threshold = 40)
PBMC2_3500<-SubsetData(PBMC2_3500,subset.name = "nCount_RNA",low.threshold = 100)

PBMC2_3500<-SCTransform(PBMC2_3500,verbose = T)
PBMC2_3500<-RunPCA(PBMC2_3500,verbose = F)
PBMC2_3500<-RunUMAP(PBMC2_3500,dims = 1:30)
PBMC2_3500<-FindNeighbors(PBMC2_3500,dims = 1:30,verbose = F)
PBMC2_3500<-FindClusters(PBMC2_3500,verbose = F,resolution = 0.8)

### Gate neutrophil
DefaultAssay(PBMC2_3500)<-"Protein"
plot<-FeatureScatter(PBMC2_3500,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
neutrophil<-CellSelector(plot,PBMC2_3500,ident = "neutrophil")
plot2<-DimPlot(PBMC2_3500,label = T)
neutrophil<-CellSelector(plot2,PBMC2_3500,ident = "neutrophil")
plot3<-FeaturePlot(PBMC2_3500,features = "MPO")
neutrophil<-CellSelector(plot3,PBMC2_3500,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(PBMC2_3500@meta.data[which(Idents(PBMC2_3500)=="8"),])
Idents(PBMC2_3500)<-"Non_neutro"
PBMC2_3500<-SetIdent(PBMC2_3500,cells=c(neutrophil),value = "neutrophil")
DimPlot(PBMC2_3500,cells.highlight = neutrophil,label = T)



### new gating
plot<-FeatureScatter(PBMC2_3500,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
base2<-CellSelector(plot,PBMC2_3500,ident = "base")
potential2<-rownames(base2@meta.data[which(Idents(base2)=="base"),])

plot2<-FeaturePlot(base2,features = "MPO",cells = potential2)
neutrophil2<-CellSelector(plot2,PBMC2_3500,ident = "neutrophil")
cell_name2<-rownames(neutrophil2@meta.data[which(Idents(neutrophil2)=="neutrophil"),])
neutrophil2<-SetIdent(base2,cells=c(cell_name2),value = "neutrophil")
FeatureScatter(neutrophil2,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data",cells = potential2)

#FeaturePlot(PBMC2_3500,features = gene.list$Up_Neutrophils)
#DoHeatmap(PBMC2_3500,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(PBMC2_3500, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### DEG
mc2.markers1<-FindMarkers(PBMC2_3500,ident.1 = 1,only.pos = T)
DotPlot(PBMC2_3500,features = head(rownames(mc2.markers1),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers2<-FindMarkers(PBMC2_3500,ident.1 = 2,only.pos = T)
DotPlot(PBMC2_3500,features = head(rownames(mc2.markers2),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers8<-FindMarkers(PBMC2_3500,ident.1 = 8,only.pos = T)
DotPlot(PBMC2_3500,features = head(rownames(mc2.markers8),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers9<-FindMarkers(PBMC2_3500,ident.1 = 9,only.pos = T)
DotPlot(PBMC2_3500,features = head(rownames(mc2.markers9),15),col.min = 0,dot.scale = 4)+RotatedAxis()

### Statistic
umisum<-colSums(as.matrix(PBMC2$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/PBMC2_3500$nCount_RNA[neutrophil]
control<-PBMC2$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/PBMC2_3500$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")

mono2<-names(PBMC2_3500$seurat_clusters[PBMC2_3500$seurat_clusters==6])
mono2<-c(mono2,names(PBMC2_3500$seurat_clusters[PBMC2_3500$seurat_clusters==13]))
mono2<-c(mono2,names(PBMC2_3500$seurat_clusters[PBMC2_3500$seurat_clusters==15]))

mono2.2<-sample6_3500$seurat_clusters[sample6_3500$seurat_clusters==14]

### PBMC3
##########
mono3.name<-read.csv("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/PBMC3_mono_mono.csv")
colnames(mono3.name) <- gsub("_3","",colnames(mono3.name))  
mono3.name<-colnames(mono3.name)

diff.cell3<-setdiff(colnames(PBMC3d$`Antibody Capture`),colnames(PBMC3_origin$`Antibody Capture`))
PBMC3 <- NormalizeData(PBMC3, assay = "Protein", normalization.method = "CLR")
PBMC3 <- ScaleData(PBMC3, assay = "Protein")

PBMC3<-SCTransform(PBMC3,verbose = T)
PBMC3<-RunPCA(PBMC3,verbose = F)
PBMC3<-RunUMAP(PBMC3,dims = 1:30)
PBMC3<-FindNeighbors(PBMC3,dims = 1:30,verbose = F)
PBMC3<-FindClusters(PBMC3,verbose = F)

FeaturePlot(PBMC3,features = "MPO")
### Gate neutrophil
plot<-FeatureScatter(PBMC3,feature1 = "CD16-TotalSeqC",feature2 = "HLA-DR-TotalSeqC",slot = "counts")
neutrophil<-CellSelector(plot,PBMC3,ident = "neutrophil")
plot2<-DimPlot(PBMC3,label = T)
neutrophil<-CellSelector(plot2,PBMC3,ident = "neutrophil")
plot3<-FeaturePlot(PBMC3,features = "MPO")
neutrophil<-CellSelector(plot3,PBMC3,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(PBMC3@meta.data[which(Idents(PBMC3)=="18"),])
Idents(PBMC3)<-"Non_neutro"
PBMC3<-SetIdent(PBMC3,cells=c(neutrophil2,neutrophil3),value = "neutrophil")
neutrophil<-c(neutrophil2,neutrophil3)
DimPlot(PBMC3,cells.highlight = neutrophil,label=T)

neutrophil.gene<-FindMarkers(PBMC3,ident.1 = c("neutrophil2","neutrophil3"),only.pos = T)
neutrophil.gene<-neutrophil.gene[which(neutrophil.gene["p_val_adj"]<0.05),]
diff_gene<-FindMarkers(PBMC3,ident.1 = "neutrophil2",ident.2 = "neutrophil3")
diff_gene<-diff_gene[which(diff_gene["p_val_adj"]<0.05),]
diff_gene_name<-head(rownames(diff_gene),10)
neutrophil_gene_name<-head(rownames(neutrophil.gene),18)

DotPlot(PBMC3,features = "PRSS57",col.min = 0,dot.scale = 4)+RotatedAxis()
#FeaturePlot(PBMC3,features = gene.list$Up_Neutrophils)
#DoHeatmap(PBMC3,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(PBMC3, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### 
mc3.markers8<-FindMarkers(PBMC3,ident.1 = 8,only.pos = T)  ### RBC
DotPlot(PBMC3,features = head(rownames(mc3.markers8),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers17<-FindMarkers(PBMC3,ident.1 = 17,only.pos = T) ### Eosinophil
DotPlot(PBMC3,features = head(rownames(mc3.markers17),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers18<-FindMarkers(PBMC3,ident.1 = 18,only.pos = T)
DotPlot(PBMC3,features = head(rownames(mc3.markers18),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers23<-FindMarkers(PBMC3,ident.1 = 23,only.pos = T)
DotPlot(PBMC3,features = head(rownames(mc3.markers23),15),col.min = 0,dot.scale = 4)+RotatedAxis()
###
umisum<-colSums(as.matrix(PBMC3d$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/PBMC3$nCount_RNA[neutrophil]
control<-PBMC3d$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/PBMC3$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")

mono3<-names(PBMC3$seurat_clusters[PBMC3$seurat_clusters==9])
mono3<-c(mono3,names(PBMC3$seurat_clusters[PBMC3$seurat_clusters==19]))
mono3<-c(mono3,names(PBMC3$seurat_clusters[PBMC3$seurat_clusters==22]))
mono3.2<-names(sample7$seurat_clusters[sample7$seurat_clusters==14])

extra.mono3<-setdiff(mono3.2,mono3)
