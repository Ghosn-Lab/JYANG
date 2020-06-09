library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/")
options(future.globals.maxSize = 4000 * 1024^2)

gene.list<-read.gmt("Geneset-Library-Granulocytes.gmt")

extra_cells<-read.csv("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/extra_cells_.csv")
sample5<-Read10X(data.dir = 'PBMC1', gene.column = 2, unique.features = TRUE)
sample6<-Read10X(data.dir = 'PBMC2', gene.column = 2, unique.features = TRUE)
sample7d<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/")
PBMC1_force<-Read10X(data.dir = 'PBMC1', gene.column = 2, unique.features = TRUE)
PBMC2_force<-Read10X(data.dir = 'PBMC2', gene.column = 2, unique.features = TRUE)
PBMC3_force<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
colnames(PBMC1_force$`Gene Expression`)<-paste(colnames(PBMC1_force$`Gene Expression`),"1",sep = "_")
colnames(PBMC2_force$`Gene Expression`)<-paste(colnames(PBMC2_force$`Gene Expression`),"2",sep = "_")
colnames(PBMC3_force$`Gene Expression`)<-paste(colnames(PBMC3_force$`Gene Expression`),"3",sep = "_")
### DSB data
DSB1_matrix<-read.csv("DSB_adt+vdj_sample5_filtered.csv",row.names = 1)
DSB2_matrix<-read.csv("DSB_adt+vdj_sample6_filtered.csv",row.names = 1)
DSB3_matrix<-read.csv("DSB_adt+vdj_sample7_filtered.csv",row.names = 1)
gene<-gene.list$Up_Neutrophils[gene.list$Up_Neutrophils%in%rownames(sample5$`Gene Expression`)]
colnames(sample5$`Gene Expression`)<-paste(colnames(sample5$`Gene Expression`),"1",sep = "_")
colnames(sample6$`Gene Expression`)<-paste(colnames(sample6$`Gene Expression`),"2",sep = "_")
colnames(sample7d$`Gene Expression`)<-paste(colnames(sample7d$`Gene Expression`),"3",sep = "_")
colnames(sample5$`Antibody Capture`)<-paste(colnames(sample5$`Antibody Capture`),"1",sep = "_")
colnames(sample6$`Antibody Capture`)<-paste(colnames(sample6$`Antibody Capture`),"2",sep = "_")
colnames(sample7d$`Antibody Capture`)<-paste(colnames(sample7d$`Antibody Capture`),"3",sep = "_")
sample5_5000<-CreateSeuratObject(counts = sample5$`Gene Expression`)
sample5_5000[['Protein']] = CreateAssayObject(counts = DSB1_matrix)
sample6_3500<-CreateSeuratObject(counts = sample6$`Gene Expression`)
sample6_3500[['Protein']] = CreateAssayObject(counts = DSB2_matrix)
sample7<-CreateSeuratObject(counts = sample7d$`Gene Expression`)
sample7[['Protein']] = CreateAssayObject(counts = DSB3_matrix)

### Raw ADT
sample5_5000<-CreateSeuratObject(counts = sample5$`Gene Expression`)
sample5_5000[['Protein']] = CreateAssayObject(counts = sample5$`Antibody Capture`)
sample6_3500<-CreateSeuratObject(counts = sample6$`Gene Expression`)
sample6_3500[['Protein']] = CreateAssayObject(counts = sample6$`Antibody Capture`)
sample7<-CreateSeuratObject(counts = sample7d$`Gene Expression`)
sample7[['Protein']] = CreateAssayObject(counts = sample7d$`Antibody Capture`)
###
sample5_5000<-AddMetaData(sample5_5000,'sample5_5000',col.name = 'batch')
sample6_3500<-AddMetaData(sample6_3500,'sample6_3500',col.name = 'batch')
sample7<-AddMetaData(sample7,'sample7',col.name = 'batch')


sample5_5000<-SubsetData(sample5_5000,subset.name = "nFeature_RNA",low.threshold = 40)
sample5_5000<-SubsetData(sample5_5000,subset.name = "nCount_RNA",low.threshold = 100)
sample6_3500<-SubsetData(sample6_3500,subset.name = "nFeature_RNA",low.threshold = 40)
sample6_3500<-SubsetData(sample6_3500,subset.name = "nCount_RNA",low.threshold = 100)
sample7<-SubsetData(sample7,subset.name = "nFeature_RNA",low.threshold = 40)
sample7<-SubsetData(sample7,subset.name = "nCount_RNA",low.threshold = 100)


Merged<-merge(sample5_5000, y = c(sample6_3500,sample7), project = "Merged", merge.data = TRUE)
Merged<-SCTransform(Merged,verbose = T)
Merged<-RunPCA(Merged,verbose = F)
Merged<-RunUMAP(Merged,dims = 1:30)
Merged<-FindNeighbors(Merged,dims = 1:30,verbose = F)
Merged<-FindClusters(Merged,verbose = F,resolution = 0.8)
### sample5_5000
###############
sample5_5000 <- NormalizeData(sample5_5000, assay = "Protein", normalization.method = "CLR")
sample5_5000 <- ScaleData(sample5_5000, assay = "Protein")

sample5_5000<-SCTransform(sample5_5000,verbose = T)
sample5_5000<-RunPCA(sample5_5000,verbose = F)
sample5_5000<-RunUMAP(sample5_5000,dims = 1:30)
sample5_5000<-FindNeighbors(sample5_5000,dims = 1:30,verbose = F)
sample5_5000<-FindClusters(sample5_5000,verbose = F,resolution = 0.8)


diff.cell1<-setdiff(names(sample5_5000$orig.ident),colnames(PBMC1_force$`Antibody Capture`))

FeaturePlot(sample5_5000,features = "MPO",label = T)
DimPlot(sample5_5000,cells.highlight = diff.cell1,label = T)
### Gate neutrophil
plot<-FeatureScatter(sample5_5000,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
neutrophil<-CellSelector(plot,sample5_5000,ident = "neutrophil")
plot2<-DimPlot(sample5_5000,label = T)
neutrophil<-CellSelector(plot2,sample5_5000,ident = "neutrophil")
plot3<-FeaturePlot(sample5_5000,label = T,features = "MPO")
neutrophil<-CellSelector(plot3,sample5_5000,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(sample5_5000@meta.data[which(Idents(sample5_5000)=="8"),])
Idents(sample5_5000)<-"Non_neutro"
sample5_5000<-SetIdent(sample5_5000,cells=c(neutrophil),value = "neutrophil")

DimPlot(sample5_5000,cells.highlight = neutrophil,label = T)
FeatureScatter(sample5_5000,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
#FeaturePlot(sample5_5000,features = gene.list$Up_Neutrophils)
#DoHeatmap(sample5_5000,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(sample5_5000, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### new gating

plot<-FeatureScatter(sample5_5000,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
base1<-CellSelector(plot,sample5_5000,ident = "base")
potential1<-rownames(base1@meta.data[which(Idents(base1)=="base"),])

plot2<-FeaturePlot(base1,features = "MPO",cells = potential1)
neutrophil1<-CellSelector(plot2,sample5_5000,ident = "neutrophil")
cell_name1<-rownames(neutrophil1@meta.data[which(Idents(neutrophil1)=="neutrophil"),])
neutrophil1<-SetIdent(base1,cells=c(cell_name1),value = "neutrophil")
FeatureScatter(neutrophil1,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data",cells = potential1)





### DEG
mc1.markers11<-FindMarkers(sample5_5000,ident.1 = 11,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers11),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers19<-FindMarkers(sample5_5000,ident.1 = 19,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers19),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers4<-FindMarkers(sample5_5000,ident.1 = 4,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers4),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers2<-FindMarkers(sample5_5000,ident.1 = 2,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers2),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers18<-FindMarkers(sample5_5000,ident.1 = 18,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers18),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc1.markers20<-FindMarkers(sample5_5000,ident.1 = 20,only.pos = T)
DotPlot(sample5_5000,features = head(rownames(mc1.markers20),15),col.min = 0,dot.scale = 4)+RotatedAxis()

###sta
Data_neutrophil<-data.frame(umi=sample5_5000$nCount_RNA[cell_name1], row.names=names(sample5_5000$nCount_RNA[cell_name1]),group="neutrohil")
mu <- ddply(Data_neutrophil, "group", summarise, grp.mean=mean(umi))
p<-ggplot(Data_neutrophil,aes(x=umi,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")


### Sta
umisum<-colSums(as.matrix(sample5$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/sample5_5000$nCount_RNA[neutrophil]
control<-sample5$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/sample5_5000$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")

### sample6_3500
#################
sample6_3500 <- NormalizeData(sample6_3500, assay = "Protein", normalization.method = "CLR")
sample6_3500 <- ScaleData(sample6_3500, assay = "Protein")

sample6_3500<-SCTransform(sample6_3500,verbose = T)
sample6_3500<-RunPCA(sample6_3500,verbose = F)
sample6_3500<-RunUMAP(sample6_3500,dims = 1:30)
sample6_3500<-FindNeighbors(sample6_3500,dims = 1:30,verbose = F)
sample6_3500<-FindClusters(sample6_3500,verbose = F,resolution = 0.8)


diff.cell2<-setdiff(names(sample6_3500$orig.ident),colnames(PBMC2_force$`Antibody Capture`))

FeaturePlot(sample6_3500,features = "MPO",label = T)
DimPlot(sample6_3500,cells.highlight = diff.cell2,label = T)
### Gate neutrophil
DefaultAssay(sample6_3500)<-"Protein"
plot<-FeatureScatter(sample6_3500,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
neutrophil<-CellSelector(plot,sample6_3500,ident = "neutrophil")
plot2<-DimPlot(sample6_3500,label = T)
neutrophil<-CellSelector(plot2,sample6_3500,ident = "neutrophil")
plot3<-FeaturePlot(sample6_3500,label = T,features = "MPO")
neutrophil<-CellSelector(plot3,sample6_3500,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(sample6_3500@meta.data[which(Idents(sample6_3500)=="8"),])
Idents(sample6_3500)<-"Non_neutro"
sample6_3500<-SetIdent(sample6_3500,cells=c(neutrophil),value = "neutrophil")
DimPlot(sample6_3500,cells.highlight = neutrophil,label = T)



### new gating
plot<-FeatureScatter(sample6_3500,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data")
base2<-CellSelector(plot,sample6_3500,ident = "base")
potential2<-rownames(base2@meta.data[which(Idents(base2)=="base"),])

plot2<-FeaturePlot(base2,features = "MPO",cells = potential2)
neutrophil2<-CellSelector(plot2,sample6_3500,ident = "neutrophil")
plot3<-FeaturePlot(sample6_3500,features = "MPO")
neutrophil2<-CellSelector(plot2,sample6_3500,ident = "neutrophil")
cell_name2<-rownames(neutrophil2@meta.data[which(Idents(neutrophil2)=="neutrophil"),])
neutrophil2<-SetIdent(base2,cells=c(cell_name2),value = "neutrophil")
FeatureScatter(neutrophil2,feature1 = "CD16-TotalSeqC",feature2 = "CD45-TotalSeqC",slot = "scale.data",cells = potential2)

#FeaturePlot(sample6_3500,features = gene.list$Up_Neutrophils)
#DoHeatmap(sample6_3500,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(sample6_3500, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### DEG
mc2.markers1<-FindMarkers(sample6_3500,ident.1 = 1,only.pos = T)
DotPlot(sample6_3500,features = head(rownames(mc2.markers1),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers2<-FindMarkers(sample6_3500,ident.1 = 2,only.pos = T)
DotPlot(sample6_3500,features = head(rownames(mc2.markers2),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers8<-FindMarkers(sample6_3500,ident.1 = 8,only.pos = T)
DotPlot(sample6_3500,features = head(rownames(mc2.markers8),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc2.markers9<-FindMarkers(sample6_3500,ident.1 = 9,only.pos = T)
DotPlot(sample6_3500,features = head(rownames(mc2.markers9),15),col.min = 0,dot.scale = 4)+RotatedAxis()

### Statistic
umisum<-colSums(as.matrix(sample6$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/sample6_3500$nCount_RNA[neutrophil]
control<-sample6$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/sample6_3500$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")
### sample7
##########
diff.cell3<-setdiff(colnames(sample7d$`Antibody Capture`),colnames(sample7_origin$`Antibody Capture`))
sample7 <- NormalizeData(sample7, assay = "Protein", normalization.method = "CLR")
sample7 <- ScaleData(sample7, assay = "Protein")

sample7<-SCTransform(sample7,verbose = T)
sample7<-RunPCA(sample7,verbose = F)
sample7<-RunUMAP(sample7,dims = 1:30)
sample7<-FindNeighbors(sample7,dims = 1:30,verbose = F)
sample7<-FindClusters(sample7,verbose = F)

diff.cell3<-setdiff(names(sample7$orig.ident),colnames(PBMC3d$`Antibody Capture`))
FeaturePlot(sample7,features = "MPO",label = T)
DimPlot(sample7,cells.highlight = diff.cell3,label=T)
### Gate neutrophil
plot<-FeatureScatter(sample7,feature1 = "CD16-TotalSeqC",feature2 = "HLA-DR-TotalSeqC",slot = "scale.data",cells = neutrophil)
neutrophil<-CellSelector(plot,sample7,ident = "neutrophil")
plot2<-DimPlot(sample7,label = T)
neutrophil<-CellSelector(plot2,sample7,ident = "neutrophil")
plot3<-FeaturePlot(sample7,features = "MPO")
neutrophil<-CellSelector(plot3,sample7,ident = "neutrophil")
neutrophil<-rownames(neutrophil@meta.data[which(Idents(neutrophil)=="neutrophil"),])
neutrophil<-rownames(sample7@meta.data[which(Idents(sample7)=="18"),])
Idents(sample7)<-"Non_neutro"
sample7<-SetIdent(sample7,cells=c(neutrophil),value = "neutrophil")
neutrophil<-c(neutrophil2,neutrophil3)
DimPlot(sample7,cells.highlight = neutrophil,label=T)

neutrophil.gene<-FindMarkers(sample7,ident.1 = c("neutrophil2","neutrophil3"),only.pos = T)
neutrophil.gene<-neutrophil.gene[which(neutrophil.gene["p_val_adj"]<0.05),]
diff_gene<-FindMarkers(sample7,ident.1 = "neutrophil2",ident.2 = "neutrophil3")
diff_gene<-diff_gene[which(diff_gene["p_val_adj"]<0.05),]
diff_gene_name<-head(rownames(diff_gene),10)
neutrophil_gene_name<-head(rownames(neutrophil.gene),18)

DotPlot(sample7,features = "PRSS57",col.min = 0,dot.scale = 4)+RotatedAxis()
#FeaturePlot(sample7,features = gene.list$Up_Neutrophils)
#DoHeatmap(sample7,cells = neutrophil,features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")
#DoHeatmap(subset(sample7, downsample = 100),features = gene.list$Up_Neutrophils,slot = "data",assay = "SCT")

### 
mc3.markers8<-FindMarkers(sample7,ident.1 = 8,only.pos = T)  ### RBC
DotPlot(sample7,features = head(rownames(mc3.markers8),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers2<-FindMarkers(sample7,ident.1 = 2,only.pos = T) ### Eosinophil
DotPlot(sample7,features = head(rownames(mc3.markers2),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers18<-FindMarkers(sample7,ident.1 = 18,only.pos = T)
DotPlot(sample7,features = head(rownames(mc3.markers18),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers22<-FindMarkers(sample7,ident.1 = 22,only.pos = T)
DotPlot(sample7,features = head(rownames(mc3.markers22),15),col.min = 0,dot.scale = 4)+RotatedAxis()
mc3.markers23<-FindMarkers(sample7,ident.1 = 23,only.pos = T)
DotPlot(sample7,features = head(rownames(mc3.markers23),15),col.min = 0,dot.scale = 4)+RotatedAxis()
###
umisum<-colSums(as.matrix(sample7d$`Gene Expression`[gene,neutrophil]))
percent_neutro<-umisum*100/sample7$nCount_RNA[neutrophil]
control<-sample7d$`Gene Expression`[gene,]

control<-colSums(as.matrix(control))[colSums(as.matrix(control))>0]
control<-control[setdiff(names(control),names(percent_neutro))]
percent_control<-control*100/sample7$nCount_RNA[names(control)]
Data_control<-data.frame(percent=percent_control, row.names=names(percent_control),group="control")
Data_neutrophil<-data.frame(percent=percent_neutro, row.names=names(percent_neutro),group="neutrohil")
Data<-rbind(Data_control,Data_neutrophil)

mu <- ddply(Data, "group", summarise, grp.mean=mean(percent))
p<-ggplot(Data,aes(x=percent,fill=group))+geom_density(alpha=0.4)
p+geom_vline(data=mu, aes(xintercept=grp.mean, color=group),linetype="dashed")


