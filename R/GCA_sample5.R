library(Seurat)
library(qusage)
library(ggplot2)
library(plyr)

### filtered matrix
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/filtered/")

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA4_CD45neg_21461_5GEX_G5",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"1",sep = "_")
E7<-CreateSeuratObject(counts = filtered)
rm(filtered)

filtered<-Read10X("GCA5_CD45pos_21478_5GEX_A6",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"2",sep = "_")
A6<-CreateSeuratObject(counts = filtered)
rm(filtered)

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA5_CD45neg_21478_5GEX_E7",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"3",sep = "_")
D7<-CreateSeuratObject(counts = filtered)
rm(filtered)

filtered<-Read10X("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Gut Cell Atlas/Matrix/unfiltered/GCA6_CD45neg_20006_5GEX_D6",gene.column = 2, unique.features = TRUE)
colnames(filtered)<-paste(colnames(filtered),"4",sep = "_")
H5<-CreateSeuratObject(counts = filtered)
rm(filtered)

E7<-SubsetData(E7,subset.name = "nFeature_RNA",low.threshold = 40)
E7<-SubsetData(E7,subset.name = "nCount_RNA",low.threshold = 100)
D7<-SubsetData(D7,subset.name = "nFeature_RNA",low.threshold = 40)
D7<-SubsetData(D7,subset.name = "nCount_RNA",low.threshold = 100)
H5<-SubsetData(H5,subset.name = "nFeature_RNA",low.threshold = 40)
H5<-SubsetData(H5,subset.name = "nCount_RNA",low.threshold = 100)

E7<-AddMetaData(E7,'Collagenase_CD45-',col.name = 'batch')
A6<-AddMetaData(A6,'Collagenase_CD45+',col.name = 'batch')
H5<-AddMetaData(H5,'Protease',col.name = 'batch')
D7<-AddMetaData(D7,'Collagenase_total',col.name = 'batch')


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



E7<-SCTransform(E7,verbose = F)
A6<-SCTransform(A6,verbose = F)
D7<-SCTransform(D7,verbose = F)
ev<-E7@assays$SCT@var.features
av<-A6@assays$SCT@var.features
dv<-D7@assays$SCT@var.features
Cov<-unique(c(ev,av,dv))



CD45n<-FindMarkers(Collagenase,only.pos = T,ident.1 = 20)### Gut cells
CD45n9<-FindMarkers(Collagenase,only.pos = T,ident.1 = 9) ### Dying cells


CD45p0<-FindMarkers(Collagenase,only.pos = T,ident.1 = 0) ### B cell
CD45p7<-FindMarkers(Collagenase,only.pos = T,ident.1 = 7) ### B cell
CD45p4<-FindMarkers(Collagenase,only.pos = T,ident.1 = 4) ### NK cells 
CD45p24<-FindMarkers(Collagenase,only.pos = T,ident.1 = 24) ### Proliferating cells
CD45p5<-FindMarkers(Collagenase,only.pos = T,ident.1 = 5) ### NKT
CD45p2<-FindMarkers(Collagenase,only.pos = T,ident.1 = 2) ### Treg
CD45p3<-FindMarkers(Collagenase,only.pos = T,ident.1 = 3) ### T cells
CD45p1<-FindMarkers(Collagenase,only.pos = T,ident.1 = 1) ###
CD45p11<-FindMarkers(Collagenase,only.pos = T,ident.1 = 11) ###

CD45p1_10<-FindMarkers(Collagenase,ident.1 = 1,ident.2 = 10) ###
CD45p25<-FindMarkers(Collagenase,only.pos = T,ident.1 = 25) ### 
CD45p8<-FindMarkers(Collagenase,only.pos = T,ident.1 = 8) ### 
CD45p23<-FindMarkers(Collagenase,only.pos = T,ident.1 = 23) ### 
CD45p14<-FindMarkers(Collagenase,only.pos = T,ident.1 = 14) ### 
CD45p21<-FindMarkers(Collagenase,only.pos = T,ident.1 = 21) ### 



### Dotplot
DotPlot(Collagenase,features = head(rownames(CD45p0),n=20),col.min = 0,dot.scale = 4)+RotatedAxis()


c0<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==0])
CD8nT<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==1])
Tcells<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==2])
Treg<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==3])
NK_cells<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==4])
NKT<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==5])
c6<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==6])
c7<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==7])
Monocytes<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==8])
c9<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==9])
CD8pT<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==10])
c11<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==11])
c12<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==12])
c13<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==13])
c14<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==14])
c15<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==15])
c16<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==16])
c17<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==17])
c18<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==18])
c19<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==19])
Gutcells<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==20])
Mastcells<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==21])
c22<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==22])
c23<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==23])
Proliferating<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==24])
c25<-names(Collagenase$seurat_clusters[Collagenase$seurat_clusters==25])

B_cells<-c(c0,c6,c7,c9,c11,c12,c13,c14,c15,c16,c17,c18,c19,c22)
sB_cells<-c(c0,c6,c11,c12,c13,c14,c15,c16,c17,c18,c19,c22)

Collagenase<-SetIdent(Collagenase,cells=c(B_cells),value = "B")
CD45pw<-FindMarkers(Collagenase,only.pos = T,ident.1 = "B") ### B cells
Collagenase<-FindClusters(Collagenase,verbose = F,resolution = 0.8)

######
###### Protease
H5<-SCTransform(H5,verbose = F)
H5<-RunPCA(H5,verbose = F)
H5<-RunUMAP(H5,verbose = F,dims = 1:30)
H5<-FindNeighbors(H5,dims = 1:30,verbose = F)
H5<-FindClusters(H5,verbose = F,resolution = 0.8)

Pro0<-FindMarkers(H5,only.pos = T,ident.1 = 0)   ### ileum 
Pro13<-FindMarkers(H5,only.pos = T,ident.1 = 13) ### ileum 
Pro16<-FindMarkers(H5,only.pos = T,ident.1 = 16)
Pro<-FindMarkers(H5,ident.1 = 0,ident.2 = 13)
Pro14<-FindMarkers(H5,only.pos = T,ident.1 = 14)
Pro11<-FindMarkers(H5,only.pos = T,ident.1 = 11)
Pro19<-FindMarkers(H5,only.pos = T,ident.1 = 19)
Pro4<-FindMarkers(H5,only.pos = T,ident.1 = 4)
Pro8<-FindMarkers(H5,only.pos = T,ident.1 = 8)
Pro1<-FindMarkers(H5,only.pos = T,ident.1 = 1)
Pro3<-FindMarkers(H5,only.pos = T,ident.1 = 3)
Pro6<-FindMarkers(H5,only.pos = T,ident.1 = 6)
Pro7<-FindMarkers(H5,only.pos = T,ident.1 = 7)
Pro18<-FindMarkers(H5,only.pos = T,ident.1 = 18)
Pro20<-FindMarkers(H5,only.pos = T,ident.1 = 20)
Pro5<-FindMarkers(H5,only.pos = T,ident.1 = 5)
Pro17<-FindMarkers(H5,only.pos = T,ident.1 = 17)
Pro15<-FindMarkers(H5,only.pos = T,ident.1 = 15)
Pro2<-FindMarkers(H5,only.pos = T,ident.1 = 2)
Pro10<-FindMarkers(H5,only.pos = T,ident.1 = 10)
Pro12<-FindMarkers(H5,only.pos = T,ident.1 = 12) ### Plasma

### Dotplot
DotPlot(H5,features = head(rownames(Pro6),n=20),col.min = 0,dot.scale = 4)+RotatedAxis()

p5<-names(H5$seurat_clusters[H5$seurat_clusters==5])
p15<-names(H5$seurat_clusters[H5$seurat_clusters==15])
p17<-names(H5$seurat_clusters[H5$seurat_clusters==17])
p12<-names(H5$seurat_clusters[H5$seurat_clusters==12])
p9<-names(H5$seurat_clusters[H5$seurat_clusters==9])
p10<-names(H5$seurat_clusters[H5$seurat_clusters==10])
P_B_cells<-c(p9,p10,p12)
P_Mono<-c(p5,p15)
P_T<-c(p17)

DimPlot(Merged,cells.highlight = list(P_B_cells,B_cells),label = T,cols.highlight = c("darkolivegreen3","mediumpurple1"),split.by = "orig.ident")
DimPlot(Merged,cells.highlight = list(P_Mono,Monocytes),label = T,cols.highlight = c("darkolivegreen3","mediumpurple1"),split.by = "orig.ident")
DimPlot(Merged,cells.highlight = list(P_T,c(Tcells,c1,c2,c3,c4,c5,c10)),label = T,cols.highlight = c("darkolivegreen3","mediumpurple1"),split.by = "orig.ident")


mB<-names(Merged$seurat_clusters[Merged$seurat_clusters==21])

Merged[["percent.mt"]] <- PercentageFeatureSet(Merged, pattern = "^MT-")
DimPlot(Merged,group.by = "batch",cols = c("green","yellow","red","mediumpurple1"))


plot2<-DimPlot(protease,label = T)
non<-CellSelector(plot2,protease,ident = "Non-immune")
a<-rownames(non@meta.data[which(Idents(non)=="Non-immune"),])

