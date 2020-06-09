sample1<-Read10X(data.dir = 'BM1', gene.column = 2, unique.features = TRUE)
A2<-CreateSeuratObject(counts = sample1$`Gene Expression`)
A2[['Protein']] = CreateAssayObject(counts = sample1$`Antibody Capture`)
sample2<-Read10X(data.dir = 'BM2', gene.column = 2, unique.features = TRUE)
B2<-CreateSeuratObject(counts = sample2$`Gene Expression`)
B2[['Protein']] = CreateAssayObject(counts = sample2$`Antibody Capture`)
sample3<-Read10X(data.dir = 'BM3', gene.column = 2, unique.features = TRUE)
C2<-CreateSeuratObject(counts = sample3$`Gene Expression`)
C2[['Protein']] = CreateAssayObject(counts = sample3$`Antibody Capture`)
C2<-AddMetaData(C2,'BM3',col.name = 'batch')
A2<-AddMetaData(A2,'BM1',col.name = 'batch')
B2<-AddMetaData(B2,'BM2',col.name = 'batch')


Merged <- merge(A2, y = c(B2,C2), project = "Merged", merge.data = TRUE)

Merged <- SCTransform(Merged, verbose = FALSE)
Merged<-RunPCA(Merged,verbose = FALSE)
Merged<-RunUMAP(Merged,dims = 1:30,verbose = FALSE)

Merged<-FindNeighbors(Merged,dims = 1:30,verbose = FALSE)
Merged<-FindClusters(Merged,verbose = FALSE,resolution =0.8)

DimPlot(Merged,label = TRUE,split.by = "batch")+NoLegend()
plots <- DimPlot(Merged, combine = FALSE,group.by = 'batch')
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(plots)