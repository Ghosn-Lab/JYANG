library(dsb)
library(Seurat)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq")
options(future.globals.maxSize = 4000 * 1024^4)

### unfiltered matrix
pbmc.htos<-Read10X(data.dir = 'C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/unfiltered/BM2_4000', gene.column = 2, unique.features = TRUE)
### filtered matrix
pbmc.umis<-Read10X(data.dir = 'BM2_4000', gene.column = 2, unique.features = TRUE)
hto<-CreateSeuratObject(counts = pbmc.htos$`Gene Expression`)

### Remove empty ADT rows
protein<- as.matrix(pbmc.htos$`Antibody Capture`)
protein<-protein[which(rowSums(protein)>0),]
protein<-protein[which(rownames(protein)!="CD38_TotalSeqC"),]


### Select negtive cells (background)
hto[['Protein']] = CreateAssayObject(counts = protein)
neg<-SubsetData(hto,subset.name = "nFeature_RNA",high.threshold = 40)
neg<-GetAssayData(hto, assay = "Protein", slot = 'counts') %>% as.matrix()

### get ADT data for real cells
real<-as.matrix(pbmc.umis$`Antibody Capture`)
### Remove empty ADT rows from real cells
real<-real[which(rowSums(real)>0),]
real<-real[which(rownames(real)!="CD38_TotalSeqC"),]

### DSB normalization
normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = real, empty_drop_matrix = neg)

### Add index to barcodes
colnames(normalized_matrix) <- paste(colnames(normalized_matrix), "3", sep = "_")
write.csv(normalized_matrix,"DSB_PBMC3_7000.csv")