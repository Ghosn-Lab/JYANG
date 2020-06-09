library(dsb)
library(Seurat)
setwd("F:/Ghosn_lab/supperseq")

Raw<-Read10X(data.dir = 'F:/Ghosn_lab/supperseq/unfiltered/PBMC3', gene.column = 2, unique.features = TRUE)
Real<-Read10X(data.dir = 'PBMC3', gene.column = 2, unique.features = TRUE)
filtered<-CreateSeuratObject(counts = Real$`Gene Expression`)
filtered[['Protein']] = CreateAssayObject(counts = Real$`Antibody Capture`)
filtered = subset(filtered, subset = nFeature_RNA > 0)




a<-setdiff(colnames(Raw$`Antibody Capture`),colnames(Real$`Antibody Capture`))

hto_data = as.matrix(Raw$`Antibody Capture`[-nrow(Raw$`Antibody Capture`), a])
rna_data = Raw$`Gene Expression`[, a]
unfiltered<-CreateSeuratObject(counts = rna_data)
unfiltered[['Protein']] = CreateAssayObject(counts = hto_data)

unfiltered <- NormalizeData(unfiltered, assay = "Protein", normalization.method = "CLR")

NegativeObject = subset(unfiltered, subset = nFeature_RNA < 40)
NegativeObject = subset(NegativeObject, subset = nCount_Protein > 0)
NegativeObject = subset(NegativeObject, subset = nFeature_RNA > 0)

neg_adt_matrix = GetAssayData(NegativeObject, assay = "Protein", slot = 'counts') %>% as.matrix()
positive_adt_matrix = GetAssayData(filtered, assay = "Protein", slot = 'counts') %>% as.matrix()
normalized_matrix = DSBNormalizeProtein(cell_protein_matrix = positive_adt_matrix,empty_drop_matrix = neg_adt_matrix, na.rm=TRUE)
