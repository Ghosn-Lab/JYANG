library(iTALK)
library(Seurat)
setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/Devon/Adding sequence/")

sample1<-Read10X(data.dir = 'LungTotalnew_01', gene.column = 2, unique.features = TRUE)
A2<-CreateSeuratObject(counts = sample1$`Gene Expression`)

A2<-SCTransform(A2)
