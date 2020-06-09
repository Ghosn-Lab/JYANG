setwd("C:/Users/yjk12/Box/Data/supperseq")
library(Seurat)
library(ggplot2)

Blood<-Read10X(data.dir = "BM1_V1", gene.column = 1, unique.features = TRUE)

b<-as.data.frame(Blood$`Antibody Capture`)
b<-t(b)
write.csv(b,"BM1_V1/ADT.csv")



test<-read.csv('PBMC3.csv',row.names = 1)
test<-t(test)
a<-rowSums(test)==0
adt<-rownames(test)
adt[a]



setwd("F:/Ghosn_lab/ADT/Fetal Samples")
Blood<-Read10X(data.dir = "./sample5", gene.column = 2, unique.features = TRUE)

b<-as.data.frame(Blood$`Antibody Capture`)