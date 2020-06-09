library(ggplot2)

setwd("C:/Users/yjk12/Box/Data/VDJ SuPERR-seq")

fb<-read.csv("./B_cell/B_BM3.csv")
ufb<-read.csv("./B_unfiltered_only/unfiltered_B_BM3.csv")
nb<-read.csv("./Non_B_unfiltered/non_B_BM3.csv")

fb<-data.frame(fb,origin=rep("filtered_B",length(row.names(fb))))
ufb<-data.frame(ufb,origin=rep("unfiltered_B",length(row.names(ufb))))
nb<-data.frame(nb,origin=rep("Non_B",length(row.names(nb))))

Big<-rbind(fb,ufb,nb)

e<-ggplot(data= Big, aes(x=CD19_TotalSeqC,y=CD20_TotalSeqC))
e+geom_point(aes(color=origin, shape=origin))+scale_shape_manual(values=c(19, 2, 1))+xlim(-3,6)+scale_color_manual(values=c("red", "blue", "seagreen3"))
