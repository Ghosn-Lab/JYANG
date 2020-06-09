library(ggplot2)

setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/")

D<-read.csv('DSB_BM2+3_vdj_CP.csv',row.names = 1)
D<-data.frame(t(D))

e<-ggplot(data= D, aes(x=CD19_TotalSeqC,y=umis+1))

e+geom_point()+scale_y_continuous(trans='log10')



kde2d(D$CD19, D$CD20, n = 25, lims = c(D$CD19, D$CD20))