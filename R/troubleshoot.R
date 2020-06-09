library(ggplot2)

B<-read.csv('F:/Ghosn_lab/ADT/VDJ SuPERR-seq/B_cell/B_PBMC3_filtered.csv',row.names = 1)
B_unfiltered_only<-read.csv('F:/Ghosn_lab/ADT/VDJ SuPERR-seq/B_unfiltered_only/B_PBMC3_unfiltered_only.csv',row.names = 1)
Non_B_filtered<-read.csv('F:/Ghosn_lab/ADT/VDJ SuPERR-seq/Non_B_filtered/Non_B_PBMC3_filtered.csv',row.names = 1)
Non_B_unfiltered<-read.csv('F:/Ghosn_lab/ADT/VDJ SuPERR-seq/Non_B_unfiltered/Non_B_PBMC3_unfiltered.csv',row.names = 1)


B<-data.frame(observation = "B",t(B),gear="B")
B_unfiltered_only<-data.frame(observation = "unfiltered_B",t(B_unfiltered_only),gear="unfiltered_B")
Non_B_unfiltered<-data.frame(observation = "Non_B_unfiltered",t(Non_B_unfiltered),gear="Non_B_unfiltered")

'out<-pheatmap(D,cluster_row = FALSE,show_colnames = F,clustering_method = "complete")
out<-pheatmap(C,cluster_row = FALSE,show_colnames = F,clustering_method = "complete")'

'e<-ggplot(data= C, aes(x=CD19,y=CD20))
f<-ggplot(data= D, aes(x=CD19,y=CD20))'

'e+geom_jitter()+xlim(0,4)+ylim(0,1.3)
f+geom_jitter()+xlim(0,4)+ylim(0,1.3)'




'ggplot() +
  geom_jitter(data = C, aes(x=CD19,y=CD20), colour = "red") +
  geom_jitter(data = D, aes(x=CD19,y=CD20), colour = "blue")+xlim(0,4)+ylim(0,1.3)+
  scale_color_discrete(name = "ltitle") '

DATA<-rbind(B,B_unfiltered_only)
DATA<-rbind(DATA,Non_B_unfiltered)


ggplot(data = DATA, aes(x=CD4,y=CD8,colour=factor(observation),shape=factor(gear)))+
  geom_jitter(width = 0.05,height = 0.05)+scale_color_discrete(name = "ltitle") + 
  scale_shape_manual(values = c(0, 8, 2))+xlim(0,4)+ylim(0,1.3)