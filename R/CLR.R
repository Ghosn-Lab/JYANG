library(compositions)
setwd('F:/Ghosn_lab/supperseq/')
D<-read.csv('concat.csv',row.names = 1)

#D<-t(D)
#D<-clr(D)
#D<-t(D)
#write.csv(D,file='CLR_log_407.csv',row.names = TRUE)
x<-t(D)
n<-t(D)
for (i in 1:nrow(x))
  n[i,]<-log1p(x = x[i,] / (exp(x = sum(log1p(x = x[i,][x[i,] > 0]), na.rm = TRUE) / length(x = x[i,]))))

write.csv(n,file='CLR_concat_pos.csv',row.names = TRUE)

#D<-t(D)
D<-clr(D)
#D<-t(D)
write.csv(D,file='CLR_concat_origin.csv',row.names = TRUE)

#exp(x = sum(log1p(x = x[i,][x[i,] > 0]), na.rm = TRUE) / length(x = x[i,]))