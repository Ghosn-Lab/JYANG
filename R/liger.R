setwd("C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/")
library(liger)
library(cowplot)
library(ggplot2)
sample5_data = "C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM1"
sample6_data = "C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM2"
sample7_data = "C:/Users/yjk12/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM3"
#Bloodsample2_data = "F:/Ghosn_lab/gene_expression/Bloodsample2"
#Bloodsample1_data = "F:/Ghosn_lab/gene_expression/Bloodsample1/BEI"
#FLsample1_data = "F:/Ghosn_lab/gene_expression/FLsample1/FL"

sample5<-read10X(sample5_data,sample.names = c('BM1'))
sample6<-read10X(sample6_data,sample.names = c('BM2'))
sample7<-read10X(sample7_data,sample.names = c('BM3'))
#Bloodsample2<-read10X(Bloodsample2_data,sample.names = c('sample2'))
#Bloodsample1<-read10X(Bloodsample1_data,sample.names = c('bloodsample1'))
#FLsample1<-read10X(FLsample1_data,sample.names = c('FLsample1'))

### Data preprocessing (normalization, feature selection, scaling)
liger10X <- createLiger(list(BM1 = sample5$`Gene Expression`, BM2 = sample6$`Gene Expression`, BM3 = sample7$`Gene Expression`))
#liger10X <- createLiger(list(sample1 = sample5, sample2 = Bloodsample2 ))

#liger10X <- createLiger(list(Blood1 = Bloodsample1$`Gene Expression`, FL = FLsample1$`Gene Expression` ))

liger10X <- normalize(liger10X)
liger10X <- selectGenes(liger10X, var.thresh = 0.1)
liger10X <- scaleNotCenter(liger10X)

### Factorization
"k.suggest <- suggestK(liger10X, num.cores = 5, gen.new = T, return.results = T, plot.log2 = F,nrep = 5)"
liger10X <- optimizeALS(liger10X, k=22, thresh = 5e-5, nrep = 3)

### Ploting before Alignment
liger10X <- runUMAP(liger10X,use.raw = T)
p1 <- plotByDatasetAndCluster(liger10X, return.plots = T)
# Plot by dataset
p1[[1]]<-p1[[1]]+facet_grid(. ~ Dataset)
p1[[2]]<-p1[[2]]+facet_grid(. ~ Dataset)
print(p1[[1]])
print(p1[[2]])
### Quantile alignment
liger10X <- quantile_norm(liger10X)



### Visulization
liger10X <- runUMAP(liger10X)

# factor plots into pdf file
pdf("plot_factors.pdf")
plotFactors(liger10X)
dev.off()

p_a <- plotByDatasetAndCluster(liger10X, return.plots = T) 
# Modify plot output slightly
p_a[[1]] <- p_a[[1]] + theme_classic() + 
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))+facet_grid(. ~ Dataset)
print(p_a[[1]])


### plot t-SNE plot colored by clusters 
p_a[[2]]<-p_a[[2]]+facet_grid(. ~ Dataset)
print(p_a[[2]])

markers <- getFactorMarkers(liger10X, dataset1='sample6', dataset2='Blood2', num.genes = 10)
head(markers$sample6[markers$sample6$factor_num == 1, ])

word_clouds <- plotWordClouds(liger10X, num.genes = 10, do.spec.plot = F, return.plots = T)
print(word_clouds[[9]])

### plot t-SNE plot colored by gene
p_g2 <- plotGene(liger10X, 'IgD', return.plots = T,axis.labels = c('tsne1','tsne2'))
plot_grid(plotlist = p_g2)

p_3<-plotGeneViolin(liger10X,'CD4',methylation.indices = NULL,
                    by.dataset = T, return.plots = F)

plotGeneLoadings(liger10X, dataset1 = 'Blood1', dataset2 = 'FL',
                 num.genes.show = 12, num.genes = 30, mark.top.genes = T,
                 factor.share.thresh = 10, log.fc.thresh = 1, umi.thresh = 30,
                 frac.thresh = 0, pval.thresh = 0.05, do.spec.plot = T,
                 max.val = 0.1, pt.size = 0.1, option = "plasma",
                 zero.color = "#F5F5F5", return.plots = F)

plotClusterProportions(liger10X, return.plot = F)



### Test quantile
test<- createLiger(list(sample1 = Y, sample2 = Z ))
test <- normalize(test)
test <- scaleNotCenter(test)
test <- quantile_norm(test)

