setwd("F:/Ghosn_lab/ADT/ADT/")
library(liger)
library(cowplot)
library(ggplot2)
BL_data = "F:/Ghosn_lab/gene_expression/Bloodsample1/BEI"
BM_data = "F:/Ghosn_lab/gene_expression/BMsample1/BMRNA/expression"
FL_data = "F:/Ghosn_lab/gene_expression/FLsample1/FL"

BL<-read10X(BL_data,sample.names = c('BL'))
BM<-read10X(BM_data,sample.names = c('BM'))
FL<-read10X(FL_data,sample.names = c('FL'))


### Data preprocessing (normalization, feature selection, scaling)
liger10X <- createLiger(list(BL = BL$`Gene Expression`, BM = BM, FL = FL$`Gene Expression`))


liger10X <- normalize(liger10X)
liger10X <- selectGenes(liger10X, var.thresh = 0.1)
liger10X <- scaleNotCenter(liger10X)

### Factorization
"k.suggest <- suggestK(liger10X, num.cores = 5, gen.new = T, return.results = T, plot.log2 = F,nrep = 5)"
liger10X <- optimizeALS(liger10X, k=22, thresh = 5e-5, nrep = 3)

### Ploting before Alignment
liger10X <- runTSNE(liger10X, use.raw = T)
p1 <- plotByDatasetAndCluster(liger10X, return.plots = T)
# Plot by dataset
print(p1[[1]])

### Quantile alignment
liger10X <- quantileAlignSNF(liger10X, resolution = 0.4, small.clust.thresh = 20)



### Visulization
liger10X <- runTSNE(liger10X)

# factor plots into pdf file
pdf("plot_factors.pdf")
plotFactors(liger10X)
dev.off()

p_a <- plotByDatasetAndCluster(liger10X, return.plots = T) 
# Modify plot output slightly
p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) + 
  guides(col=guide_legend(title = '', override.aes = list(size = 4)))
print(p_a[[1]])


### plot t-SNE plot colored by clusters 

print(p_a[[2]])

markers <- getFactorMarkers(liger10X, dataset1='sample6', dataset2='Blood2', num.genes = 10)
head(markers$sample6[markers$sample6$factor_num == 1, ])

word_clouds <- plotWordClouds(liger10X, num.genes = 10, do.spec.plot = F, return.plots = T)
print(word_clouds[[9]])

### plot t-SNE plot colored by gene
p_g2 <- plotGene(liger10X, 'CD4', return.plots = T,axis.labels = c('tsne1','tsne2'))
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

