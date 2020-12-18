#座標軸となる変数は列
library(tagcloud)
data = t(normcounts(NormSceobj)[1:976,])

genotype = colData(NormSceobj)[1]$condition[1:976]
arrnames = c(WTMainData$Replicate) %>% unique()#, MutantMainData$Replicate) %>% unique()

library(RColorBrewer)
n <- length(arrnames)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, n) %>% setNames(arrnames)

rep = c(WTMainData$Replicate) %>% 
  sapply(function(x){
    return(colors[x])
  }) %>% as.array()
  
pcaobj <- prcomp(data, scale = FALSE)

PC1 <- pcaobj$x[, 1]
PC2 <- pcaobj$x[, 2]
PC3 <- pcaobj$x[, 3]
PC4 <- pcaobj$x[, 4]

plot(PC1, PC2)
plot(PC2, PC3)

plot(PC1, PC2, pch=genotype, col=col)
plot(PC2, PC3, pch=genotype, col=col)
plot(PC3, PC4, pch=genotype, col=col)


library(dbscan)
data("DS3")
cl <- sNNclust(data.frame(PC1,PC2), k = 20, eps = 7, minPts = 16)
plot(data.frame(PC1,PC2), col = cl$cluster + 1L, cex = .5, pch=genotype)

cl <- sNNclust(data.frame(PC2,PC3), k = 20, eps = 7, minPts = 16)
plot(data.frame(PC2,PC3), col = cl$cluster + 1L, cex = .5, pch=genotype)

cl <- sNNclust(data.frame(PC3,PC4), k = 20, eps = 7, minPts = 16)
plot(data.frame(PC3,PC4), col = cl$cluster + 1L, cex = .5, pch=genotype)


