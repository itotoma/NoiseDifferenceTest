###Turoeial of scran Normalization

#####
"Reference" 
'https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#8_Detecting_correlated_genes'
#####

library(scran)
#import scuttle
for(f in list.files("/Users/itoutouma/Lab_Analysis/scuttle/R/")){
  source(paste("/Users/itoutouma/Lab_Analysis/scuttle/R/",f,sep=""))
}

data(scDatExList)
condition <- c(rep(1, ncol(scDatExList$C1)), rep(2, ncol(scDatExList$C2)))
names(condition) <- c(colnames(scDatExList$C1), colnames(scDatExList$C2))
sce <- SingleCellExperiment(
  assays=list(counts= cbind(scDatExList$C1, scDatExList$C2)), 
  colData=data.frame(condition))
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
summary(sizeFactors(sce))
sce <- logNormCounts(sce)




dec <- modelGeneVar(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
# Get the top 10% of genes.
top.havgs <- getTopHVGs(dec, prop=0.1)

# Get the top 2000 genes.
top.hvgs2 <- getTopHVGs(dec, n=2000)

# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(dec, var.threshold=0)

# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)


