library(scDD)
library(dplyr)
###TURORIAL of SCDD
data(scDatExSim)
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
scDatExSim <- scDD(scDatExSim, prior_param=prior_param, testZeroes=FALSE)

results(scDatExSim)
#Clusters.combined is number of mode in c1 c2 combind distribution
#Clusters.c1 is number of mode in c1
#Clusters.c1 is number of mode in c2

PARTITION.C1 <- results(scDatExSim, type="Zhat.c1")
PARTITION.C2 <- results(scDatExSim, type="Zhat.c2")
#PARTITION.C1 indicates that which sample is belongs to the modes

i = 15
sideViolin(normcounts(scDatExSim)[i,], scDatExSim$condition,
           MAP = list(PARTITION.C1[i,], PARTITION.C2[i,]),
           title.gene=results(scDatExSim)$DDcategory[i])
#Reference: https://bioconductor.org/packages/release/bioc/manuals/scDD/man/scDD.pdf

