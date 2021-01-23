library(tidyverse)
library(scran)
library(scDD)
library(BASiCS)
setwd("/Users/itoutouma/Lab_Analysis/scRNAseq")
dir = "./"
dirlist_0.99 = list.files(dir)[grep("N0", list.files(dir))] %>% sort()

build_sceObj_from_readcount <- function(Multi.Count, geno1, geno2){
  Condition <- 
    Multi.Count$Genotype_Group %>% 
    sapply(function(geno){if(geno=="WT(ho)"){return(1)}else{return(2)}}) %>% 
    as.matrix()
  rownames(Condition) <- rownames(Multi.Count)
  
  Multi.Count.Matrix = Multi.Count %>% select(-Genotype_Group, -Replicate) %>% t()
  SingleCellExperiment(assays=list(counts=Multi.Count.Matrix), 
                         colData=data.frame(condition)) %>% return()
}


detect_DD_from_normsceObj <- function(norm_sce){
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  scDatExSim <- scDD(norm_sce, prior_param=prior_param, testZeroes=TRUE)
  return(scDatExSim)
}
CreateChainObject <- function(FilteredData, RunName, dir){
  Batch <- c(1:length(unique(FilteredData$Replicate))) %>% 
    set_names(unique(FilteredData$Replicate))
  CountRep <- FilteredData %>% select(Replicate) %>% 
    nth(1) %>% sapply(function(x){return(Batch[x])})
  
  FilteredCount = FilteredData %>% select(-Replicate) %>% as.matrix() %>% t()
  SCE <- SingleCellExperiment(assays = list(counts = FilteredCount), 
                              colData = data.frame(BatchInfo = CountRep))
  Regression <- BASiCS_MCMC(Data = SCE, N = 20000, 
                            Thin = 20, Burn = 10000, 
                            Regression = TRUE , 
                            PrintProgress = TRUE,
                            WithSpikes = FALSE,
                            StoreChains = TRUE,
                            StoreDir = dir,
                            RunName = RunName)
}


#dirlist_0.99 <- c("YPD_WTsSTP2_zh0.99", "YPD_WTsSTP2_zh0.99_c1", "YPD_WTsSTP2_zh0.99_c2",
#                 "YPD_WTsSTP2_zh0.99_c3", "YPD_WTsSTP2_zh0.99_c4", "YPD_WTsSTP2_zh0.99_c5")

#Data loading

##processed data (Doublets and low-count cells have already been removed )
##BASiCS use pre-processed data. For set condition identical, scDD also use this data
Processed.count.data <- read.table("../elife-51254-code2.tsv", header=T, sep="\t")
#Jackson et al. Elife. 2020 

Yeastract.interact = read.csv("/Users/itoutouma/Lab_Analysis/GRN/RegulationTwoColumnTable_both.tsv", sep=';', header = FALSE)

BIOGRID.interact = read.csv("../DFGs/BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.59.tab.txt", sep="\t")

Yeast.gene.name = read.csv("/Users/itoutouma/Lab_Analysis/GRN/orftogene.txt", sep="\t", header=T) %>% 
  mutate(SGD = ORF)
Yeast.gene.name$ORF = apply(Yeast.gene.name, MARGIN=1, 
                        function(row){ if(length(grep("-", row[2]))){ return(gsub("-", "\\.", row[2]))}else{ return(row[2])}})


