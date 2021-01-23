if(!exists("Processed.count.data")){source("function_data_setting_CS1.R")}

media <<- "YPD" 
genotype <<- c("WT(ho)", "stp2")
Create.BASICS <- TRUE

Geno.Count.List = list()
for(i in c(1:2)){
  Geno.Count.List[[i]] <-  
    Processed.count.data %>% 
    filter(Genotype_Group == genotype[i]) %>% 
    filter(Condition == media) %>% 
    select(-KANMX, -NATMX, -Genotype, -Genotype_Group, -Replicate, -Condition, -tenXBarcode) %>%
    t()
}
names(Geno.Count.List) = genotype

Multi.Count = 
  Processed.count.data %>% 
  filter(Condition == media) %>% 
  filter(Genotype_Group == genotype[1] | Genotype_Group == genotype[2]) %>% 
  select(-KANMX, -NATMX, -Condition, -tenXBarcode)
Multi.Count$Genotype_Group <- factor(Multi.Count$Genotype_Group, levels=genotype)
Multi.Count = arrange(Multi.Count, Genotype_Group)

Cell.id <-
  rownames(Multi.Count)
Condition.id <- 
  Multi.Count$Genotype_Group %>% 
  sapply(function(geno){if(geno=="WT(ho)"){return(1)}else{return(2)}}) %>% 
  as.matrix()

rownames(Condition.id) <- rownames(Multi.Count)
Multi.Count.Matrix = Multi.Count %>% select(-Genotype, -Genotype_Group, -Replicate) %>% t()
sce = SingleCellExperiment(assays=list(counts=Multi.Count.Matrix), 
                     colData=data.frame(condition=Condition.id))
sce$Genotype_Replicate = Multi.Count$Genotype
sce$Genotype_Group = Multi.Count$Genotype_Group
sce$Replicate = Multi.Count$Replicate

#Normalization
Norm.sce.zh1 <- preprocess(sce, zero.thresh=1, scran_norm=TRUE)
Norm.sce.zh.99 <- preprocess(sce, zero.thresh=0.99, scran_norm=TRUE)


SNNGraph <- buildSNNGraph(normcounts(Norm.sce.zh1))
clusteredDF = 
  data.frame(
    rowname = Cell.id , 
    Genotype = Condition.id, 
    Cluster = as.factor(igraph::cluster_louvain(SNNGraph)$membership)
  )

if(0){source("UMAP_Mapper_CS3.R")}

Non.sum.zero <- function(df){
  sum = colSums(select(Geno1FilteredData, -Replicate))
  return(names(sum[sum != 0]))
}

if(Create.BASICS){
  #Non.Sum.zero.genes = rownames(Norm.sce.zh.99)
  for(index in c(1:6)){
    ClusterIndex = index - 1
    dir = dirlist_0.99[index]
    print(dir)
    if(ClusterIndex==0){All.cluster = TRUE}else{All.cluster = FALSE}
    if(All.cluster){
      ExtractedCells1 = clusteredDF %>%  filter(Genotype == 1) %>% nth(1)
      ExtractedCells2 = clusteredDF %>%  filter(Genotype == 2) %>% nth(1)
    }else{
      ExtractedCells1 = clusteredDF %>% filter(Genotype == 1) %>% 
        filter(Cluster == ClusterIndex) %>% nth(1)
      ExtractedCells2 = clusteredDF %>% filter(Genotype == 2) %>% 
        filter(Cluster == ClusterIndex) %>% nth(1)
    }
    Geno1FilteredData = Multi.Count[ExtractedCells1, c(Non.sum.zero(Multi.Count[ExtractedCells1,]), 'Replicate')]
    Geno2FilteredData = Multi.Count[ExtractedCells2, c(Non.sum.zero(Multi.Count[ExtractedCells2,]), 'Replicate')]
    print(dim(Geno1FilteredData))
    print(dim(Geno2FilteredData))
    #Create BASiCS Chain Object
    Geno1ChainName = paste(genotype[1], "Regression", sep="")
    Geno2ChainName = paste(genotype[2], "Regression", sep="")
    CreateChainObject(Geno1FilteredData, Geno1ChainName, dir)
    CreateChainObject(Geno2FilteredData, Geno2ChainName, dir)
  }
}


