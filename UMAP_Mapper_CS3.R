library(scran)
library(umap)

#Eliminate cell state variety
if(!exists("clusteredDF")){source("BASiCS_DataPrepareCS2.R")}

#G1-phase specific marker: PIR1 (YKL164C)
PIR1Level = normcounts(Norm.sce.zh1)['YKL164C',] %>% 
  as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "PIR1Level"))

#G1-phase daughter-cell specific marker DSE2,: DSE2 (YHR143W)
DSE2Level = normcounts(Norm.sce.zh1)['YHR143W', ] %>%
  as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "DSE2Level"))

#GS-phase specific marker histone 2B: HTB (HTB1 / YDR224C, HTB2 / YBL002W)
HTBLevel = (normcounts(Norm.sce.zh1)['YDR224C',] + normcounts(Norm.sce.zh1)['YBL002W',]) %>%
  as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "HTBLevel"))

#Decompose dimension by UMAP
umap_score <- umap(t(normcounts(Norm.sce.zh1)))
UMAP_result = data.frame(PC1=umap_score$layout[,1], PC2=umap_score$layout[,2]) %>% rownames_to_column()
#WARNING sNNclust parameters strongly affects the clustering result
#Cluster Mapping
g <- UMAP_result %>% left_join(clusteredDF) %>% 
  ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color=Cluster, alpha=0.8)) + 
  ggtitle("SNN Clustering") + scale_alpha(guide='none')
print(g)
ggsave(file.path(dir, "SNNclustering.png"), g, width = 5.6, height = 4.3)

#Batch Mapping
g <- Processed.count.data %>% select(Replicate) %>% 
  rownames_to_column() %>% right_join(UMAP_result) %>% 
  ggplot(aes(x=PC1, y=PC2)) + geom_point(aes(color=Replicate, alpha=0.8)) + 
  ggtitle("BatchMapping") + scale_alpha(guide='none')
ggsave(file.path(dir, "BatchMapping.png"), g, width = 5.6, height = 4.3)


g <- PIR1Level %>% left_join(UMAP_result) %>% 
  ggplot(aes(x=PC1, y=PC2))+geom_point(aes(color=PIR1Level, alpha=0.8)) + 
  scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012") +
  ggtitle("PIR1 Mapping")
print(g)
ggsave(file.path(dir, "PIR1Mapping.png"), g, width = 5.6, height = 4.3)

g <- DSE2Level %>% left_join(UMAP_result) %>% 
  ggplot(aes(x=PC1, y=PC2))+geom_point(aes(color=DSE2Level, alpha=0.8)) + 
  scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012") +
  ggtitle("DSE2 Mapping")
print(g)
ggsave(file.path(dir, "DSE2Mapping.png"), g, width = 5.6, height = 4.3)

g <- HTBLevel %>% left_join(UMAP_result) %>% 
  ggplot(aes(x=PC1, y=PC2))+geom_point(aes(color=HTBLevel, alpha=0.8)) + 
  scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012") +
  ggtitle("HTB Mapping")
print(g)
ggsave(file.path(dir, "HTBMapping.png"), g, width = 5.6, height = 4.3)

#8 HTA2 well localized
#10 HHF1 (chromatin assembly)
#11 HHT1 (chromatin assembly)

g <- clusteredDF %>% left_join(UMAP_result) %>% ggplot(aes(x=PC1, y=PC2))+geom_point(aes(color=as.factor(Genotype), alpha=0.8)) + 
  scale_alpha(guide='none') + scale_color_hue(name="Genotype", labels=c('WT', 'ΔSTP2')) + 
  ggtitle("Genotype Mapping")
print(g)
ggsave(file.path(dir, "GenotypeMapping.png"), g, width = 5.6, height = 4.3)
  
#For test plot of arbitrary gene 
see.well.categorized.gene = FALSE
if(see.well.categorized.gene){
  source("AnovaOnCluster_CS2.5.R")
  print(arrange(Cluster.comp.2, Anova.pval)[1:5, ])
}

normcounts(Norm.sce.zh1)['YNL030W', ] %>%
  as.data.frame() %>% rownames_to_column() %>% set_names(c("rowname", "ExLevel")) %>% 
  left_join(UMAP_result) %>% 
  ggplot(aes(x=PC1, y=PC2))+geom_point(aes(color=ExLevel, alpha=0.8)) + 
  scale_alpha(guide='none') + scale_color_gradient(low = "black", high = "#E60012")



