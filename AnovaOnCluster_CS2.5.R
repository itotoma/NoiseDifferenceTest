C1.cells = clusteredDF %>% filter(Cluster == 1) %>% nth(1)
C2.cells = clusteredDF %>% filter(Cluster == 2) %>% nth(1)
C3.cells = clusteredDF %>% filter(Cluster == 3) %>% nth(1)
C4.cells = clusteredDF %>% filter(Cluster == 4) %>% nth(1)
C5.cells = clusteredDF %>% filter(Cluster == 5) %>% nth(1)

Anova.df <- normcounts(Norm.sce.zh.99) %>% t() %>% 
  as.data.frame() %>% rownames_to_column() %>% left_join(clusteredDF, by="rowname") %>% 
  column_to_rownames("rowname")

non.sum.zero.genes = rownames(normcounts(Norm.sce.zh.99))

clust.by.gene.meandf <- data.frame(
  C1 = Anova.df[C1.cells,] %>% select(-Cluster, -Genotype) %>% colMeans(),
  C2 = Anova.df[C2.cells,] %>% select(-Cluster, -Genotype) %>% colMeans(),
  C3 = Anova.df[C3.cells,] %>% select(-Cluster, -Genotype) %>% colMeans(),
  C4 = Anova.df[C4.cells,] %>% select(-Cluster, -Genotype) %>% colMeans(),
  C5 = Anova.df[C5.cells,] %>% select(-Cluster, -Genotype) %>% colMeans()
)

Anova.pval = sapply(non.sum.zero.genes, function(gene){
  anova.res = anova(aov(Anova.df[,gene]~Anova.df$Cluster))
  return(anova.res$`Pr(>F)`[1])
}) %>% as.data.frame() %>% set_names(c("Anova.pval"))

Cluster.comp = cbind(clust.by.gene.meandf, Anova.pval)
Cluster.comp = Cluster.comp %>% 
  mutate(maxval = rowMaxs(as.matrix(select(Cluster.comp, -Anova.pval))))
Cluster.comp.2 = Cluster.comp %>% 
  mutate(Maxclust = ifelse(C1 == maxval, "C1", "")) %>% 
  mutate(Maxclust = ifelse(C2 == maxval, "C2", Maxclust)) %>% 
  mutate(Maxclust = ifelse(C3 == maxval, "C3", Maxclust)) %>% 
  mutate(Maxclust = ifelse(C4 == maxval, "C4", Maxclust)) %>% 
  mutate(Maxclust = ifelse(C5 == maxval, "C5", Maxclust)) %>% 
  select(-maxval)
rownames(Cluster.comp.2) = non.sum.zero.genes

Cluster.comp.2 %>% arrange(Anova.pval)


