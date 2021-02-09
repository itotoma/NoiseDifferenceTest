library(DESeq2)
if(!exists("clusteredDF")){source("BASiCS_DataPrepareCS2.R")}
Create.DSEeq2.result = TRUE
if(Create.DSEeq2.result){
  ClusterMeanDiff.DESeq2 <- 
    lapply(c(0:5), function(ci){
      print(ci)
      if(ci != 0){
        filter.group = dplyr::filter(clusteredDF, Cluster==ci)$rowname
      }else{
        filter.group = Cell.id
      }
      filter.count = Multi.Count.Matrix[,filter.group]
      filter.condition.id = data.frame(con=factor(Condition.id[filter.group,]))
      dds <- DESeqDataSetFromMatrix(countData = filter.count, colData = filter.condition.id, design = ~ con)
      dds <- DESeq(dds)
      res <- results(dds)
      print(summary(res))
      res %>% as.data.frame() %>% rownames_to_column() %>% 
      return()
  })
} %>% set_names(c("All", "C1", "C2", "C3", "C4", "C5"))
ClusterMeanDiff.DESeq2 %>% enframe() %>% unnest() %>%  
  setNames(c("Cluster", "ORF", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) %>% 
  write.csv("MeanDiff_DESeq_YPD.csv")

MeanDiff.DESeq2.df = read.csv("MeanDiff_DESeq_YPD.csv") %>% 
  left_join(Yeast.gene.name) %>% 
  mutate(star = ifelse(padj < 0.01, "*", "")) %>% 
  select(Name, Cluster, log2FoldChange, star) %>% 
  add_row(Name="Joint", Cluster=c("All", "C1", "C2", "C3", "C4", "C5"), log2FoldChange = 0, star = "")

TF.Target = Yeastract.interact %>% dplyr::filter(V1 == "RTG1") %>% dplyr::select(V2)
gene.idx_joint = c(TF.Target$V2, "Joint", setdiff(unique(MeanDiff.DESeq2.df$Name),c(TF.Target$V2, "Joint")))
MeanDiff.DESeq2.df$Name  <- factor(MeanDiff.DESeq2.df$Name, levels = gene.idx_joint)
#df$Group <- factor(df$Group, levels = rev(c("MinimalEtOH", "MinimalGlucose", "Glutamine", "Urea", "YPDRapa", "YPDDiauxic", "Proline", "AmmoniumSulfate", "YPEtOH", "YPD", "CStarve")))
MeanDiff.DESeq2.df$Cluster <- factor(MeanDiff.DESeq2.df$Cluster, levels = c("C5", "C4", "C3", "C2", "C1", "All"))

#LYS20 #POR1 

ghm <- ggplot(na.omit(MeanDiff.DESeq2.df), aes(x = Name, y = Cluster, fill = log2FoldChange)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label=star),color="black", size=5) + 
  geom_segment(x='Joint', xend='Joint', y=0, yend=51.5, size=.9, alpha=0.4) +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_blank()
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.5)
  )
ghm <- ghm + xlab("Genes") + ylab("") + 
  scale_fill_gradient2(low="#00A0E9", mid="white",high="#E60012") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  ggtitle("YPD media | DESeq Mean | WT vs ΔRTG1 ")
ghm
ggsave('DESeqMeanRTG1TargetsSeg.png', ghm, width=12.80, height=2.50)



#####
#DSE2 for differential expression analysis
DDSres= list()
genotypelist = unique(Processed.count.data$Genotype_Group)

DDSresult <- 
  lapply(genotypelist[-8], function(gene){
    print(gene)
    Filtered.data <- 
      Processed.count.data %>% filter(Condition == "YPD") %>% 
      filter((Genotype_Group == "WT(ho)")|(Genotype_Group == gene)) 
    Filtered.condition.id <- 
      Filtered.data %>% select(Genotype_Group) %>% 
      mutate(id = ifelse(Genotype_Group=="WT(ho)", 1, 2))
    Filtered.count <-
      Filtered.data %>% select(-Genotype, -Genotype_Group, -Replicate, -Condition, -tenXBarcode) %>% 
      t() %>% as.data.frame()
    filter.condition.id = data.frame(con=factor(Filtered.condition.id$id))
    
    dds <- DESeqDataSetFromMatrix(countData = Filtered.count, colData = filter.condition.id, design = ~ con)
    dds <- DESeq(dds)
    res <- results(dds)
    return(res)
  })
DDSresultDF <- enframe(DDSresult) %>% unnest()

DDSresultDF <- lapply(DDSresult, function(df){df %>% as.data.frame() %>% 
    na.omit() %>% rownames_to_column() %>% return()})
names(DDSresultDF) <- genotypelist[-8]
DDSresultDF %>% enframe() %>% unnest() %>% write.csv("DESeq2_all_gene_in_YPD.csv")

DESeq.result <- read.csv("DESeq2_all_gene_in_YPD.csv")
DESeq.result.clust <- read.csv("Clust_YPD_WTvsSTP2_Result/MeanDiff_DESeq_YPD.csv")

