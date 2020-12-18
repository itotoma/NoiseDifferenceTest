library(ggplot2)
library(reshape2)
#Reference: https://github.com/andrewheiss/Attitudes-in-the-Arab-World/blob/master/figure12.R

###mean & noise version
#HeatmapDF from MainAnalysis
DfMode = 0
if(DfMode){
  data = HeatmapDf %>% column_to_rownames("Gene") %>% select(noise, mean) %>% as.matrix()
  df  <- melt(data) %>% setNames(c("Gene", "Group", "Pvalue")) %>% 
    left_join(HeatmapDf, by=c("Gene")) %>% 
    mutate(stars = ifelse((r1 != "NoDiff") & (Group == "mean"), "*", "")) %>% 
    mutate(stars = ifelse((r2 != "NoDiff") & (Group == "noise"), "*", "")) %>% 
    select(Gene, Group, Pvalue, stars)
  #Clustering and rearrange
}


###noise by each media version
ListMode = 1
NaOmit = 1
if(ListMode){
  RedundautTargets = data.frame(V2 = c("AGP1", "BAP2", "BAP3", "DIP5", "GNP1", "MUP1", "TAT1", "TAT2", "PTR2"))
  
  if(1){
    #Mean Analysis
    data = MeanList %>% lapply(function(df){
      df %>% right_join(yeast_gene_name, by=c("GeneName"="ORF")) %>% 
        mutate(star = ifelse(ResultDiffMean != "NoDiff", "*", "")) %>% 
        #filter(is.element(Name, RedundautTargets$V2)) %>% 
        select(Name, MeanLog2FC, star) %>% 
        setNames(c("Gene", "Value", "star")) %>% 
        return()})
  }
  if(0){
    #Noise Analysis
    data = ResDispList %>% lapply(function(df){
      df %>% right_join(yeast_gene_name, by=c("GeneName"="ORF")) %>% 
        mutate(star = ifelse(ResultDiffResDisp != "NoDiff", "*", "")) %>% 
        filter(is.element(Name, RedundautTargets$V2)) %>% 
        select(Name, ResDispDistance, star) %>% 
        setNames(c("Gene", "Value", "star")) %>% 
        return()})
  }
  
  PreClueterdData = 
      data %>% do.call(cbind.data.frame,.) %>% 
      column_to_rownames("Glutamine_WTvsSTP2.Gene")
  ClueterdData = 
      PreClueterdData[,-grep("(Gene|star)",names(PreClueterdData))] %>% setNames(dirlist) %>% 
      as.matrix()
  df = data %>% do.call(rbind.data.frame,.) %>% rownames_to_column() %>% 
    setNames(c("Group", "Gene", "Value", "star")) %>% 
    mutate(Group = gsub("_WTvsSTP2.[0-9]+", "", Group))
  if(NaOmit){
    ClueterdData = na.omit(ClueterdData)
    df = na.omit(df)
  }else{
    ClueterdData[is.na(ClueterdData)]=0
  }
}
dev.new()
clr <- heatmap(ClueterdData, scale = "none")
dev.off()
gene.idx  <- rownames(ClueterdData)[clr$rowInd]
group.idx <- colnames(ClueterdData)[clr$colInd]
df$Gene  <- factor(df$Gene, levels = gene.idx)
df$Group <- factor(df$Group, levels = rev(c("MinimalEtOH", "MinimalGlucose", "Glutamine", "Urea", "YPDRapa", "YPDDiauxic", "Proline", "AmmoniumSulfate", "YPEtOH", "YPD", "CStarve")))
#df$Group <- factor(df$Group, levels = group.idx)

#generate heatmap
ghm <- ggplot(df, aes(x = Gene, y = Group, fill = Value)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label=star),color="black", size=5) + 
  theme(plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   #axis.text.x = element_blank()
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
                  )
ghm <- ghm + xlab("Genes") + ylab("") + 
  scale_fill_gradient2(low="#00A0E9", mid="white",high="#E60012") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  ggtitle("Mean Difference Probability")
ghm


