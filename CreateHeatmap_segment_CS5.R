library(ggplot2)
library(ggtext)
library(reshape2)

#Reference: https://github.com/andrewheiss/Attitudes-in-the-Arab-World/blob/master/figure12.R

RedundautTargets = data.frame(V2 = c("AGP1", "BAP2", "BAP3", "DIP5", "GNP1", "MUP1", "TAT1", "TAT2", "PTR2"))
Cluster.id = c("All", "C1", "C2", "C3", "C4", "C5")
STP2.Target = Yeastract.interact %>% filter(V1 == "RTG1") %>% select(V2)
#Mean Analysis
ClusterMeanDiff= enframe(ClusterMeanList) %>% unnest() %>% 
  right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
  mutate(star = ifelse(grepl("\\+", ResultDiffMean), "*", "")) %>%
  select(name, Name, MeanLog2FC, star, ResultDiffMean) %>% 
  setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
  mutate(comp = "Mean")

ClusterNoiseDiff = enframe(ClusterResDispList) %>% unnest() %>% 
  right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
  mutate(star = ifelse(grepl("\\+", ResultDiffResDisp), "*", "")) %>%
  select(name, Name, ResDispDistance, star, ResultDiffResDisp) %>% 
  setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
  mutate(comp = "Noise")

Split.gene.by.interaction = FALSE
Extract.interact.gene = FALSE
if(Split.gene.by.interaction){
  ClusterMeanDiff = ClusterMeanDiff %>% 
    add_row(Group=Cluster.id, Gene=rep('Joint',length(Cluster.id)), Value=0, star="", diff="NoDiff", comp="Mean")
  ClusterNoiseDiff = ClusterNoiseDiff %>% 
    add_row(Group=Cluster.id, Gene=rep('Joint',length(Cluster.id)), Value=0, star="", diff="NoDiff", comp="Noise")
}else if(Extract.interact.gene){
  ClusterMeanDiff = ClusterMeanDiff %>% filter(is.element(Gene, STP2.Target$V2))
  ClusterNoiseDiff = ClusterNoiseDiff %>% filter(is.element(Gene, STP2.Target$V2))
}
Mean.Noise.Bind.df <- rbind(ClusterMeanDiff, ClusterNoiseDiff) %>% na.omit()
gene.idx_joint = c(STP2.Target$V2, "Joint", setdiff(unique(Mean.Noise.Bind.df$Gene),c(STP2.Target$V2, "Joint")))
Mean.Noise.Bind.df$Gene  <- factor(Mean.Noise.Bind.df$Gene, levels = gene.idx_joint)
Mean.Noise.Bind.df$Group <- factor(Mean.Noise.Bind.df$Group, levels = c("C5", "C4", "C3", "C2", "C1", "All"))

#generate heatmap
ghm <- ggplot(Mean.Noise.Bind.df, aes(x = Gene, y = Group, fill = Value)) + 
  geom_tile() + theme_bw() + 
  geom_text(aes(label=star),color="black", size=5) + 
  geom_segment(x='Joint', xend='Joint', y=0, yend=51.5, size=.9, alpha=0.4) +
  facet_wrap(~comp, scale="free_y", nrow=2) +
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
  ggtitle("YPD media | WT vs ΔSTP2 | STP2 downstream genes")
ghm

#ggsave('MeanNoiseCompSTP2Target.png', ghm, width=12.80, height=5.00)


#Cross plot
CrossPlotDf <- ClusterMeanDiff %>% left_join(ClusterNoiseDiff, by=c("Gene", "Group")) %>% 
  na.omit() %>% mutate(signif=ifelse(star.x == "*" | star.y == "*", 1, 0.1)) %>% 
  mutate(STP2TG = ifelse(is.element(Gene, STP2.Target$V2), TRUE, FALSE))

g <- ggplot(CrossPlotDf,aes(x=Value.x, y=Value.y, label=Gene)) + geom_point(aes(colour=STP2TG, alpha=signif)) + 
  xlab("Δmean")+ylab("Δnoise") + guides(alpha=FALSE) + 
  facet_wrap(~Group, scale="free_y", ncol=2)
#ggsave("MeanNoiseCrossFaCluster.png", g, width=12.0, height=7.0)
g <- ggplot(CrossPlotDf,aes(x=Value.x, y=Value.y, label=Gene)) + geom_point(aes(colour=Group, alpha=signif)) + 
  xlab("Δmean")+ylab("Δnoise") + guides(alpha=FALSE) + 
  scale_fill_manual(values = c("brown", "tomato", "yellow4", "mediumaquamarine", "skyblue2", "darkorchid1"))
g
 
#BIOGRID data
BIOGRID.interact.STP2 = BIOGRID.interact %>% 
  filter((OFFICIAL_SYMBOL_A == "STP2") | (OFFICIAL_SYMBOL_B == "STP2")) %>% 
  select(OFFICIAL_SYMBOL_A, OFFICIAL_SYMBOL_B, EXPERIMENTAL_SYSTEM) %>% 
  pivot_longer(c(OFFICIAL_SYMBOL_A, OFFICIAL_SYMBOL_B), names_to = "AB", values_to = "Gene") %>% 
  filter(Gene != "STP2")
CrossPlotDf %>% left_join(BIOGRID.interact.STP2, by="Gene") %>% na.omit()


#Enriched GO term on cross plot
#WT+ W  STP2+ S NoDiff 0
#Menn Noise
Create.cross.gene.csv = FALSE
if(Create.cross.gene.csv){
  file.name.rule = c("W", "S", "N") %>% set_names(c("WT+", "STP2+", "NoDiff"))
  tmpdf <- CrossPlotDf %>% left_join(Yeast.gene.name, by=c("Gene"="Name"))
  for(mean_cond in c("WT+", "STP2+", "NoDiff")){
    for(noise_cond in c("WT+", "STP2+", "NoDiff")){
      filterdf <- tmpdf %>% filter(diff.x==mean_cond) %>% filter(diff.y == noise_cond)
      if(dim(filterdf)[1] != 0){
        if((mean_cond != "NoDiff")||(noise_cond != "NoDiff")){
          file.name = paste(file.name.rule[mean_cond],file.name.rule[noise_cond],".csv",sep="")
          print(file.name)
          filterdf %>% select(SGD) %>% distinct() %>% 
            write.table(file.name, row.names = F, col.names = F)
        }
      }
    }
  }
  rm(tmpdf, filterdf, file.name, file.name.rule)
}

NW_GO = read.csv("NW_GO.txt", sep="\t") %>% mutate(label = paste(GOID, ": ", TERM, "<br>", sep=""))
SN_GO = read.csv("SN_GO.txt", sep="\t")
SS_GO = read.csv("SS_GO.txt", sep="\t") %>% mutate(label = paste(GOID, ": ", TERM, "<br>", sep=""))
WN_GO = read.csv("WN_GO.txt", sep="\t")
SN_GO_lab = SN_GO %>% mutate(num = 1:length(SN_GO$TERM)) %>% 
  mutate(label = ifelse(num%%4==0, paste(GOID, "<br>", sep=""),paste(GOID," ",sep="")))
WN_GO_lab = WN_GO %>% mutate(num = 1:length(WN_GO$TERM)) %>% 
  mutate(label = ifelse(num%%4==0, paste(GOID, "<br>", sep=""),paste(GOID," ",sep="")))

dammuy_size = dim(CrossPlotDf)[1]-8
ggtext_df <- tibble(
  label = c(
    str_sub(paste(NW_GO$label, collapse = ""),end=-5),
    paste(SN_GO_lab$label, collapse = ""),
    str_sub(paste(SS_GO$label, collapse = ""),end=-5),
    paste(WN_GO_lab$label, collapse = ""),
    "No significant GO terms",
    "No sample",
    "No sample",
    "No significant GO terms",
    rep("",dammuy_size)
  ),
  x = c(0, -4, -4, 4, 0, 4, -4, 4, rep(-10,dammuy_size)),
  y = c(5, 0, -4.5, 0, -4.5, -4.5, 5, 5,  rep(-10, dammuy_size)),
  hjust = c(.5, 0, 0, 1, .5, 1, 0, 1, rep(-10, dammuy_size)),
  vjust = c(0, .5, 1, .5, 1, 1, 0, 0, rep(-10, dammuy_size)),
  angle = c(0, 0, 0, 0, 0, 0, 0, 0, rep(-10, dammuy_size)),
  color = c("black", "black", "black", "black", "blue", "blue", "blue", "blue", rep("black", dammuy_size)),
  fill = c("cornsilk", "cornsilk", "cornsilk", "cornsilk", "white", "white",  "white", "white", rep("cornsilk", dammuy_size))
  )
df = cbind(ggtext_df, CrossPlotDf)
g <- 
  ggplot(df) +
  aes(
    x, y, label = label, angle = angle, color = color, fill = fill,
    hjust = hjust, vjust = vjust, alpha=0.99
  ) +
  geom_point(aes(x=Value.x, y=Value.y, color='skyblue')) +
  geom_richtext(size=5) +
  geom_point(color = "cornsilk", size = 1) +
  scale_color_identity() +
  scale_fill_identity() + xlim(-4,4) + ylim(-6,6) + 
  guides(alpha=FALSE) + xlab("Δmean")+ylab("Δnoise")
g  




