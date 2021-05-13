library(stringr)
library(plyr)
library(tidyverse)
library(scran)
library(BASiCS)
Yeastract.interact = read.csv("/Users/itoutouma/Lab_Analysis/GRN/RegulationTwoColumnTable_both.tsv", sep=';', header = FALSE)
Yeast.gene.name = read.csv("/Users/itoutouma/Lab_Analysis/GRN/orftogene.txt", sep="\t", header=T)
Yeast.gene.name = Yeast.gene.name %>% 
  mutate(SGD = ORF) %>% 
  mutate(ORF = apply(Yeast.gene.name, MARGIN=1, 
                     function(row){ if(length(grep("-", row[2]))){ return(gsub("-", "\\.", row[2]))}
                       else{ return(row[2])}}))

Calc.Ficher.pval <-
  function(PropDF, strain.list){
    sapply(strain.list, function(clust){
      tmp = PropDF %>% filter(Group == clust)
      print(tmp)
      CategoryCount = c(0, 0, 0, 0) %>% set_names(c("FT", "TT", "FF", "TF"))
      CategoryCount[unfactor(tmp$Category)] = tmp$count
      fisher.res = fisher.test(matrix(CategoryCount, nrow=2))
      print(fisher.res)
      return(fisher.res$p.value)
    }) %>% as.data.frame() %>% set_names(c("Pvalue")) %>% 
      mutate(pvalue.label = ifelse(Pvalue < 0.05, "p < 0.05", "")) %>% 
      #mutate(pvalue.label = ifelse(Pvalue < 0.01, "p < 0.01", pvalue.label)) %>%
      #mutate(pvalue.label = ifelse(Pvalue < 0.005, "p < 0.005", pvalue.label)) %>%
      #mutate(pvalue.label = ifelse(Pvalue < 0.001, "p < 0.001", pvalue.label)) %>% 
      select(pvalue.label) %>% nth(1) %>% as.array() %>% 
      set_names(strain.list) %>% return()
  }


DeriveProportion <- function(DiffResultDF, ExcludeMeanChange, interest.TF.Target){
  if(interest.TF.Target == "GATA"){
    print("GATA mode")
    interest.TF.Target = Extract.GATA.common(DiffResultDF$Group[1])
  }
  if(ExcludeMeanChange){
    print(dim(DiffResultDF))
    if(dim(DiffResultDF)[2]<6){
      print("DF don't much format")
      return()
    }else{
      DiffValueMean = DiffResultDF %>%  na.omit() %>% 
        mutate(TFTarget = ifelse(is.element(Gene, interest.TF.Target$V2), TRUE, FALSE)) %>% 
        mutate(starnum = ifelse(star.x != "*" & star.y == "*", 1, 0)) %>% 
        select(Group, Gene, Value.y, star.y, diff.y, comp.y, starnum, TFTarget) %>% 
        set_names(c("Group", "Gene", "Value", "star", "diff", "comp", "starnum", "TFTarget")) %>% 
        dplyr::group_by(TFTarget, Group) %>% 
        dplyr::summarise(signif.num = sum(starnum), nonsignif.num = length(Gene)-sum(starnum), significant = c(TRUE, FALSE)) %>% 
        mutate(count = ifelse(significant, signif.num, nonsignif.num)) %>% 
        select(Group, TFTarget, significant, count) %>% 
        mutate(id = as.numeric(str_extract_all(Group, "[0-9.]+"))+1) %>% 
        mutate(id = ifelse(Group=="All", 1, id)) %>% arrange(Group)
    }
  }else{
    DiffValueMean = DiffResultDF %>%  na.omit() %>% 
      mutate(TFTarget = ifelse(is.element(Gene, interest.TF.Target$V2), TRUE, FALSE)) %>% 
      mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
      dplyr::group_by(TFTarget, Group) %>% 
      dplyr::summarise(signif.num = sum(starnum), nonsignif.num = length(Gene)-sum(starnum), significant = c(TRUE, FALSE)) %>% 
      mutate(count = ifelse(significant, signif.num, nonsignif.num)) %>% 
      select(Group, TFTarget, significant, count) %>% 
      mutate(id = as.numeric(str_extract_all(Group, "[0-9.]+"))+1) %>% 
      mutate(id = ifelse(Group=="All", 1, id)) %>% arrange(Group)
  }
  DiffValueMean %>% 
    mutate(Category = ifelse(!TFTarget & significant, "FT", "")) %>% 
    mutate(Category = ifelse(TFTarget & significant, "TT", Category)) %>% 
    mutate(Category = ifelse(!TFTarget & !significant, "FF", Category)) %>% 
    mutate(Category = ifelse(TFTarget & !significant, "TF", Category)) %>% 
    mutate(significant = ifelse(significant, "Estimated", "NotEstimated")) %>% 
    dplyr::group_by(significant, Group) %>% 
    dplyr::summarise(Group = Group, Category = Category, significant=significant, count=count, prop = count/sum(count)) %>% 
    ddply(.(significant, Group), transform, pos = cumsum(prop) - (0.5 * prop)) %>% 
    mutate(alphav = ifelse(significant == "Estimated", 1, 0.5)) %>% 
    mutate(Category =  factor(Category, levels = c("TT", "FT", "TF", "FF"))) %>% 
    mutate(significant = factor(significant, levels = c("Estimated", "NotEstimated"))) %>% 
    return()
}

Extract.GATA.common <- function(target){
  GATAfamilyTargetList = list()
  GATAfamilyTargetList[[1]] = Yeastract.interact %>% filter(V1 == "GAT1")
  GATAfamilyTargetList[[2]] = Yeastract.interact %>% filter(V1 == "GZF3")
  GATAfamilyTargetList[[3]] = Yeastract.interact %>% filter(V1 == "DAL80")
  GATAfamilyTargetList[[4]] = Yeastract.interact %>% filter(V1 == "GLN3")
  GATAfamily = c("GAT1", "GZF3", "DAL80", "GLN3")
  names(GATAfamilyTargetList) = GATAfamily
  others = setdiff(GATAfamily, target)
  a = intersect(GATAfamilyTargetList[[target]]$V2, GATAfamilyTargetList[[others[1]]]$V2)
  b = intersect(GATAfamilyTargetList[[target]]$V2, GATAfamilyTargetList[[others[2]]]$V2)
  c = intersect(GATAfamilyTargetList[[target]]$V2, GATAfamilyTargetList[[others[2]]]$V2)
  Common.target = data.frame(V2 = unique(c(a, b, c)))
  return(Common.target)
}

#setwd("/Users/itoutouma/Lab_Analysis/scRNAseq")
strain.list = list.files(getwd())[grep("YPD_WTvs", list.files(getwd()))] %>% 
  sapply(function(chr){strsplit(chr, "vs")[[1]][2]}) %>% sort(decreasing = T)

#GCN4 is excluded because of it regulate most of whole genes
strain.list = strain.list[strain.list != "GCN4"]
ClusterResDispList = list()
ClusterMeanList = list()
for(target in strain.list){
  print(target)
  dir = paste("YPD_WTvs", target, "/All_YPD_WTvs", target, "_zh0.99", sep="")
  DelChainName = paste(tolower(target), "Regression", sep="")
  WT_Regression <- BASiCS_LoadChain("WT(ho)Regression", StoreDir = dir) 
  Del_Regression <- BASiCS_LoadChain(DelChainName, StoreDir = dir) 
  TestWTSTP2 <- BASiCS_TestDE(Chain1 = WT_Regression, Chain2 = Del_Regression,
                              GroupLabel1 = "WT", GroupLabel2 = target,
                              EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                              EpsilonR= log2(1.5)/log2(exp(1)),
                              EFDR_M = 0.01, EFDR_D = 0.01, EFDR_R = 0.01,
                              #EFDR_M = 0.10, EFDR_D = 0.10,
                              Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)
  ResDisp = as.data.frame(TestWTSTP2, Parameter = "ResDisp", Filter=FALSE) %>% 
    mutate(Group = target)
  ClusterResDispList = c(ClusterResDispList, list(ResDisp))
  Mean = as.data.frame(TestWTSTP2, Parameter = "Mean", Filter=FALSE) %>% 
    mutate(Group = target)
  ClusterMeanList = c(ClusterMeanList, list(Mean))
}


ClusterMeanList2 = 
  ClusterMeanList %>% lapply(function(df){
    df %>% right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
      mutate(star = ifelse(grepl("\\+", ResultDiffMean), "*", "")) %>%
      dplyr::select(Group, Name, MeanLog2FC, star, ResultDiffMean) %>% 
      setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
      mutate(comp = "Mean") %>% return()
  }) %>% set_names(strain.list)

ClusterNoiseList2 = 
  ClusterResDispList %>% lapply(function(df){
    df %>% right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
      mutate(star = ifelse(grepl("\\+", ResultDiffResDisp), "*", "")) %>%
      dplyr::select(Group, Name, ResDispDistance, star, ResultDiffResDisp) %>% 
      setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
      mutate(comp = "Noise") %>% return()
  }) %>% set_names(strain.list)

MeanNoiseList = lapply(names(ClusterMeanList2), 
                       function(target){
                         ClusterMeanList2[[target]] %>% 
                           left_join(ClusterNoiseList2[[target]], by=c("Gene", "Group")) %>% return()                        
                       }) %>% set_names(c(strain.list))

RTG1.target = Yeastract.interact %>% filter(V1 == "STP1")
RTG3.target = Yeastract.interact %>% filter(V1 == "STP2")
RTG.common = intersect(RTG1.target$V2, RTG3.target$V2) %>% 
  as.data.frame() %>% set_names(c("V2"))

#RTG.common = "GATA"
EnrichmentTarget = RTG.common

MeanProp = ClusterMeanList2 %>% lapply(DeriveProportion, FALSE, EnrichmentTarget) %>% 
  do.call(rbind.data.frame,.)

NoiseProp = ClusterNoiseList2 %>% lapply(DeriveProportion, FALSE, EnrichmentTarget) %>% 
  do.call(rbind.data.frame,.)

NoiseOnlyProp = MeanNoiseList %>% lapply(DeriveProportion, ExcludeMeanChange=TRUE, RTG.common) %>% 
  do.call(rbind.data.frame,.)

target = "STP"
Ficher.res.pval <- Calc.Ficher.pval(NoiseOnlyProp, strain.list)
g <- 
  NoiseOnlyProp %>% 
  mutate(count = ifelse(count == 0, "", count)) %>% 
  ggplot(aes(fill = Category, y = prop, x = significant, alpha=alphav)) +
  facet_wrap(~Group, scale="free_y", nrow=1) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(labels=c("Common targets(Estimated)", "Non common targets(Estimated)",
                             "Common targets(Not Estimated)", "Non common targets(Not Estimated)"), 
                    values = c("#CB181D", "#2171B5", "#CB181D", "#2171B5")) + 
  scale_alpha(range=c(0.4,1)) + 
  geom_text(aes(label = count, y=pos, alpha=1), size = 6, hjust = 0.5) + 
  guides(alpha=FALSE) + 
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,0.1,0.1)))) + 
  geom_text(aes(x=1.5, y=1.05, label=Ficher.res.pval[Group])) + 
  ylab("Proportion") + xlab("") + 
  ggtitle(paste()) 
print(g)
dir = getwd()
figtitle = file.path(dir, paste(target, "_CommonDE_NoiseOnly.png", sep=""))
print(figtitle)
ggsave(figtitle, g, width=5.00, height=5.00)



Ficher.res.pval <- Calc.Ficher.pval(NoiseProp, strain.list)
g <- 
  NoiseProp %>% 
  mutate(count = ifelse(count == 0, "", count)) %>% 
  ggplot(aes(fill = Category, y = prop, x = significant, alpha=alphav)) +
  facet_wrap(~Group, scale="free_y", nrow=2) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(labels=c("Common targets(Estimated)", "Non common targets(Estimated)",
                             "Common targets(Not Estimated)", "Non common targets(Not Estimated)"), 
                    values = c("#CB181D", "#2171B5", "#CB181D", "#2171B5")) + 
  scale_alpha(range=c(0.4,1)) + 
  geom_text(aes(label = count, y=pos, alpha=1), size = 6, hjust = 0.5) + 
  guides(alpha=FALSE) + 
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,0.1,0.1)))) + 
  geom_text(aes(x=1.5, y=1.05, label=Ficher.res.pval[Group])) + 
  ylab("Proportion") + xlab("") + 
  ggtitle(paste()) 
print(g)


Ficher.res.pval <- Calc.Ficher.pval(MeanProp, strain.list)
g <- 
  MeanProp %>% 
  mutate(count = ifelse(count == 0, "", count)) %>% 
  ggplot(aes(fill = Category, y = prop, x = significant, alpha=alphav)) +
  facet_wrap(~Group, scale="free_y", nrow=1) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(labels=c("Common targets(Estimated)", "Non common targets(Estimated)",
                             "Common targets(Not Estimated)", "Non common targets(Not Estimated)"), 
                    values = c("#CB181D", "#2171B5", "#CB181D", "#2171B5")) + 
  scale_alpha(range=c(0.4,1)) + 
  geom_text(aes(label = count, y=pos, alpha=1), size = 6, hjust = 0.5) + 
  guides(alpha=FALSE) + 
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,0.1,0.1)))) + 
  geom_text(aes(x=1.5, y=1.05, label=Ficher.res.pval[Group])) + 
  ylab("Proportion") + xlab("") + 
  ggtitle(paste()) 

print(g)
dir = getwd()
figtitle = file.path(dir, paste(target, "_CommonDE_mean.png", sep=""))
print(figtitle)
#ggsave(figtitle, g, width=7.00, height=5.00)
ggsave(figtitle, g, width=13.00, height=5.00)
#ggsave("GATACommonDE_mean.png", g, width=7.00, height=5.00)


Ficher.res.pval <- Calc.Ficher.pval(NoiseProp, strain.list)
g <- 
  NoiseProp %>% 
  mutate(count = ifelse(count == 0, "", count)) %>% 
  ggplot(aes(fill = Category, y = prop, x = significant, alpha=alphav)) +
  facet_wrap(~Group, scale="free_y", nrow=1) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(labels=c("Common targets(Estimated)", "Non common targets(Estimated)",
                             "Common targets(Not Estimated)", "Non common targets(Not Estimated)"), 
                    values = c("#CB181D", "#2171B5", "#CB181D", "#2171B5")) + 
  scale_alpha(range=c(0.4,1)) + 
  geom_text(aes(label = count, y=pos, alpha=1), size = 6, hjust = 0.5) + 
  guides(alpha=FALSE) + 
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,0.1,0.1)))) + 
  geom_text(aes(x=1.5, y=1.05, label=Ficher.res.pval[Group])) + 
  ylab("Proportion") + xlab("") + 
  ggtitle(paste()) 

print(g)
dir = getwd()
figtitle = file.path(dir, paste(target, "_CommonDE_Noise.png", sep=""))
print(figtitle)
ggsave(figtitle, g, width=13.00, height=5.00)
#ggsave(figtitle, g, width=7.00, height=5.00)












print(g)
figtitle = paste("MeanAllSampleAllStrain.png", sep="")
ggsave(figtitle, g, width=13.40, height=7.00)

#STP2 0.03168
#STP1 0.006876

#Extract novel target candidate genes for GO enrichment analysis
target = "STP1"
RTG1.target = Yeastract.interact %>% filter(V1 == "GLN3")
RTG3.target = Yeastract.interact %>% filter(V1 == "DAL80")
STP.common = intersect(RTG1.target$V2, RTG3.target$V2) %>% 
  as.data.frame() %>% set_names(c("V2"))

STP.common = Extract.GATA.common("GLN3")
EstimatedGene = ClusterNoiseList2[[target]] %>% filter(star == "*") %>% select(Gene)
EstimatedGene = ClusterNoiseDiffMerge %>% filter(star == "*") %>% select(Gene)
Known.target = Yeastract.interact %>% filter(V1 == target)
Novel.gene = EstimatedGene %>% filter(!is.element(Gene, STP.common$V2)) 
Novel.gene %>% write.table("STP1NovelRedundantFromNoiseCluster.csv", row.names = F, col.names = F)



