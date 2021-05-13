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
setwd("/Users/itoutouma/Lab_Analysis/scRNAseq")
#AllGenotype = list.files(getwd())[grep("YPD_WTvs", list.files(getwd()))] %>% 
#  sapply(function(chr){strsplit(chr, "vs")[[1]][2]}) %>% sort(decreasing = T)
#AllGenotype = AllGenotype[AllGenotype != "GCN4"]
#AllGenotype = AllGenotype[AllGenotype != "GLN3"]

NoiseOnlyMergePropList = data.frame()
MeanMergePropList = data.frame()
NoiseMergePropList = data.frame()

RTG1.target = Yeastract.interact %>% filter(V1 == "RTG1")
RTG3.target = Yeastract.interact %>% filter(V1 == "RTG3")
STP.common = intersect(RTG1.target$V2, RTG3.target$V2) %>% 
  as.data.frame() %>% set_names(c("V2"))
Common.target <- STP.common
strain.list = c("RTG1", "RTG3")
strain.list = c("GAT1", "GZF3", "DAL80")
setwd("/Users/itoutouma/Lab_Analysis/scRNAseq")
for(target in strain.list){
  print(target)
  directory = paste("YPD_WTvs", target, sep="")
  setwd(directory)
  print(getwd())
  dir = "./"
  #dirlist_0.99 = paste("All_YPD_WTvs", target, "_zh0.99", sep="")
  dirlist_0.99 = list.files(getwd())[grep("C[0-9]_YPD_WTvs", list.files(getwd()))] 
  print(dirlist_0.99)
  Cluster.id = sapply(dirlist_0.99, function(chr){return(strsplit(chr, "_")[[1]][1])})
  ClusterMeanList = list()
  ClusterResDispList = list()
  DelChainName = paste(tolower(target), "Regression", sep="")
  for(dir in dirlist_0.99){
    WT_Regression <- BASiCS_LoadChain("WT(ho)Regression", StoreDir = dir) 
    Del_Regression <- BASiCS_LoadChain(DelChainName, StoreDir = dir) 
    TestWTSTP2 <- BASiCS_TestDE(Chain1 = WT_Regression, Chain2 = Del_Regression,
                                GroupLabel1 = "WT", GroupLabel2 = target,
                                EpsilonM = log2(1.5), EpsilonD = log2(1.5), 
                                EpsilonR= log2(1.5)/log2(exp(1)),
                                EFDR_M = 0.01/length(dirlist_0.99 ), EFDR_D = 0.01/length(dirlist_0.99 ), EFDR_R = 0.01/length(dirlist_0.99 ),
                                #EFDR_M = 0.10, EFDR_D = 0.10,
                                Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)
    ResDisp = as.data.frame(TestWTSTP2, Parameter = "ResDisp", Filter=FALSE)
    ClusterResDispList = c(ClusterResDispList, list(ResDisp))
    Mean = as.data.frame(TestWTSTP2, Parameter = "Mean", Filter=FALSE)
    ClusterMeanList = c(ClusterMeanList, list(Mean))
  }
  names(ClusterResDispList) = Cluster.id
  names(ClusterMeanList) = Cluster.id
  #STP2.Target = Yeastract.interact %>% dplyr::filter(V1 == target) %>% dplyr::select(V2)
  ClusterMeanDiff= enframe(ClusterMeanList) %>% unnest() %>% 
    right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
    mutate(star = ifelse(grepl("\\+", ResultDiffMean), "*", "")) %>%
    dplyr::select(name, Name, MeanLog2FC, star, ResultDiffMean) %>% 
    setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
    mutate(comp = "Mean")
  
  ClusterNoiseDiff = enframe(ClusterResDispList) %>% unnest() %>% 
    right_join(Yeast.gene.name, by=c("GeneName"="ORF")) %>% 
    mutate(star = ifelse(grepl("\\+", ResultDiffResDisp), "*", "")) %>%
    dplyr::select(name, Name, ResDispDistance, star, ResultDiffResDisp) %>% 
    setNames(c("Group", "Gene", "Value", "star", "diff")) %>% 
    mutate(comp = "Noise")
  
  ClusterMeanDiffMerge = 
    ClusterMeanDiff %>%  
    mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
    group_by(Gene) %>% 
    dplyr::summarise(Group = target, Gene = Gene, Value = 0, diff = "Merge", comp = "Mean", starcount = sum(starnum)) %>% 
    filter(!duplicated(Gene)) %>% 
    mutate(star = ifelse(starcount > 0, "*", "")) %>% select(-starcount)
  
  ClusterNoiseDiffMerge = 
    ClusterNoiseDiff %>%  
    mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
    group_by(Gene) %>% 
    dplyr::summarise(Group = target, Gene = Gene, Value = 0, diff = "Merge", comp = "Noise", starcount = sum(starnum)) %>% 
    filter(!duplicated(Gene)) %>% 
    mutate(star = ifelse(starcount > 0, "*", "")) %>% select(-starcount)
 
  ClusterNoiseOnlyDiffMerge = 
    left_join(ClusterMeanDiff, ClusterNoiseDiff, by=c("Gene", "Group")) %>% 
    mutate(starnum = ifelse(star.x != "*" & star.y == "*", 1, 0)) %>% 
    group_by(Gene) %>% 
    dplyr::summarise(Group = target, Gene = Gene, Value = 0, diff = "Merge", comp = "Noise", starcount = sum(starnum)) %>% 
    filter(!duplicated(Gene)) %>% 
    mutate(star = ifelse(starcount > 0, "*", "")) %>% select(-starcount)
  
  MeanMergeProp = DeriveProportion(ClusterMeanDiffMerge, FALSE, Common.target)
  MeanMergePropList = rbind(MeanMergePropList,  MeanMergeProp)
  
  NoiseMergeProp = DeriveProportion(ClusterNoiseDiffMerge, FALSE, Common.target)
  NoiseMergePropList = rbind(NoiseMergePropList,  NoiseMergeProp)
  
  NoiseOnlyMergeProp = DeriveProportion(ClusterNoiseOnlyDiffMerge, FALSE, Common.target)
  NoiseOnlyMergePropList = rbind(NoiseOnlyMergePropList,  NoiseOnlyMergeProp)
  setwd("../")
}


target = "STP"

Ficher.res.pval <- Calc.Ficher.pval(NoiseOnlyMergePropList, strain.list)
g <- 
  NoiseOnlyMergePropList %>% 
  mutate(count = ifelse(count == 0, "", count)) %>% 
  ggplot(aes(fill = Category, y = prop, x = significant, alpha=alphav)) +
  facet_wrap(~Group, scale="free_y", nrow=1) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(labels=c("Known targets(Estimated)", "Unknown targets(Estimated)",
                             "Known targets(Not Estimated)", "Unknown targets(Not Estimated)"), 
                    values = c("#CB181D", "#2171B5", "#CB181D", "#2171B5")) + 
  scale_alpha(range=c(0.4,1)) + 
  geom_text(aes(label = count, y=pos, alpha=1), size = 6, hjust = 0.5) + 
  guides(alpha=FALSE) + 
  guides(fill = guide_legend(override.aes = list(alpha = c(1,1,0.1,0.1)))) + 
  geom_text(aes(x=1.5, y=1.05, label=Ficher.res.pval[Group])) + 
  ylab("Proportion") + xlab("")# +  theme(legend.position = c(0.9, 0.3))

print(g)
dir = getwd()
figtitle = file.path(dir, paste(target, "_CommonDE_NoiseOnly_MergeCluster.png", sep=""))
print(figtitle)
ggsave(figtitle, g, width=7.00, height=5.00)
#ggsave(figtitle, g, width=10.00, height=7.00)




Ficher.res.pval <- Calc.Ficher.pval(MeanMergePropList, strain.list)
g <- 
  MeanMergePropList %>% 
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
  ylab("Proportion") + xlab("") #+ theme(legend.position = c(0.9, 0.3))
print(g)
dir = "/Users/itoutouma/Lab_Analysis/scRNAseq"
figtitle = file.path(dir, paste(target, "_CommonDE_Mean_MergeCluster.png", sep=""))
print(figtitle)
#ggsave(figtitle, g, width=7.00, height=5.00)
ggsave(figtitle, g, width=10.00, height=7.00)



Ficher.res.pval <- Calc.Ficher.pval(NoiseMergePropList, strain.list)
g <- 
  NoiseMergePropList %>% 
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
  ylab("Proportion") + xlab("") #+ theme(legend.position = c(0.9, 0.3))
print(g)
dir = "/Users/itoutouma/Lab_Analysis/scRNAseq"
figtitle = file.path(dir, paste(target, "_CommonDE_Noise_MergeCluster.png", sep=""))
print(figtitle)
ggsave(figtitle, g, width=7.00, height=5.00)







#figtitle = paste("/Users/itoutouma/Lab_Analysis/scRNAseq/NosieOnlyFromClusterAllStrain.png", sep="")
figtitle = paste("/Users/itoutouma/Lab_Analysis/scRNAseq/DFDDNosieOnlyFromCluster.png", sep="")
ggsave(figtitle, g, width=13.40, height=7.00)

  