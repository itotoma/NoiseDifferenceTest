dirlist_0.99 <- c("YPD_WTsSTP2_zh0.99", "YPD_WTsSTP2_zh0.99_c1", "YPD_WTsSTP2_zh0.99_c2",
                  "YPD_WTsSTP2_zh0.99_c3", "YPD_WTsSTP2_zh0.99_c4", "YPD_WTsSTP2_zh0.99_c5")
ClusterMeanList = list()
ClusterResDispList = list()
ClusterDispList = list()

for(dir in dirlist_0.99){
  WT_Regression <- BASiCS_LoadChain("WT(HO)Regression", StoreDir = dir) 
  STP2_Regression <- BASiCS_LoadChain("RTG1Regression", StoreDir = dir) 
  TestWTSTP2 <- BASiCS_TestDE(Chain1 = WT_Regression, Chain2 = STP2_Regression,
                              GroupLabel1 = "WT", GroupLabel2 = "STP2",
                              EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                              EFDR_M = 0.10, EFDR_D = 0.10,
                              Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)
  ResDisp = as.data.frame(TestWTSTP2, Parameter = "ResDisp", Filter=FALSE)
  ClusterResDispList = c(ClusterResDispList, list(ResDisp))
  Disp = as.data.frame(TestWTSTP2, Parameter = "Disp", Filter=FALSE)
  ClusterDispList = c(ClusterDispList, list(Disp))
  Mean = as.data.frame(TestWTSTP2, Parameter = "Mean", Filter=FALSE)
  ClusterMeanList = c(ClusterMeanList, list(Mean))
}
if(0){source(CreateHeatmap_segment.R)}
names(ClusterResDispList) = c("All", "C1", "C2", "C3", "C4", "C5")
names(ClusterDispList) = c("All", "C1", "C2", "C3", "C4", "C5")
names(ClusterMeanList) = c("All", "C1", "C2", "C3", "C4", "C5")




ClusterResDispDf = enframe(ClusterResDispList) %>% unnest() %>% na.omit()
LowResDispTh <- quantile(filter(ClusterDispDf, name=="All")$Disp1)[[2]]
HighResDispTh <- quantile(filter(ClusterDispDf, name=="All")$Disp1)[[4]]
ClusterDispDf2 = ClusterDispDf %>% 
  mutate(cate = ifelse(Disp1 <= LowResDispTh, "Low", "Mid")) %>% 
  mutate(cate = ifelse(Disp1 >= HighResDispTh, "High", cate))
ClusterDispDf2$cate = factor(ClusterDispDf2$cate, levels = c("High", "Mid", "Low"))
g1 <- ClusterDispDf2 %>% 
  ggplot(aes(x=name, y=Disp1, color=name)) + geom_boxplot() + 
  facet_wrap(~cate, scale="free_y", nrow=3) +
  #geom_jitter(size = 0.8, alpha=0.1) +
  scale_color_manual(values = c("brown", "tomato", "yellow4", "mediumaquamarine", "skyblue2", "darkorchid1")) +
  ylab("ResDisp")+xlab("Cluster group") + 
  theme(legend.position = 'none')
g1
#0 proportion がこれだと覗かれている


ClusterID = c(0:5) %>% set_names(c("All","C1","C2","C3","C4","C5"))
Multi.count.stat <-
  lapply(c("All","C1","C2","C3","C4","C5"), function(ci){
    if(ci == "All"){
      clusterdf = Multi.Count.Matrix
    }else{
      clusterdf =  Multi.Count.Matrix[,filter(clusteredDF, Cluster == ClusterID[[ci]])$rowname] 
    }
    data.frame(var = rowVars(clusterdf), sum = rowSums(clusterdf), mean = rowMeans(clusterdf),Gene = rownames(clusterdf)) %>%  
      mutate(cluster = ci)  %>% return()
  }) %>% do.call(rbind.data.frame,.)


g <- 
  Multi.count.stat %>% 
  mutate(cat = ifelse(var>=quantile(Multi.count.stat$var)[[3]],"high","mid" )) %>% 
  mutate(cat = ifelse(var<=quantile(Multi.count.stat$var)[[2]], "low", cat)) %>% 
  ggplot(aes(x=cluster, y=var/mean)) + geom_boxplot() + 
  facet_wrap(~cat, scale="free_y", nrow=3)
g




#GO slim mapper
groupname = c("All", "C1", "C2", "C3", "C4", "C5") %>% set_names(dirlist_0.99)
GOSlimRes = lapply(dirlist_0.99, function(dir){
  read.csv(file.path(dir, "GenericGOslimProcess.txt"), sep="\t") %>% 
    select(GOID, TERM, CLUSTER_FREQUENCY, TERM, NUM_LIST_ANNOTATIONS) %>% 
    mutate(CLUSTER_FREQUENCY = gsub("%", "", CLUSTER_FREQUENCY)) %>% 
    mutate(Group = groupname[dir]) %>% 
    mutate(CLUSTER_FREQUENCY = as.numeric(CLUSTER_FREQUENCY)) %>% 
    mutate(NUM_LIST_ANNOTATIONS = as.numeric(NUM_LIST_ANNOTATIONS)) %>% 
    mutate(TERMLENGTH = nchar(TERM)) %>% 
    return()}
) %>% do.call('rbind.data.frame',.) %>% arrange(TERMLENGTH)
#[is.element(GOSlimRes$GOID, unique(GOSlimRes$GOID)[1:27]),]

GOSlimRes$TERM = factor(GOSlimRes$TERM, levels = unique(GOSlimRes$TERM))
g1 <- 
  GOSlimRes %>%  
  ggplot(aes(x = TERM, y = CLUSTER_FREQUENCY, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme(legend.position = c(0.01, 0.9), legend.justification = c(0, 1))  +
  #theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 48, hjust = 1)) + 
  scale_fill_manual(values = c("brown", "tomato", "yellow4", "mediumaquamarine", "skyblue2", "darkorchid1")) + 
  xlab("") + ylab("#HVG Genes")
#ggsave("GenericGOslimProcess.png", g1, width = 12.0, height = 6.0)
g1
