Fisher.exact <- function(DiffValueMean, target){
  
  Ficher.res.pval <-
    sapply(Cluster.id, function(clust){
    tmp = DiffValueMean %>% filter(Group == clust)
    fisher.res = fisher.test(matrix(tmp$count, nrow=2))
    # [Sig, nonsig] x [down, nondown]
    print(fisher.res)
    return(fisher.res$p.value)
  })
  print(Ficher.res.pval)
  
  TFtarget.TRUE <- filter(DiffValueMean, TFTarget) %>% 
    group_by(Group) %>% arrange(significant) %>% 
    mutate(pos = cumsum(count) - count / 2) %>% 
    mutate(pos = ifelse(significant, pos+150, pos)) %>% 
    mutate(signif2 = ifelse(significant, "A", "B"))
  TFtarget.TRUE$Group <- factor(TFtarget.TRUE$Group, levels = Cluster.id)
  TFtarget.TRUE$significant <- factor(TFtarget.TRUE$significant, levels=c(TRUE, FALSE))
  TFtarget.FALSE <- filter(DiffValueMean, !TFTarget) %>% 
    group_by(Group) %>% arrange(significant) %>% 
    mutate(pos = cumsum(count) - count / 2) %>% 
    mutate(signif2 = ifelse(significant, "C", "D"))
  TFtarget.FALSE$Group <- factor(TFtarget.FALSE$Group, levels = Cluster.id)
  TFtarget.FALSE$significant <- factor(TFtarget.FALSE$significant, levels=c(TRUE, FALSE))
  barwidth = 0.35
  Ficher.res.pval.df = 
    data.frame(xpos = unique(TFtarget.FALSE$id + barwidth/2),
               label = Ficher.res.pval, label2 = c(rep(NA,6))) %>% 
    mutate(label2 = ifelse(label < 0.05/5, "padj\n*", label2)) %>%
    mutate(label2 = ifelse(label < 0.01/5, "padj\n**", label2)) %>% 
    mutate(label2 = ifelse(label < 0.005/5, "padj\n***", label2)) %>% 
    mutate(label2 = ifelse(label < 0.001/5, "padj\n****", label2))
  Ficher.res.pval.df[1,] = 
    Ficher.res.pval.df[1,] %>%   
    mutate(label2 = ifelse(label < 0.05, "*", label2)) %>%
    mutate(label2 = ifelse(label < 0.01, "**", label2)) %>% 
    mutate(label2 = ifelse(label < 0.005, "***", label2)) %>% 
    mutate(label2 = ifelse(label < 0.001, "****", label2))
  fig.title = paste("Left: ", target, "TG ", "Right: Non", target, "TG", sep="")
  g1 <- 
    ggplot() + 
    geom_bar(data = TFtarget.TRUE, 
             mapping = aes(x = Group, y = count, fill = significant), 
             stat="identity", 
             position='stack', 
             width = barwidth) + 
    geom_text(data = TFtarget.TRUE, size=4,
              aes(x = Group, y = pos, label = count)) + 
    geom_bar(data = TFtarget.FALSE, 
             mapping = aes(x = id+ barwidth + 0.01, y = count, fill = significant), 
             stat="identity", 
             position='stack' , 
             width = barwidth) + 
    geom_text(data = TFtarget.FALSE, size=4,
              aes(x = id + barwidth + 0.01, y = pos, label = count)) + 
    labs(color="Category")+labs(shape="Category")+
    xlab("Cell group") + ylab("#gene") + ggtitle(paste(comparison, fig.title)) +  
    guides(size=FALSE) + guides(fill=FALSE) + 
    scale_fill_manual(values=c("#fafad2", "#d3d3d3")) +
    theme_bw()
  #Add point
  g <-
    g1 + 
    geom_point(data = TFtarget.TRUE, size=5,
               mapping = aes(x = Group, y = pos-100, shape=signif2, color=signif2)) + 
    geom_point(data = TFtarget.FALSE, size=5,
               mapping = aes(x = id+ barwidth + 0.01, y = pos-100, shape=signif2, color=signif2)) +
    scale_shape_manual(labels=c("Estimated(In Prior)", "Not Estimated(In Prior)\n", "Estimated(Not In Prior)", "Not Estimated(Not In Prior)"), values=c(8, 2, 8, 2)) + 
    scale_color_manual(labels=c("Estimated(In Prior)", "Not Estimated(In Prior)\n", "Estimated(Not In Prior)", "Not Estimated(Not In Prior)"), values=c("red", "red", "blue", "blue")) + 
    geom_text(data = Ficher.res.pval.df,
              aes(x = xpos, y = -250, label = label2, size=4)) 
  
  return(g)
}

DeriveProportion <- function(comparison, target){
  interest.TF.Target = Yeastract.interact %>% dplyr::filter(V1 == target)
  if(comparison == "bind"){
    DiffValueMean = Mean.Noise.Bind.df %>%  na.omit() %>% 
      mutate(starnum.split = ifelse(star == "*", 1, 0)) %>% 
      dplyr::group_by(Group, Gene) %>% 
      dplyr::summarise(starnum = sum(starnum.split)) %>% 
      mutate(starnum = ifelse(starnum >= 1, 1, 0)) %>% 
      mutate(TFTarget = ifelse(is.element(Gene, interest.TF.Target$V2), TRUE, FALSE)) %>% 
      dplyr::group_by(TFTarget, Group) %>%
      dplyr::summarise(signif.num = sum(starnum), nonsignif.num = length(Gene)-sum(starnum), significant = c(TRUE, FALSE)) %>% 
      mutate(count = ifelse(significant, signif.num, nonsignif.num)) %>% 
      select(Group, TFTarget, significant, count) %>% 
      mutate(id = as.numeric(str_extract_all(Group, "[0-9.]+"))+1) %>% 
      mutate(id = ifelse(Group=="All", 1, id)) %>% arrange(Group)
  }else{
    print("No bind")
    DiffValueMean = Mean.Noise.Bind.df %>%  na.omit() %>% 
      filter(comp==comparison) %>% 
      mutate(TFTarget = ifelse(is.element(Gene, interest.TF.Target$V2), TRUE, FALSE)) %>% 
      mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
      dplyr::group_by(TFTarget, Group) %>% 
      dplyr::summarise(signif.num = sum(starnum), nonsignif.num = length(Gene)-sum(starnum), significant = c(TRUE, FALSE)) %>% 
      mutate(count = ifelse(significant, signif.num, nonsignif.num)) %>% 
      select(Group, TFTarget, significant, count) %>% 
      mutate(id = as.numeric(str_extract_all(Group, "[0-9.]+"))+1) %>% 
      mutate(id = ifelse(Group=="All", 1, id)) %>% arrange(Group)
  }
  return(DiffValueMean)
}
genotype = "STP1"
comparison = "Mean"
DiffValueMean = DeriveProportion("Mean", genotype)
g <- Fisher.exact(DiffValueMean, genotype)
print(g)
#ggsave("FicherTypeBarMeanBASiCS2.png", g, width=7.00, height=6.80)

DiffValueMean = DeriveProportion("Noise", genotype)
comparison = "Noise"
g <- Fisher.exact(DiffValueMean, genotype)
print(g)
#ggsave("FicherTypeBarNoise.png", g, width=7.00, height=6.80)

DiffValueMean = DeriveProportion("bind", genotype)
g <- Fisher.exact(DiffValueMean, genotype)
print(g)
