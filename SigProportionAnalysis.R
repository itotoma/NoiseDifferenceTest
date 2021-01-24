if(!exists('Mean.Noise.Bind.df')){cat("Execute CreateHeatmap_segment_CS5.R first!")}

interest.TF.Target = Yeastract.interact %>% dplyr::filter(V1 == "RTG1")
  
DiffValueMean = Mean.Noise.Bind.df %>%  na.omit() %>% 
  mutate(TFTarget = ifelse(is.element(Gene, interest.TF.Target$V2), TRUE, FALSE)) %>% 
  mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
  dplyr::group_by(TFTarget, Group, comp) %>% dplyr::summarise(starnum = sum(starnum), num = length(Gene)) 

Create.pval.bar <- function(df, ylist, margin, plist, ymax, title){
  non.signif.index <- c(1:6)[plist > 0.05]
  plist[plist < 0.005] = "p < 0.005"
  plist[plist < 0.01] = "p < 0.01"
  plist[plist < 0.05] = "p < 0.05"
  xcor = c(1:6)
  xcor[non.signif.index] = -10
  g1 <- ggplot(df, aes(x = Group, y = prop, fill = TFTarget)) + 
    ylab("Proportion of significance") + 
    geom_bar(stat = "identity", position = "dodge") + 
    ylim(0,ymax) + theme_classic(base_size = 10, base_family = "Helvetica") +
    scale_fill_manual(values=c("gray10","grey70")) + 
    geom_text(x = xcor[1],  y = ylist[1]+margin,  label = plist[1], colour = "black") +
    geom_text(x = xcor[2],  y = ylist[2]+margin,  label = plist[2], colour = "black") +
    geom_text(x = xcor[3],  y = ylist[3]+margin,  label = plist[3], colour = "black") +
    geom_text(x = xcor[4],  y = ylist[4]+margin,  label = plist[4], colour = "black") +
    geom_text(x = xcor[5],  y = ylist[5]+margin,  label = plist[5], colour = "black") +
    geom_text(x = xcor[6],  y = ylist[6]+margin,  label = plist[6], colour = "black") +
    geom_segment(x=xcor[1]-0.3, xend=xcor[1]+0.3, y=ylist[1], yend=ylist[1]) + 
    geom_segment(x=xcor[2]-0.3, xend=xcor[2]+0.3, y=ylist[2], yend=ylist[2]) + 
    geom_segment(x=xcor[3]-0.3, xend=xcor[3]+0.3, y=ylist[3], yend=ylist[3]) + 
    geom_segment(x=xcor[4]-0.3, xend=xcor[4]+0.3, y=ylist[4], yend=ylist[4]) + 
    geom_segment(x=xcor[5]-0.3, xend=xcor[5]+0.3, y=ylist[5], yend=ylist[5]) + 
    geom_segment(x=xcor[6]-0.3, xend=xcor[6]+0.3, y=ylist[6], yend=ylist[6]) +
    xlab("") + ggtitle(title)
  return(g1)
}


compair = "Noise"
fisher.res.list = sapply(c("All", "C1", "C2", "C3", "C4", "C5"), 
                         function(ci){
                           tmpdf <- DiffValueMean %>% 
                             dplyr::filter(comp==compair) %>% 
                             dplyr::filter(Group == ci)
                           TG = tmpdf %>% dplyr::filter(TFTarget)
                           NoTG = tmpdf %>% dplyr::filter(!TFTarget)
                           fisher.res <- fisher.test(matrix(c(TG$starnum, NoTG$starnum, (TG$num-TG$starnum), (NoTG$num-NoTG$starnum)),nrow=2))
                           print(fisher.res)
                           return(fisher.res$p.value)
                         })
#＊　p<0.05、＊＊　P<0.01、　＊＊＊P<0.005
test.df <- DiffValueMean %>% dplyr::filter(comp==compair) %>% mutate(prop = starnum/num)
test.df$Group = factor(test.df$Group, levels = c("All", "C1", "C2", "C3", "C4", "C5"))
ylist = test.df %>% dplyr::group_by(Group) %>% summarise(max = max(prop))
g1 <- Create.pval.bar(test.df, ylist = ylist$max+0.01, 
                margin=0.01, fisher.res.list, ymax = 0.25, title=compair)
print(g1)
#ggsave("NoisePropSig.png", g1, width=5.6, height=4.0)


compair = "Mean"
fisher.res.list = sapply(c("All", "C1", "C2", "C3", "C4", "C5"), 
                         function(ci){
                           tmpdf <- DiffValueMean %>% 
                             dplyr::filter(comp==compair) %>% 
                             dplyr::filter(Group == ci)
                           TG = tmpdf %>% dplyr::filter(TFTarget)
                           NoTG = tmpdf %>% dplyr::filter(!TFTarget)
                           fisher.res <- fisher.test(matrix(c(TG$starnum, NoTG$starnum, (TG$num-TG$starnum), (NoTG$num-NoTG$starnum)),nrow=2))
                           print(fisher.res)
                           return(fisher.res$p.value)
                         })
#＊　p<0.05、＊＊　P<0.01、　＊＊＊P<0.005
test.df <- DiffValueMean %>% dplyr::filter(comp==compair) %>% mutate(prop = starnum/num)
test.df$Group = factor(test.df$Group, levels = c("All", "C1", "C2", "C3", "C4", "C5"))
ylist = test.df %>% dplyr::group_by(Group) %>% summarise(max = max(prop))
g1 <- Create.pval.bar(test.df, ylist = ylist$max+0.01, 
                margin=0.01, fisher.res.list, ymax = 0.25, title=compair)
print(g1)
ggsave("MeanPropSig.png", g1, width=5.6, height=4.0)


#Analyze accordance of results among clusters
ShannonEntropy <- function(num1, num2){
  num3 = 1-num1-num2
  if(num1 == 0 & num2 != 0){
    return(-1*(num2*log(num2)+(1-num2)*log(1-num2)))
  }else if(num1 != 0 & num2 == 0){
    return(-1*(num1*log(num1)+(1-num1)*log(1-num1)))
  }else if(num1 == 0 & num2 == 0){
    return(0)
  }else{
    return(-1*(num1*log(num1)+num2*log(num2)+num3*log(num3)))
  }
}
SignificanceEntropy =
  df %>% na.omit() %>% dplyr::group_by(Gene, comp) %>% 
  mutate(WTsig = ifelse(star == "*" & Value > 0, 1, 0)) %>%
  mutate(STP2sig = ifelse(star == "*" & Value < 0, 1, 0)) %>%
  dplyr::summarise(SE = ShannonEntropy(sum(WTsig)/6, sum(STP2sig)/6),
                   WTp = sum(WTsig)/6, STP2p = sum(STP2sig)/6)
ggplot(SignificanceEntropy, aes(x = Gene, y = "", fill = SE)) + 
  geom_tile() + theme_bw() + facet_wrap(~comp, scale="free_y", nrow=2) +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        #axis.text.x = element_blank()
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=4.5)
  ) + ggtitle("Shannon Entropy")


ggplot(SignificanceEntropy, aes(x=comp, y=SE)) + geom_boxplot() + 
  xlab("") + ylab("Shannon entropy")
