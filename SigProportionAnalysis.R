if(!exists('Mean.Noise.Bind.df')){cat("Execute CreateHeatmap_segment_CS5.R first!")}

DiffValueMean = Mean.Noise.Bind.df %>%  na.omit() %>% 
  mutate(STP2Target = ifelse(is.element(Gene, STP2BindGenes$V2), TRUE, FALSE)) %>% 
  mutate(starnum = ifelse(star == "*", 1, 0)) %>% 
  dplyr::group_by(STP2Target, Group, comp) %>% dplyr::summarise(starnum = sum(starnum), num = length(Gene)) 

#Test
for(cons in c("All", "C1", "C2", "C3", "C4", "C5")){
  print(cons)
  tmpdf <- DiffValueMean %>% filter(comp=="Noise") %>% filter(Group == cons)
  TG = tmpdf %>% filter(STP2Target)
  NoTG = tmpdf %>% filter(!STP2Target)
  fisher.test(matrix(c(TG$starnum, NoTG$starnum, (TG$num-TG$starnum), (NoTG$num-NoTG$starnum)),nrow=2)) %>% 
    print()
}
#＊　p<0.05、＊＊　P<0.01、　＊＊＊P<0.005
DiffValueMean$Group <- factor(DiffValueMean$Group, levels = c("All", "C1", "C2", "C3", "C4", "C5"))
g1 <- ggplot(filter(DiffValueMean, comp=="Mean"), aes(x = Group, y = starnum/num, fill = STP2Target)) + 
  ylab("Proportion of significance") + 
  geom_bar(stat = "identity", position = "dodge") + 
  ylim(0,0.21) + theme_classic(base_size = 10, base_family = "Helvetica") +
  scale_fill_manual(values=c("gray10","grey70")) + 
  geom_text(x = 5,  y = 0.14,  label = "*", colour = "black") +
  geom_segment(x=4.7, xend=5.3, y=0.13, yend=0.13) + 
  xlab("") + ggtitle("Mean")
g1
#ggsave("MeanPropSig.png", g1, width=5.6, height=4.0)

g2 <- ggplot(filter(DiffValueMean, comp=="Noise"), aes(x = Group, y = starnum/num, fill = STP2Target)) + 
  ylab("Proportion of significance") + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_classic(base_size = 10, base_family = "Helvetica") +
  ylim(0,0.21)+
  scale_fill_manual(values=c("gray10","grey70")) + 
  ##All label
  geom_text(x = 1,  y = 0.06,  label = "***", colour = "black") +
  geom_segment(x=0.7, xend=1.3, y=0.05, yend=0.05) + 
  #C1 label
  geom_text(x = 2,  y = 0.03,  label = "*", colour = "black") +
  geom_segment(x=1.7, xend=2.3, y=0.02, yend=0.02) + 
  #C2 label
  geom_text(x = 3,  y = 0.04,  label = "***", colour = "black") +
  geom_segment(x=2.7, xend=3.3, y=0.03, yend=0.03) + 
  #C3 is not significant
  #C4 label
  geom_text(x = 5,  y = 0.07,  label = "***", colour = "black") +
  geom_segment(x=4.7, xend=5.3, y=0.06, yend=0.06) + 
  #C5 label
  geom_text(x = 6,  y = 0.04,  label = "***", colour = "black") +
  geom_segment(x=5.7, xend=6.3, y=0.03, yend=0.03) + 
  xlab("") + ggtitle("Noise")
g2
#ggsave("NoisePropSig.png", g2, width=5.6, height=4.0)




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
