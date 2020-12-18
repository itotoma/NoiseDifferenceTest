STP2BindGenes

hist_df = data.frame()
DataList = MeanList
for(i in c(1:length(DataList))){
  df <- DataList[[i]]
  media <- strsplit(names(DataList)[i], "_")[[1]][1]
  print(media)
  commondf = left_join(df, yeast_gene_name, by=c("GeneName"="ORF")) %>% 
    select(Name, ResultDiffMean)
  AllAnalyzed = commondf %>% filter(is.element(Name,STP2BindGenes$V2))
  Lacks = nrow(STP2BindGenes) - nrow(AllAnalyzed)
  InferredPrior = AllAnalyzed %>% filter(ResultDiffMean != "NoDiff") %>% nrow()
  NonInferredPrior = nrow(AllAnalyzed) - InferredPrior
  tmpdf <- rbind(c(media, "Lack", Lacks),c(media, "Detected", InferredPrior), c(media, "Not Detected", NonInferredPrior))
  hist_df = rbind(hist_df, tmpdf)
}

colnames(hist_df) = c("group", "subgroup", "value")
hist_df$value = as.numeric(hist_df$value)
hist_df$subgroup2 <- factor(hist_df$subgroup, levels = c("Detected", "Not Detected"))
g <- ggplot(hist_df, aes(x = group, y = value, fill = subgroup2)) 
g <- g + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_blank())
g <- g + ggtitle("Mean difference test")
g

