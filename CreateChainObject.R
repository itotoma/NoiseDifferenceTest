library(BASiCS)
library(tidyverse)

Batch <- c(1:length(unique(FilteredData$Replicate))) %>% 
         set_names(unique(FilteredData$Replicate))
CountRep <- FilteredData %>% select(Replicate) %>% 
              nth(1) %>% sapply(function(x){return(Batch[x])})

FilteredCount = FilteredData %>% select(-Replicate) %>% as.matrix() %>% t()
SCE <- SingleCellExperiment(assays = list(counts = FilteredCount), 
                               colData = data.frame(BatchInfo = CountRep))

RunName = "STP2Regression"
Regression <- BASiCS_MCMC(Data = SCE, N = 20000, 
                             Thin = 20, Burn = 10000, 
                             Regression = TRUE , 
                             PrintProgress = TRUE,
                             WithSpikes = FALSE,
                             StoreChains = TRUE,
                             StoreDir = getwd(),
                             RunName = RunName)


