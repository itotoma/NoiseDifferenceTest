library(BASiCS)
#Reference http://bioconductor.org/packages/release/bioc/vignettes/BASiCS/inst/doc/BASiCS.html#alternative-implementation-modes

################################################################################
#Construct BASiCS_Chain object from single cell expreriment object
################################################################################
set.seed(1)
CountsNoSpikes1 <- matrix(rpois(50*40, 2), ncol = 40)
rownames(CountsNoSpikes1) <- paste0("Gene", 1:50)
DataExampleNoSpikes1 <- SingleCellExperiment(assays = list(counts = CountsNoSpikes1), 
                                            colData = data.frame(BatchInfo = rep(c(1,2), each = 20)))
set.seed(2)
CountsNoSpikes2 <- matrix(rpois(50*40, 2), ncol = 40)
rownames(CountsNoSpikes2) <- paste0("Gene", 1:50)
DataExampleNoSpikes2 <- SingleCellExperiment(assays = list(counts = CountsNoSpikes2), 
                                             colData = data.frame(BatchInfo = rep(c(1,2), each = 20)))
#RECOMMENDED SETTING IS N=20000, Thin=20 and Burn=10000 (it takes long time)

#Please ensure the acceptance rates displayed in the console 
#output of BASiCS_MCMC are around 0.44. If they are too far from 
#this value, you should increase N and Burn

#BASiCS assumes that a pre-processing quality control step has been applied

ChainRegression1 <- BASiCS_MCMC(Data = DataExampleNoSpikes1, N = 1000, 
                               Thin = 10, Burn = 500, 
                               Regression = TRUE , 
                               PrintProgress = FALSE,
                               WithSpikes = FALSE,
                               StoreChains = TRUE,  # This three options
                               StoreDir = getwd(),  # are required if you
                               RunName = "Example") # want to export 

ChainRegression2 <- BASiCS_MCMC(Data = DataExampleNoSpikes2, N = 1000, 
                                Thin = 10, Burn = 500, 
                                Regression = TRUE , 
                                PrintProgress = FALSE,
                                WithSpikes = FALSE)
TestRegression <- BASiCS_TestDE(Chain1 = ChainRegression1, Chain2 = ChainRegression2,
                      GroupLabel1 = "SC", GroupLabel2 = "PaS",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EFDR_M = 0.10, EFDR_D = 0.10,
                      Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)

class(ChainRegression1)
class(TestRegression)
#BASiCS_Chain object

#Decompose total variability of gene expression into biological and technical components.
#The regression is considered when you create the BASiCS_Chain
BASiCS_VarianceDecomp(ChainRegression1)

########################################
#Detech highly or lowely variable gene 
########################################

#For regression data
#EpsilonThreshold: Threshold for residual overdispersion
DetectHVG <- BASiCS_DetectHVG(ChainRegression1, EpsilonThreshold = log(2), Plot = TRUE)
DetectLVG <- BASiCS_DetectLVG(ChainRegression1, EpsilonThreshold = -log(2), Plot = TRUE)
as.data.frame(DetectHVG)
as.data.frame(DetectLVG)

#Non regression data
data(ChainSC)
HVG <- BASiCS_DetectHVG(ChainSC, VarThreshold = 0.6)
LVG <- BASiCS_DetectLVG(ChainSC, VarThreshold = 0.2)
BASiCS_PlotVG(HVG, "VG")
BASiCS_PlotVG(LVG, "VG")

################################################################################
#Differential over-dispersion tests between two cell populations
################################################################################
data(ChainSC)
data(ChainRNA)
Test <- BASiCS_TestDE(Chain1 = ChainSC, Chain2 = ChainRNA,
                      GroupLabel1 = "SC", GroupLabel2 = "PaS",
                      EpsilonM = log2(1.5), EpsilonD = log2(1.5),
                      EFDR_M = 0.10, EFDR_D = 0.10,
                      Offset = TRUE, PlotOffset = TRUE, Plot = TRUE)
#EpsilonM threshold for mean
#EpsilonD threshold for over-dispersion
#EpsilonR threshold for residual over-dispersion

#Specify mean, dispersion and residual overdispersion
head(as.data.frame(Test, Parameter = "Mean"))
head(as.data.frame(Test, Parameter = "Disp"))
BASiCS_PlotDE(Test, Plots = "MA", Parameters = "Mean")
BASiCS_PlotDE(Test, Plots = "MA", Parameters = "Disp")

data("ChainRNAReg")
BASiCS_ShowFit(ChainRNAReg)



