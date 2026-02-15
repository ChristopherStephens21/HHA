#------DORC Gene Dynamics Analysis-----#
#In this script, we will analyze DORC gene dynamics across pseudotime for the cortex and IRS cuticle

#-----Loading Packages-----#
library(Seurat)
library(scCustomize)
library(scDblFinder)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(presto)
library(AUCell)
library(clusterProfiler)
library(ggpubr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Bulb_Seurat_Plots/Figure_6/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#----Helper Functions----#
#Extracts raw count matrix and runs AUCell
RunAUCell <- function(SeuratObj, GeneLists, assay = "RNA") {
  #----Calculating AUCell Scores-----#
  #Pulling out raw count data
  cat("Pulling Count Matrix \n") #Modifying to support Assay Argument 
  SeuratObj.counts <- GetAssayData(object = SeuratObj, slot = "counts", assay = assay)
  #computing gene rank order by cell
  cat("Computing Gene Rank Order by Cell \n")
  SeuratObj.rankings <- AUCell_buildRankings(SeuratObj.counts, plotStats=F)
  cat("Scoring Regulons \n")
  #computing auc scores for each set
  SeuratObj.scores <- AUCell_calcAUC(GeneLists, SeuratObj.rankings, aucMaxRank = ceiling(0.05 * nrow(SeuratObj.rankings)))
  #converting to cell x gene DF
  SeuratObj.AUCresults <- getAUC(SeuratObj.scores) %>% t() %>% as.data.frame()
  colnames(SeuratObj.AUCresults) <- names(GeneLists) %>% str_replace_all(pattern = " ", replacement = "_")
  return(SeuratObj.AUCresults)}
#Custom function for Ordering Genes by timepoint of max expression. 
OrderGenes <- function(Smoothed) {
  colnames(Smoothed) <- 1:ncol(Smoothed)
  #Getting max expression values per gene
  MaxExp <- rowMaxs(Smoothed)
  names(MaxExp) <- rownames(Smoothed)
  #Getting row indices for max expression values
  MaxIndices <- c()
  for (i in 1:length(MaxExp)) {
    #Pulling smoothed expression into vector
    SmoothedExp <- Smoothed[i, ]
    MaxIndex <- SmoothedExp[SmoothedExp == MaxExp[i]]
    MaxIndices <- c(MaxIndices, MaxIndex)}
  Sorting.DF <- data.frame(Gene = rownames(Smoothed), MaxIndex = MaxIndices)
  GeneOrder <- Sorting.DF %>% arrange(MaxIndex) %>% pull(Gene)
  return(GeneOrder)}

#--------------Plotting--------------#
#------Palettes-----#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
#From MetBrewer::met.brewer("Signac", 10) "#92c051" "#2b9b81" "#1f6e9c" "#633372" "#9f5691" "#e87b89" "#de597c" "#d8443c" "#fe9b00" "#f4c40f"
SampleIDPal <- c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f")
Metadata.Pal <- scale_fill_manual(values = SampleIDPal)

#--------Loading Checkpoint-------#
#load(paste0(SeuratDirectory, "HHA_Matrix_DORC_Gene_Dynamics_2_1_26.RData"))

#--------Loading Image------#
load(paste0(SeuratDirectory, file = "HHA_Matrix_DORC_Scoring_1_7_26.RData"))

#----------Examining Gene Trends Across Pseudotime: Cortex-------#
#Reading in Pseudotime Results: Cortex
#Pseudotime
Cortex.Pseudotime <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_Cortex_Cuticle_Pseudotime.csv")  %>% 
  select(barcode = "barcode...1", "palantir_pseudotime_CoCu", "palantir_entropy_CoCu")
Cortex.Pseudotime
#Branch Assignments 
Cortex.Branches <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_Cortex_Cuticle_Branches.csv") %>% 
  select(barcode = "barcode...1", "Cortex_Terminal", "Cuticle_Terminal")
Cortex.BranchMetadata <- Cortex.Branches %>% mutate(Barcode = barcode) %>% left_join(HHA.Matrix[[]], by = "Barcode")
Cortex.BranchMetadata
#Selecting Multiome nuclei within the IRS Cuticle Branch
Cortex.Branch <- Cortex.BranchMetadata %>% filter(Cortex_Terminal) %>% filter(Platform == "Multiome") %>% select(barcode, UnifiedBarcode)
#Filtering Pseudotime for cells in branch and ordering by Pseudotime
Cortex.BranchPseudotime <- Cortex.Pseudotime %>% filter(barcode %in% Cortex.Branch$barcode == T) %>% arrange(palantir_pseudotime_CoCu)
Cortex.BranchPseudotime #3,338 nuclei
Cortex.BranchCells <- Cortex.BranchPseudotime$barcode

#------Pulling Matrix Associated with the IRS Trajectory-----#
colnames(Matrix.Multiome) <- Matrix.Multiome$Barcode

#----------Examining Gene Trends Across Pseudotime: IRS Cuticle-------#
#Reading in Pseudotime Results: IRS
#Pseudotime
IRS.Pseudotime <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_IRS_Pseudotime.csv")
IRS.Pseudotime
#Branch Assignments 
IRS.Branches <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_IRS_Branches.csv")
IRS.BranchMetadata <- IRS.Branches %>% mutate(Barcode = barcode) %>% left_join(HHA.Matrix[[]], by = "Barcode")
IRS.BranchMetadata
#Selecting Multiome nuclei within the IRS Cuticle Branch
IRSCuticle.Branch <- IRS.BranchMetadata %>% filter(IRS_Cuticle_Terminal == T) %>% 
  filter(Platform == "Multiome") %>% select(barcode, UnifiedBarcode)
#Filtering Pseudotime for cells in branch and ordering by Pseudotime
IRSCuticle.BranchPseudotime <- IRS.Pseudotime %>% filter(barcode %in% IRSCuticle.Branch$barcode == T) %>% arrange(palantir_pseudotime_IRS)
IRSCuticle.BranchPseudotime #1650 nuclei
IRSCuticle.BranchCells <- IRSCuticle.BranchPseudotime$barcode
#Selecting Multiome nuclei within the IRS Cuticle Branch
IRSCuticle.Branch <- IRS.BranchMetadata %>% filter(IRS_Cuticle_Terminal == T) %>% 
  filter(Platform == "Multiome") %>% select(barcode, UnifiedBarcode)
#Filtering Pseudotime for cells in branch and ordering by Pseudotime
IRSCuticle.BranchPseudotime <- IRS.Pseudotime %>% filter(barcode %in% IRSCuticle.Branch$barcode == T) %>% arrange(palantir_pseudotime_IRS)
IRSCuticle.BranchPseudotime #987 nuclei
IRSCuticle.BranchCells <- IRSCuticle.BranchPseudotime$barcode

#------Strategy 2: Min-Max Normalize Expression and DORCs across all genes, then plot against Pseudotime-------#
DefaultAssay(Matrix.Multiome) <- "RNA"
#-------Gene Expression-----#
TotalCells <- c(Cortex.BranchCells, IRSCuticle.BranchCells)
#Grabbing Log Normalized Expression
Total.GeneMat <- Matrix.Multiome %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
#Getting Expression for all genes across all cells
Total.GeneMat <- Total.GeneMat[HRG.Genes, TotalCells] %>% t() %>% as.data.frame()
str(Total.GeneMat)
#MinMax Normalizing Expression across all cells
MinMaxNormalize <- function(x) {(x - min(x)) / (max(x) - min(x))}
TotalGeneMat.Norm <- as.data.frame(lapply(Total.GeneMat, MinMaxNormalize))
str(TotalGeneMat.Norm)
rownames(TotalGeneMat.Norm) <- rownames(Total.GeneMat)
colMaxs(as.matrix(TotalGeneMat.Norm))
#Smoothing Expression along Cortex/Cuticle Trajectory
CortexGeneMat.Norm <- TotalGeneMat.Norm[Cortex.BranchCells, ]
str(CortexGeneMat.Norm)
#Spline Smoothing Matrix
CortexGeneMat.NormSmoothed <- as.data.frame(lapply(CortexGeneMat.Norm, function(x){smooth.spline(x,df=3)$y}))
#Loess Smoothing Expression
#CortexGeneMat.NormSmoothed <- as.data.frame(lapply(CortexGeneMat.Norm, function(x){predict(loess(x ~ Cortex.BranchPseudotime$palantir_pseudotime_CoCu, span = 0.5))}))
#Smoothing Expression along IRS Trajectory
IRSCuticleGeneMat.Norm <- TotalGeneMat.Norm[IRSCuticle.BranchCells, ]
str(IRSCuticleGeneMat.Norm)
#Spline Smoothing Matrix
IRSCuticleGeneMat.NormSmoothed <- as.data.frame(lapply(IRSCuticleGeneMat.Norm, function(x){smooth.spline(x,df=3)$y}))
#Loess Smoothing Expression
#IRSCuticleGeneMat.NormSmoothed <- as.data.frame(lapply(IRSCuticleGeneMat.Norm, function(x){predict(loess(x ~ IRSCuticle.BranchPseudotime$palantir_pseudotime_IRS, span = 0.5))}))
str(IRSCuticleGeneMat.NormSmoothed)

#--------DORC Scores--------#
#Grabbing DORC Scores
Total.DORCMat <- Matrix.Multiome %>% JoinLayers() %>% SeuratObject::LayerData(assay = "DORC_Scores", layer = 'data')
#Subsetting for Cells on the Cortex trajectory, keeping all genes
Total.DORCMat <- Total.DORCMat[HRG.Genes, TotalCells] %>% t() %>% as.data.frame()
str(Total.DORCMat)
#MinMax Normalizing Accessibility across all cells
TotalDORCMat.Norm <- as.data.frame(lapply(Total.DORCMat, MinMaxNormalize))
str(TotalDORCMat.Norm)
rownames(TotalDORCMat.Norm) <- rownames(Total.DORCMat)
#Smoothing Expression along Cortex/Cuticle Trajectory
CortexDORCMat.Norm <- TotalDORCMat.Norm[Cortex.BranchCells, ]
str(CortexDORCMat.Norm)
#Spline Smoothing Matrix
CortexDORCMat.NormSmoothed <- as.data.frame(lapply(CortexDORCMat.Norm, function(x){smooth.spline(x,df=3)$y}))
#Loess Smoothing Expression
#CortexDORCMat.NormSmoothed <- as.data.frame(lapply(CortexDORCMat.Norm, function(x){predict(loess(x ~ Cortex.BranchPseudotime$palantir_pseudotime_CoCu, span = 0.5))}))
str(CortexDORCMat.NormSmoothed)
#Smoothing Expression along IRS Trajectory
IRSCuticleDORCMat.Norm <- TotalDORCMat.Norm[IRSCuticle.BranchCells, ]
str(IRSCuticleDORCMat.Norm)
#Spline Smoothing Matrix
IRSCuticleDORCMat.NormSmoothed <- as.data.frame(lapply(IRSCuticleDORCMat.Norm, function(x){smooth.spline(x,df=3)$y}))
#Loess Smoothing Expression
#IRSCuticleDORCMat.NormSmoothed <- as.data.frame(lapply(IRSCuticleDORCMat.Norm, function(x){predict(loess(x ~ IRSCuticle.BranchPseudotime$palantir_pseudotime_IRS, span = 0.5))}))
str(IRSCuticleDORCMat.NormSmoothed)

#--------Computing Lags: Cortex---------#
#Min Max Normalizing Smoothed Curves: Gene Expression
CortexGeneMat.NormSmoothedScaled <- as.data.frame(lapply(CortexGeneMat.NormSmoothed, MinMaxNormalize)) %>% t() %>% as.matrix()
str(CortexGeneMat.NormSmoothedScaled)
rowMaxs(CortexGeneMat.NormSmoothedScaled)
#Min Max Normalizing Smoothed Curves: DORC Accessibility
CortexDORCMat.NormSmoothedScaled <- as.data.frame(lapply(CortexDORCMat.NormSmoothed, MinMaxNormalize)) %>% t() %>% as.matrix()
str(CortexDORCMat.NormSmoothedScaled)
rowMaxs(CortexDORCMat.NormSmoothedScaled)
#Calculating Lag
Cortex.LagMat <- CortexDORCMat.NormSmoothedScaled - CortexGeneMat.NormSmoothedScaled
Cortex.LagMat <- as.data.frame(t(Cortex.LagMat))
rownames(Cortex.LagMat) <- rownames(CortexDORCMat.NormSmoothed)
str(Cortex.LagMat)

#--------Computing Lags: IRS Cuticle---------#
#Min Max Normalizing Smoothed Curves: Gene Expression
IRSCuticleGeneMat.NormSmoothedScaled <- as.data.frame(lapply(IRSCuticleGeneMat.NormSmoothed, MinMaxNormalize)) %>% t() %>% as.matrix()
str(IRSCuticleGeneMat.NormSmoothedScaled)
rowMaxs(IRSCuticleGeneMat.NormSmoothedScaled)
#Min Max Normalizing Smoothed Curves: DORC Accessibility
IRSCuticleDORCMat.NormSmoothedScaled <- as.data.frame(lapply(IRSCuticleDORCMat.NormSmoothed, MinMaxNormalize)) %>% t() %>% as.matrix()
str(IRSCuticleDORCMat.NormSmoothedScaled)
rowMaxs(IRSCuticleDORCMat.NormSmoothedScaled)
#Calculating Lag
IRSCuticle.LagMat <- IRSCuticleDORCMat.NormSmoothedScaled - IRSCuticleGeneMat.NormSmoothedScaled
IRSCuticle.LagMat <- as.data.frame(t(IRSCuticle.LagMat))
rownames(IRSCuticle.LagMat) <- rownames(IRSCuticleDORCMat.NormSmoothed)
str(IRSCuticle.LagMat)

#----------Getting Combined Results-------#
#Normalizing trajectories to same starts and ends
Cortex.NormPseudotime <- MinMaxNormalize(Cortex.BranchPseudotime$palantir_pseudotime_CoCu)
IRSCuticle.NormPseudotime <- MinMaxNormalize(IRSCuticle.BranchPseudotime$palantir_pseudotime_IRS)

#---------Gene Expression---------#
#Combining Gene Expression Across both clusters in tidy format
Combined.NormGenes <- rbind(CortexGeneMat.NormSmoothed, IRSCuticleGeneMat.NormSmoothed)
#Min-Max Normalizing Across Trajectories
Combined.NormGenes <- as.data.frame(lapply(Combined.NormGenes, MinMaxNormalize))
str(Combined.NormGenes)
#Adding Trajectory Information
Combined.NormGenes$Trajectory <- c(rep("Cortex", nrow(CortexGeneMat.NormSmoothed)), rep("IRS Cuticle", nrow(IRSCuticleGeneMat.NormSmoothed)))
Combined.NormGenes$Pseudotime <- c(Cortex.NormPseudotime, IRSCuticle.NormPseudotime)
Combined.NormGenes$Barcode <- TotalCells

#-------DORC Accessibility-------#
#Combining DORC Accessibility Across both clusters in tidy format
Combined.NormDORCs <- rbind(CortexDORCMat.NormSmoothed, IRSCuticleDORCMat.NormSmoothed)
#Min-Max Normalizing Across Trajectories
Combined.NormDORCs <- as.data.frame(lapply(Combined.NormDORCs, MinMaxNormalize))
str(Combined.NormDORCs)
#Adding Trajectory Information
Combined.NormDORCs$Trajectory <- c(rep("Cortex", nrow(CortexDORCMat.NormSmoothed)), rep("IRS Cuticle", nrow(IRSCuticleDORCMat.NormSmoothed)))
Combined.NormDORCs$Pseudotime <- c(Cortex.NormPseudotime, IRSCuticle.NormPseudotime)
Combined.NormDORCs$Barcode <- TotalCells

#----------Formatting for Spatial Plotting------#
#-------Cortex------#
#Gene Expression
Cortex.PlottingGenes <- Combined.NormGenes %>% select(!c(Trajectory, Pseudotime)) %>% column_to_rownames("Barcode")
Cortex.PlottingGenes <- Cortex.PlottingGenes[Cortex.BranchCells,]
colnames(Cortex.PlottingGenes) <- paste0(colnames(Cortex.PlottingGenes), "_GEX")
#DORC Accessibility
Cortex.PlottingDORCs <- Combined.NormDORCs %>% select(!c(Trajectory, Pseudotime)) %>% column_to_rownames("Barcode")
Cortex.PlottingDORCs <- Cortex.PlottingDORCs[Cortex.BranchCells,]
colnames(Cortex.PlottingDORCs) <- paste0(colnames(Cortex.PlottingDORCs), "_DORC")
#Latency
Cortex.PlottingLags <- Cortex.LagMat
rownames(Cortex.PlottingLags) <- Cortex.BranchCells
colnames(Cortex.PlottingLags) <- paste0(colnames(Cortex.PlottingLags), "_Lag")
#Binding into 1 DF
Cortex.Plotting <- cbind(Cortex.PlottingGenes, Cortex.PlottingDORCs, Cortex.PlottingLags) %>% 
  mutate(Pseudotime = Cortex.NormPseudotime, Barcode = rownames(Cortex.PlottingLags))
#Exporting for Plotting
write_csv(Cortex.Plotting, file = "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Cortex_Gene_Trends_Pseudotime.csv")
Cortex.Metadata <- read_csv("/projects/b1217/HHA/Bulb_Spatial/Cortex_Multiome_Subset.csv") %>% column_to_rownames("barcode...1")
Cortex.Metadata
Cortex.PlottingReordered <-Cortex.Plotting[rownames(Cortex.Metadata),]
write_csv(Cortex.PlottingReordered, file = "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Cortex_Gene_Trends_Pseudotime_Reordered.csv")

#-------IRS Cuticle------#
#Gene Expression
IRSCuticle.PlottingGenes <- Combined.NormGenes %>% select(!c(Trajectory, Pseudotime)) %>% column_to_rownames("Barcode")
IRSCuticle.PlottingGenes <- IRSCuticle.PlottingGenes[IRSCuticle.BranchCells,]
colnames(IRSCuticle.PlottingGenes) <- paste0(colnames(IRSCuticle.PlottingGenes), "_GEX")
#DORC Accessibility
IRSCuticle.PlottingDORCs <- Combined.NormDORCs %>% select(!c(Trajectory, Pseudotime)) %>% column_to_rownames("Barcode")
IRSCuticle.PlottingDORCs <- IRSCuticle.PlottingDORCs[IRSCuticle.BranchCells,]
colnames(IRSCuticle.PlottingDORCs) <- paste0(colnames(IRSCuticle.PlottingDORCs), "_DORC")
#Latency
IRSCuticle.PlottingLags <- IRSCuticle.LagMat
rownames(IRSCuticle.PlottingLags) <- IRSCuticle.BranchCells
colnames(IRSCuticle.PlottingLags) <- paste0(colnames(IRSCuticle.PlottingLags), "_Lag")
#Binding into 1 DF
IRSCuticle.Plotting <- cbind(IRSCuticle.PlottingGenes, IRSCuticle.PlottingDORCs, IRSCuticle.PlottingLags) %>% 
  mutate(Pseudotime = IRSCuticle.NormPseudotime, Barcode = rownames(IRSCuticle.PlottingLags))
#Exporting for Plotting
write_csv(IRSCuticle.Plotting, file = "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/IRSCuticle_Gene_Trends_Pseudotime.csv")
IRSCuticle.Metadata <- read_csv("/projects/b1217/HHA/Bulb_Spatial/IRS_Multiome_Subset.csv") %>% column_to_rownames("barcode...1")
IRSCuticle.Metadata
IRSCuticle.PlottingReordered <-IRSCuticle.Plotting[rownames(IRSCuticle.Metadata),]
write_csv(IRSCuticle.PlottingReordered, file = "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/IRSCuticle_Gene_Trends_Pseudotime_Reordered.csv")

#--------Plotting Combined Results-------#
#---------GATA3---------#
GATA3.DORCPlot <- ggplot(Combined.NormDORCs, aes(x = Pseudotime, y = GATA3, color = Trajectory)) + geom_line() +
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("DORC Accessibility")
GATA3.DORCPlot
GATA3.GEXPlot <- ggplot(Combined.NormGenes, aes(x = Pseudotime, y = GATA3, color = Trajectory)) + geom_line() + 
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("Gene Expression")
GATA3.GEXPlot
GATA3.LagPlot <-ggplot(IRSCuticle.Plotting, aes(x = Pseudotime, y = GATA3_Lag)) + geom_line(color = "darkseagreen4") + 
  theme_bw() + labs(title = "'Chromatin-RNA Latency", y = "GATA3")
GATA3.LagPlot

#---------KRT36---------#
KRT36.DORCPlot <- ggplot(Combined.NormDORCs, aes(x = Pseudotime, y = KRT36, color = Trajectory)) + geom_line() +
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("DORC Accessibility")
KRT36.DORCPlot
KRT36.GEXPlot <- ggplot(Combined.NormGenes, aes(x = Pseudotime, y = KRT36, color = Trajectory)) + geom_line() + 
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("Gene Expression")
KRT36.GEXPlot
KRT36.LagPlot <-ggplot(Cortex.Plotting, aes(x = Pseudotime, y = KRT36_Lag)) + geom_line(color = "darkseagreen4") + 
  theme_bw() + labs(title = "'Chromatin-RNA Latency", y = "KRT36")
KRT36.LagPlot

#---------TP63---------#
TP63.DORCPlot <- ggplot(Combined.NormDORCs, aes(x = Pseudotime, y = TP63, color = Trajectory)) + geom_line() +
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("DORC Accessibility")
TP63.DORCPlot
TP63.GEXPlot <- ggplot(Combined.NormGenes, aes(x = Pseudotime, y = TP63, color = Trajectory)) + geom_line() + 
  scale_color_manual(values = c("cornflowerblue", "firebrick")) + theme_bw() + ggtitle("Gene Expression")
TP63.GEXPlot
TP63.LagPlot <-ggplot(Cortex.Plotting, aes(x = Pseudotime, y = TP63_Lag)) + geom_line(color = "darkseagreen4") + 
  theme_bw() + labs(title = "'Chromatin-RNA Latency", y = "TP63")
TP63.LagPlot

#-------Writing csvs-------#
#DORCs
write_csv(Combined.NormDORCs, file = "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/DORC_Trends.csv")
#Genes
write_csv(Combined.NormGenes, file = "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/Gene_Trends.csv")
#Lags
write_csv(Cortex.Plotting, file = "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/Cortex_Gene_Dynamics_Trends.csv")
write_csv(IRSCuticle.Plotting, file = "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/IRS_Cuticle_Gene_Dynamics_Trends.csv")



#----------Examining Gene Dynamics Trends Across Pseudotime: Cortex-------#
#------------Plotting Corticle Cell Distribution Along Pseudotime----------#
Cortex.PseudotimeMetadata <- Cortex.BranchPseudotime %>% left_join(Cortex.BranchMetadata, by = "barcode")
Cortex.PseudotimeMetadata
#Plotting Cluster Pseudotime Distributions
ggplot(Cortex.PseudotimeMetadata, aes(x = palantir_pseudotime_CoCu, y = MatrixAnnotationFine, color = MatrixAnnotationFine)) + ggbeeswarm::geom_quasirandom() + theme_bw() +
  scale_color_manual(values = c("#3BBCA8", "#9ECAE1", "#4292C6", "#08306B")) + NoLegend() +
  labs(x = 'Pseudotime', y = "Annotation", title = "Cortical Pseudotime Distributions") + theme(plot.title = element_text(face = "bold", hjust = 0.5))
#Plotting Entropy Predictions over time
ggplot() + geom_point(data = Cortex.PseudotimeMetadata, aes(x = palantir_pseudotime_CoCu, y = palantir_entropy_CoCu, color = MatrixAnnotationFine), alpha = 0.5) + theme_bw() +
  scale_color_manual(values = c("#3BBCA8", "#9ECAE1", "#4292C6", "#08306B")) + geom_smooth(data = Cortex.PseudotimeMetadata, aes(x = palantir_pseudotime_CoCu, y = palantir_entropy_CoCu), color = "red")

#--------------Grouping Cortex into Deciles----------------#
#Using DE test between first and last deciles to identify basal and late differentiation genes
Cortex.Terciles <- quantile(Cortex.PseudotimeMetadata$palantir_pseudotime_CoCu, c(0.1, 0.9))
Cortex.PseudotimeMetadata$CortexPseudotimeGrouping <- case_when(Cortex.PseudotimeMetadata$palantir_pseudotime_CoCu <= Cortex.Terciles[1] ~ "Early",
                                                                Cortex.PseudotimeMetadata$palantir_pseudotime_CoCu > Cortex.Terciles[1] & Cortex.PseudotimeMetadata$palantir_pseudotime_CoCu <= Cortex.Terciles[2]  ~ "Middle",
                                                                Cortex.PseudotimeMetadata$palantir_pseudotime_CoCu > Cortex.Terciles[2] ~ "Late") %>% factor(levels = c("Early", "Middle", 'Late'))
table(Cortex.PseudotimeMetadata$CortexPseudotimeGrouping)
#Checking Results
ggplot(Cortex.PseudotimeMetadata, aes(x = palantir_pseudotime_CoCu, y = CortexPseudotimeGrouping, color = MatrixAnnotationFine)) + ggbeeswarm::geom_quasirandom() + theme_bw() +
  scale_color_manual(values = c("#3BBCA8", "#9ECAE1", "#4292C6", "#08306B")) + NoLegend() +
  labs(x = 'Pseudotime', y = "Tercile", title = "Cortical Pseudotime Distributions") + theme(plot.title = element_text(face = "bold", hjust = 0.5))

#--------------Formatting for Pseudotime Differential Expression Analysis----------------#
#Subsetting for cells along Cortex Branch
Cortex.BranchSub <- Matrix.Multiome %>% JoinLayers() %>% subset(Barcode %in% Cortex.BranchCells)
Cortex.BranchSub 
#Ordering to Match Multiome Data
Cortex.BranchMetadataOrdered <- as.data.frame(Cortex.PseudotimeMetadata) %>% column_to_rownames("Barcode") 
Cortex.BranchMetadataOrdered <- Cortex.BranchMetadataOrdered[Cortex.BranchSub$Barcode,]
identical(Cortex.BranchMetadataOrdered$UnifiedBarcode,colnames(Cortex.BranchSub))
#Adding Pseudotime and Pseudotime Grouping
Cortex.BranchSub$palantir_pseudotime_CoCu <- Cortex.BranchMetadataOrdered$palantir_pseudotime_CoCu
Cortex.BranchSub$CortexPseudotimeGrouping<- Cortex.BranchMetadataOrdered$CortexPseudotimeGrouping
DimPlot(Cortex.BranchSub, group.by = "CortexPseudotimeGrouping", reduction = "UMAP_WNN") 

#--------------Running Differential Expression Analysis--------------#
DefaultAssay(Cortex.BranchSub) <- 'RNA'
Idents(Cortex.BranchSub) <- "CortexPseudotimeGrouping"
#Running DE
Cortex.EarlyVsLate <- FindMarkers(Cortex.BranchSub, ident.1 = "Early", ident.2 = "Late")
Cortex.EarlyVsLate$gene <- rownames(Cortex.EarlyVsLate)
Cortex.EarlyvLHRGs <-  Cortex.EarlyVsLate %>% filter(gene %in% HRG.Genes) %>% mutate(Early = avg_log2FC > 1.75 & p_val_adj < 0.01,
                                                                                     Late = avg_log2FC < -1.75 & p_val_adj < 0.01,
                                                                                     PercentDiff = abs(pct.1 - pct.2)) %>%
  filter(Early == T | Late == T)
Cortex.EarlyvLHRGs
#Filtering for Early Genes
Cortex.EarlyHRGs <- Cortex.EarlyvLHRGs %>% filter(Early == T) %>% pull(gene)
Cortex.EarlyHRGs %>% sort()
Cortex.LateHRGs <- Cortex.EarlyvLHRGs %>% filter(Late == T) %>% pull(gene)
Cortex.LateHRGs %>% sort()
#Pulling genes for Heatmap
Cortex.HeatmapGenes <- c(Cortex.EarlyHRGs, Cortex.LateHRGs)
Cortex.HeatmapGenes <- Cortex.HeatmapGenes[Cortex.HeatmapGenes %in% rownames(CortexGeneMat.NormSmoothedScaled)]

#-------------Getting Heatmaps for Early/Late DORC Genes-------------#
#Plotting Lagged Scores Across Pseudotime
Cortex.GeneHeatmat <- as.matrix(CortexGeneMat.NormSmoothedScaled)[Cortex.HeatmapGenes,]
Cortex.DORCHeatmat <- as.matrix(CortexDORCMat.NormSmoothedScaled)[Cortex.HeatmapGenes,]
Cortex.LagHeatmat <- t(as.matrix(Cortex.LagMat))[Cortex.HeatmapGenes,]
  
#--------------Clustering Based on GeneExp, DORCAcc, and Lag Across Pseudotime--------#
#Pulling everything to one matrix
Cortex.Dynamics <- cbind(Cortex.GeneHeatmat, Cortex.DORCHeatmat, Cortex.LagHeatmat)
#Splitting early/late genes
Heatmap.GeneSplits <- factor(ifelse(rownames(Cortex.GeneHeatmat) %in% Cortex.EarlyHRGs, "Early", "Late"), levels = c("Early", "Late"))
#Clustering rows using Cortex.Dynamics, not the heatmap matrix
cluster_on_dynamics <- function(mat) {
  idx <- match(rownames(mat), rownames(Cortex.Dynamics))
  hc  <- stats::hclust(stats::dist(Cortex.Dynamics[idx, , drop = FALSE]), method = "ward.D2")
  hc}

#------------Building Heatmap------------#
#Gene Expression
Full.Heatmap <- Heatmap(Cortex.DORCHeatmat,
                        name = "DORC Accessibility",
                        row_split = Heatmap.GeneSplits,        
                        cluster_rows = cluster_on_dynamics,   
                        cluster_row_slices = TRUE,             
                        row_dend_reorder = FALSE,
                        show_row_dend = TRUE,                
                        cluster_columns = FALSE,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        column_title = "DORC Accessibility",
                        heatmap_legend_param = list(title = "DORC Accessibility"),
                        col = colorRamp2(seq(0, 1, length=9), viridisLite::magma(9, direction=1))) +
  #DORC Accessibility
  Heatmap(Cortex.GeneHeatmat,
          name = "Gene Expression",
          row_split = Heatmap.GeneSplits,
          cluster_rows = FALSE,                 
          cluster_columns = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          column_title = "Gene Expression",
          heatmap_legend_param = list(title = "Gene Expression"),
          col = colorRamp2(seq(0, 1, length=9), viridisLite::magma(9, direction=1))) +
  #DORC-RNA Lag
  Heatmap(Cortex.LagHeatmat,
          name = "Residual",
          row_split = Heatmap.GeneSplits,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = T,
          row_names_gp = gpar(fontsize = 4),
          show_column_names = FALSE,
          column_title = "Residual",
          heatmap_legend_param = list(title = "Residual"),
          col = colorRamp2(seq(-0.5, 0.5, length=9), rev(RColorBrewer::brewer.pal(9, "RdBu"))))
#Displaying
draw(Full.Heatmap, padding = unit(c(1, 1, 1, 1), "cm"))

#---------------Drawing with Selected Row Labels-------------------#
Heatmap.LabelGenes <- c("CNN3", "TGM3", 'FOXO3', "KRT83", "PSORS1C2", "WNT3", "KRT36", "IVL", #Late
                        "FGFR2", "IGFBP2", "CHSY1", "WNT5A", "KRT5", 'TP63', 'OAF') #Early
#Drawing to learn final row order per slice 
FullHeatmap.Temp <- draw(Full.Heatmap)
#Verifying identical row ordering 
identical(row_order(FullHeatmap.Temp, "DORC Accessibility"), row_order(FullHeatmap.Temp, "Gene Expression"))
identical(row_order(FullHeatmap.Temp, "DORC Accessibility"), row_order(FullHeatmap.Temp, "Residual"))
#Pulling row order from drawn heatmap 
Heatmap.RowOrder <- row_order(FullHeatmap.Temp, "DORC Accessibility")
Heatmat.RowNames  <- rownames(Cortex.DORCHeatmat)
#Building slice-aware positions for label annotation
Heatmap.IndexList <- lapply(Heatmap.RowOrder, function(index) index[Heatmat.RowNames[index] %in% Heatmap.LabelGenes])
Heatmap.Index <- unlist(Heatmap.IndexList, use.names = FALSE)
Heatmap.LabelList <- lapply(Heatmap.RowOrder, function(index) Heatmat.RowNames[index][Heatmat.RowNames[index] %in% Heatmap.LabelGenes])
Heatmap.Label <- unlist(Heatmap.LabelList, use.names = FALSE)
#Ordering 
Heatmap.AnnotOrder <- order(Heatmap.Index)
Heatmap.Index  <- Heatmap.Index[Heatmap.AnnotOrder]
Heatmap.Label <- Heatmap.Label[Heatmap.AnnotOrder]
#Rebuilding with Annotated Rows
Full.HeatmapLabeled <- Full.Heatmap + rowAnnotation(mark = anno_mark(at = Heatmap.Index, labels = Heatmap.Label, labels_gp = gpar(fontsize = 6), link_gp   = gpar(lwd = 0.5)))
draw(Full.HeatmapLabeled, padding = unit(c(1, 1, 1, 1), "cm"))

#----------------Saving as pdf------------------#
HeatmapDirectory <- "/projects/b1217/HHA/Revisions/Matrix_DORC_Heatmaps/"
pdf(paste0(HeatmapDirectory, "Matrix_GeneDynamics_Heatmap_Labeled.pdf"), width = 10, height = 8, useDingbats = FALSE)
draw(Full.HeatmapLabeled, padding = unit(c(1,1,1,1), "cm"))
dev.off()

#-------------------Saving Workspace Image---------------------#
save.image(paste0(SeuratDirectory, "HHA_Matrix_DORC_Gene_Dynamics_2_1_26.RData"))


