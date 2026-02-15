#-------Clustering and Integration of HHA scRNA-seq Data-----#
#This is part 1 of our single cell workflow, comprising the initial QC filtering. 
#We will load in our data and use linear mixture modeling (implemented with miQC) to identify putative living/dead cells in each dataset. 
#This is separated from the rest of the workflow due to issues with environment compatability with Seurat. 
#Input data are count matrices corrected by CellBender. 

#-----Loading in Packages----#
library(Seurat) #scRNA-seq
library(EnsDb.Hsapiens.v86) #Genome for ClusterProfiler
library(clusterProfiler) #Gene set overrepresentation analysis 
library(tidyverse) #Data handling
library(scCustomize) #Reading in CellBender output
library(scDblFinder) #Doublet Removal 
library(flexmix) #Linear mixture modeling
library(miQC) #Linear mixture modeling for detection of dead cells 

#-----Directories-----#
SeuratDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/SeuratData/"
PlotsDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Plots/"

#-----Wrapper Functions-----# 
#Adds mitochondrial and ribosomal read percentages
AddMitoRibo <- function(SeuratObj) {
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-")
  SeuratObj[["percent.ribo"]] <- PercentageFeatureSet(SeuratObj, pattern = "^RP[SL]")
  return(SeuratObj)
}
#Parses the metadata file and adds relevant dataset information, converting numerical variables to factors  
AddDatasetMetadata <- function(SeuratObj, Metadata = Dataset.Metadata, Sample = SeuratObj@project.name) {
  Sample.Metadata <- Dataset.Metadata %>% filter(SampleID == Sample)
  SeuratObj[["SampleID"]] <-  Sample.Metadata$SampleID
  SeuratObj[["DonorID"]] <- as.factor(Sample.Metadata$DonorID)
  SeuratObj[["Platform"]] <- factor(Sample.Metadata$Platform, levels = c("V3.1", "V4"))
  SeuratObj[["Isolation"]] <- factor(Sample.Metadata$Isolation, levels = c("4mm", "FUE", "Bulb"))
  SeuratObj[["Sex"]] <- factor(Sample.Metadata$Sex)
  SeuratObj[["Age"]] <- as.factor(Sample.Metadata$Age)
  SeuratObj[["Source"]] <- Sample.Metadata$Source
  return(SeuratObj)
}
#Plots in viridis color schemes
PlotGenesViridis <- function(SeuratObj, reduction = "UMAP", genes, viridis_option = "viridis", ...) {
  FeaturePlot(SeuratObj, reduction = reduction, features = genes, ...) & scale_color_viridis_c(option = viridis_option) & 
    labs(x = "UMAP 1", y = "UMAP 2")
}
#Loads in data from CellBender outputted h5 files
SeurAuto_h5 <- function(Path, project, ...) {
  SeuratObj <- scCustomize::Read_CellBender_h5_Mat(file_name = Path) %>%
    CreateSeuratObject(project = project, min.cells = 3, ...) %>%
    AddMitoRibo()
  return(SeuratObj)
}
#Prepares Seurat Object for flexmix linear mixture modeling with a miQC wrapper. 
PrepForModeling <- function(SeuratObj) {
  SCE <- as.SingleCellExperiment(SeuratObj)
  #These cannot be specified in the function without an error so must be named explicitly.
  colData(SCE)$detected <- colData(SCE)$nFeature_RNA
  colData(SCE)$subsets_mito_percent <- colData(SCE)$percent.mt
  return(SCE)
}

#Generates an SCE object from a Seurat count matrix and runs scDblFinder using the default parameters (clustering workflow)
RunSCDblFinder <- function(SeuratObj) {
  set.seed(1701)
  #Cluster specific information is not really worth keeping as it is not cell-type annotated and varies across datasets. 
  KeptCols <- c("scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted", "scDblFinder.difficulty", "scDblFinder.cxds_score")
  #Running scDblFinder using default parameters
  SDF <-scDblFinder(GetAssayData(SeuratObj, layer ="counts", assay = "RNA"), clusters = T)
  SeuratObj@meta.data <- cbind(SeuratObj@meta.data, as.data.frame(colData(SDF)[, KeptCols])) #Adding output to Seurat Metadata
  return(SeuratObj)
}

#-----Palettes-----#
SampleIDPal <- MetBrewer::met.brewer("Signac", 12)

#-----Dataset Metadata-----#
#Relevant metadata for each (scRNAseq) sample included in atlas
Dataset.Metadata <- read_csv("/projects/p31989/HHA_Spatial_SC_Atlas/Dataset_Metadata/HHA_scRNAseq_Metadata_9_11_24.csv")

#-----Datasets-----#
EL_A5 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_A5/EL_A5_cellbender_corr_matrix_filtered.h5" #FUE
EL_A6 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_A6/EL_A6_cellbender_corr_matrix_filtered.h5" #4mm Punch
EL_B6 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_B6/EL_B6_cellbender_corr_matrix_filtered.h5" #FUE
EL_B7 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_B7/EL_B7_cellbender_corr_matrix_filtered.h5" #FUE 
#These two datasets were resequenced due to lower than expected count depth, and reaggregated with the initial sequencing run
EL_B8 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_B8_RESEQ/EL_B8_RESEQ_cellbender_corr_matrix_filtered.h5" #Bulb 
EL_E1 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_E1_RESEQ/EL_E1_RESEQ_cellbender_corr_matrix_filtered.h5" #4mm Punch
EL_E3 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_E3/EL_E3_cellbender_corr_matrix_filtered.h5" #4mm Punch
EL_E5 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_E5/EL_E5_cellbender_corr_matrix_filtered.h5" #Bulb 
EL_H3 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_H3/EL_H3_cellbender_corr_matrix_filtered.h5" #4mm Punch
EL_H5 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_H5/EL_H5_cellbender_corr_matrix_filtered.h5" #4mm Punch
EL_H6 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/EL_H6/EL_H6_cellbender_corr_matrix_filtered.h5" #4mm Punch

#-----Creating Seurat Objects and adding relevant metadata-----#
#200 min features excludes the high number of dead/empty droplets included by CellBender
EL_A5 <- SeurAuto_h5(EL_A5, project = "EL_A5", min.features = 200) %>% AddDatasetMetadata()
EL_A6 <- SeurAuto_h5(EL_A6, project = "EL_A6", min.features = 200) %>% AddDatasetMetadata()
EL_B6 <- SeurAuto_h5(EL_B6, project = "EL_B6", min.features = 200) %>% AddDatasetMetadata()
EL_B7 <- SeurAuto_h5(EL_B7, project = "EL_B7", min.features = 200) %>% AddDatasetMetadata()
EL_B8 <- SeurAuto_h5(EL_B8, project = "EL_B8", min.features = 200) %>% AddDatasetMetadata()
EL_E1 <- SeurAuto_h5(EL_E1, project = "EL_E1", min.features = 200) %>% AddDatasetMetadata()
EL_E3 <- SeurAuto_h5(EL_E3, project = "EL_E3", min.features = 200) %>% AddDatasetMetadata()
EL_E5 <- SeurAuto_h5(EL_E5, project = "EL_E5", min.features = 200) %>% AddDatasetMetadata()
EL_H3 <- SeurAuto_h5(EL_H3, project = "EL_H3", min.features = 200) %>% AddDatasetMetadata()
EL_H5 <- SeurAuto_h5(EL_H5, project = "EL_H5", min.features = 200) %>% AddDatasetMetadata()
EL_H6 <- SeurAuto_h5(EL_H6, project = "EL_H6", min.features = 200) %>% AddDatasetMetadata()

#-----Mixture Modeling for Detection of Dead/Permeabilized Cells-----#
#-----EL_A5-----#
#Setting seed 
set.seed(1701) #Setting seed each time to prevent accidental reproducibility issues
#Calculating LMM 
EL_A5.sce <- PrepForModeling(EL_A5)
EL_A5.lm <- mixtureModel(EL_A5.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_A5.LMMPlot <- plotModel(EL_A5.sce, EL_A5.lm) & ggsci::scale_color_gsea() & geom_hline(yintercept = 12.5, linetype = "dashed")  & theme_bw()  & ggtitle("EL_A5: Posterior Probabilities") |
  plotFiltering(EL_A5.sce, EL_A5.lm) & geom_hline(yintercept = 12.5, linetype = "dashed") & theme_bw() & ggtitle("EL_A5: Filtering at 75% Posterior Probability Threshold")
EL_A5.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_A5_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_A5.LMMPlot
dev.off()
#Adding to Seurat Object
#miQC does not label clusters identified from mixture model, so sometimes posteriors are computed for the live or dead cluster. The 1 - posterior corrects for this so all posteriors are calculated for the dead group. 
EL_A5@meta.data <- EL_A5@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_A5.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75) 
#Verifying Correct Assignment Direction 
EL_A5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_A5$Likely_Ruptured_075)
(EL_A5[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_A5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_A6-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_A6.sce <- PrepForModeling(EL_A6)
EL_A6.lm <- mixtureModel(EL_A6.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_A6.LMMPlot <- plotModel(EL_A6.sce, EL_A6.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_A6: Posterior Probabilities") |
  plotFiltering(EL_A6.sce, EL_A6.lm) & theme_bw() & ggtitle("EL_A6: Filtering at 75% Posterior Probability Threshold")
EL_A6.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_A6_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_A6.LMMPlot
dev.off()
#Adding to Seurat Object
EL_A6@meta.data <- EL_A6@meta.data %>% mutate(Posterior_Ruptured = 1 - posterior(EL_A6.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_A6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_A6$Likely_Ruptured_075)
(EL_A6[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_A6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_B6-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_B6.sce <- PrepForModeling(EL_B6)
EL_B6.lm <- mixtureModel(EL_B6.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_B6.LMMPlot <- plotModel(EL_B6.sce, EL_B6.lm)  & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_B6: Posterior Probabilities") |
  plotFiltering(EL_B6.sce, EL_B6.lm) & theme_bw() & ggtitle("EL_B6: Filtering at 75% Posterior Probability Threshold")
EL_B6.LMMPlot 
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_B6_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_B6.LMMPlot 
dev.off()
#Adding to Seurat Object
EL_B6@meta.data <- EL_B6@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_B6.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_B6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_B6$Likely_Ruptured_075)
(EL_B6[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_B6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_B7-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_B7.sce <- PrepForModeling(EL_B7)
EL_B7.lm <- mixtureModel(EL_B7.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_B7.LMMPlot <- plotModel(EL_B7.sce, EL_B7.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_B7: Posterior Probabilities") |
  plotFiltering(EL_B7.sce, EL_B7.lm) & theme_bw() & ggtitle("EL_B7: Filtering at 75% Posterior Probability Threshold")
EL_B7.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_B7_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_B7.LMMPlot
dev.off()
#Adding to Seurat Object
EL_B7@meta.data <- EL_B7@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_B7.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_B7[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_B7$Likely_Ruptured_075)
(EL_B7[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_B7[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_B8-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_B8.sce <- PrepForModeling(EL_B8)
EL_B8.lm <- mixtureModel(EL_B8.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_B8.LMMPlot <- plotModel(EL_B8.sce, EL_B8.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_B8: Posterior Probabilities") |
  plotFiltering(EL_B8.sce, EL_B8.lm) & theme_bw() & ggtitle("EL_B8: Filtering at 75% Posterior Probability Threshold")
EL_B8.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_B8_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_B8.LMMPlot
dev.off()
#Adding to Seurat Object
EL_B8@meta.data <- EL_B8@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_B8.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_B8[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_B8$Likely_Ruptured_075)
(EL_B8[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_B8[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_E1-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_E1.sce <- PrepForModeling(EL_E1)
EL_E1.lm <- mixtureModel(EL_E1.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_E1.LMMPlot <- plotModel(EL_E1.sce, EL_E1.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_E1: Posterior Probabilities") |
  plotFiltering(EL_E1.sce, EL_E1.lm) & theme_bw() & ggtitle("EL_E1: Filtering at 75% Posterior Probability Threshold")
EL_E1.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_E1_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_E1.LMMPlot
dev.off()
#Adding to Seurat Object
EL_E1@meta.data <- EL_E1@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_E1.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_E1[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_E1$Likely_Ruptured_075)
(EL_E1[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_E1[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_E3-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_E3.sce <- PrepForModeling(EL_E3)
EL_E3.lm <- mixtureModel(EL_E3.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_E3.LMMPlot <- plotModel(EL_E3.sce, EL_E3.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_E3: Posterior Probabilities") |
  plotFiltering(EL_E3.sce, EL_E3.lm) & theme_bw() & ggtitle("EL_E3: Filtering at 75% Posterior Probability Threshold")
EL_E3.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_E3_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_E3.LMMPlot
dev.off()
#Adding to Seurat Object
EL_E3@meta.data <- EL_E3@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_E3.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_E3[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_E3$Likely_Ruptured_075)
(EL_E3[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_E3[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_E5-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_E5.sce <- PrepForModeling(EL_E5)
EL_E5.lm <- mixtureModel(EL_E5.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_E5.LMMPlot <- plotModel(EL_E5.sce, EL_E5.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_E5: Posterior Probabilities") |
  plotFiltering(EL_E5.sce, EL_E5.lm) & theme_bw() & ggtitle("EL_E5: Filtering at 75% Posterior Probability Threshold")
EL_E5.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_E5_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_E5.LMMPlot
dev.off()
#Adding to Seurat Object
EL_E5@meta.data <- EL_E5@meta.data %>% mutate(Posterior_Ruptured = 1- posterior(EL_E5.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_E5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_E5$Likely_Ruptured_075)
(EL_E5[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_E5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_H3-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_H3.sce <- PrepForModeling(EL_H3)
EL_H3.lm <- mixtureModel(EL_H3.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_H3.LMMPlot <- plotModel(EL_H3.sce, EL_H3.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_H3: Posterior Probabilities") |
  plotFiltering(EL_H3.sce, EL_H3.lm) & theme_bw() & ggtitle("EL_H3: Filtering at 75% Posterior Probability Threshold")
EL_H3.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_H3_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_H3.LMMPlot
dev.off()
#Adding to Seurat Object
EL_H3@meta.data <- EL_H3@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_H3.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_H3[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_H3$Likely_Ruptured_075)
(EL_H3[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_H3[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_H5-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_H5.sce <- PrepForModeling(EL_H5)
EL_H5.lm <- mixtureModel(EL_H5.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_H5.LMMPlot <- plotModel(EL_H5.sce, EL_H5.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_H5: Posterior Probabilities") |
  plotFiltering(EL_H5.sce, EL_H5.lm) & theme_bw() & ggtitle("EL_H5: Filtering at 75% Posterior Probability Threshold")
EL_H5.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_H5_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_H5.LMMPlot
dev.off()
#Adding to Seurat Object
EL_H5@meta.data <- EL_H5@meta.data %>% mutate(Posterior_Ruptured = posterior(EL_H5.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_H5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_H5$Likely_Ruptured_075)
(EL_H5[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_H5[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----EL_H6-----#
#Setting seed 
set.seed(1701)
#Calculating LMM 
EL_H6.sce <- PrepForModeling(EL_H6)
EL_H6.lm <- mixtureModel(EL_H6.sce, model_type = "linear")
#Plotting Model with 75% Posterior Cutoff
EL_H6.LMMPlot <- plotModel(EL_H6.sce, EL_H6.lm) & ggsci::scale_color_gsea() & theme_bw() & ggtitle("EL_H6: Posterior Probabilities") |
  plotFiltering(EL_H6.sce, EL_H6.lm) & theme_bw() & ggtitle("EL_H6: Filtering at 75% Posterior Probability Threshold")
EL_H6.LMMPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/EL_H6_Posterior_ThresholdPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
EL_H6.LMMPlot
dev.off()
#Adding to Seurat Object
EL_H6@meta.data <- EL_H6@meta.data %>% mutate(Posterior_Ruptured = 1- posterior(EL_H6.lm)[,1], Likely_Ruptured_075 = Posterior_Ruptured > 0.75)
#Verifying Correct Assignment Direction 
EL_H6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5)
table(EL_H6$Likely_Ruptured_075)
(EL_H6[[]] %>% ggplot(aes(x = percent.mt, y = "")) + ggridges::geom_density_ridges(alpha = 0.5)) /
  EL_H6[[]] %>% ggplot(aes(x = percent.mt, y = "", fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & geom_vline(xintercept = 12.5)

#-----Merging Data Together------#
#Merging into one Seurat Object
HHA <- merge(EL_A5, list(EL_A6, EL_B6, EL_B7, EL_B8, EL_E1,EL_E3, EL_E5, EL_H3, EL_H5, EL_H6))
HHA
#Removing individual datasets
rm(EL_A5, EL_A6, EL_B6, EL_B7, EL_B8, EL_E1, EL_E3, EL_E5, EL_H3, EL_H5, EL_H6,
   EL_A5.sce, EL_A6.sce, EL_B7.sce, EL_B8.sce, EL_E1.sce, EL_E3.sce, EL_E5.sce, EL_H3.sce, EL_H5.sce, EL_H6.sce)
HHA$DonorID <- as.factor(HHA$DonorID)
HHA$Isolation <- factor(HHA$Isolation, levels = c("4mm", "FUE", "Bulb"))

#----Plotting LMM Results-----#
#Plotting Ridgeplots of mitochondrial read percentage in likely dead and likely viable cells. 
HHA[[]] %>% ggplot(aes(x = percent.mt, y = SampleID, fill = Likely_Ruptured_075)) + ggridges::geom_density_ridges(alpha = 0.5) & 
  geom_vline(xintercept = 10, linetype = "dashed") & theme_bw() & 
  labs(x = "Mitochondrial UMI Percentage", y = "SampleID", fill = "Predicted Rupture Status", 
       title = "Mitochondrial UMI Percentage in Likely Live/Dead Cells by Sample") &
  theme(axis.title = element_text(face = "bold"), plot.title = element_text(face = "bold"))
table(HHA$Likely_Ruptured_075) #113,028 likely viable cells  19,058 dead cells
#Plotting percent mt as violin plots 
CellDeathMitoPlot <- VlnPlot(HHA, features = "percent.mt", group.by = "SampleID", split.by = "Likely_Ruptured_075", alpha = 0.01, split.plot = T, raster = F) +
  labs(x = "", y = "Mitochondrial UMI Percentage", title = "Mitochondrial UMI Percentage Distribution in Likely Live/Dead Cells by Sample", fill = "Predicted Rupture Status")
CellDeathMitoPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/PercentMTbyPrediction_10_10_24.png"), res = 300, height = 6, width = 12, units = "in")
CellDeathMitoPlot
dev.off()
#Plotting posterior probability of rupture as violin plots 
CellDeathPosteriorPlot <- VlnPlot(HHA, features = "Posterior_Ruptured", group.by = "SampleID", split.by = "Likely_Ruptured_075", split.plot = T, alpha = 0, raster = F) +
  labs(x = "", y = "Posterior Probability of Rupture", title = "Posterior Probabilities by Sample", fill = "Predicted Rupture Status")
CellDeathPosteriorPlot 
ragg::agg_tiff(filename = paste0(PlotsDirectory, "LMM_QC/CellDeathPosteriorsPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
CellDeathPosteriorPlot
dev.off()
#Calculating cell death rates by dataset and saving to csv
Sample.DeathRates <- HHA[[]] %>% group_by(SampleID) %>% summarize(Percent_Viable = mean(Posterior_Ruptured < 0.75) * 100)
Sample.DeathRates
write_csv(Sample.DeathRates, file = "/projects/p31989/HHA_Spatial_SC_Atlas/Dataset_Metadata/HHA_Sample_Estimated_CellDeath_Rates_10_1_24.csv")

#-------Quality Control and Filtering-----#
PreQCMetadata <- HHA[[]] 
PreQCMetadata$Filtering = "Pre-Filtering"
#Getting median QC params per dataset
PreQCParams <- PreQCMetadata %>% 
  group_by(orig.ident) %>%
  summarize(across(.cols = c("nCount_RNA", "nFeature_RNA", "percent.mt"), .fn = median, .names = {"median_{.col}"}), 
            NCells = n(),
            nFeature25 = quantile(nFeature_RNA, prob = 0.25), nFeature95 = quantile(nFeature_RNA, prob = 0.95), nMito90 = quantile(percent.mt, prob = 0.90), nCount95 = quantile(nCount_RNA, prob = 0.95))
PreQCParams 
#Computing proportion of cells retained following filtering. 
PreQCCellProps <- PreQCMetadata %>% group_by(orig.ident) %>% summarize(PropRetained = mean(nFeature_RNA > 800 & nFeature_RNA < 8000 & nCount_RNA < 65000 & Posterior_Ruptured < 0.75), CellsRetained = PropRetained * n())
PreQCCellProps
PreQCCellProps %>% pull(CellsRetained) %>% sum() #103,755 cells retained following filtering 

#-----Plotting Initial QC Metrics-----#
#Violin Plots
PreQCVlns <- VlnPlot(HHA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), cols = SampleIDPal, pt.size = 0, ncol = 4, raster = F)
PreQCVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCVlns_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCVlns
dev.off()
#Cell Count by Sample
PreQCSampleCountPlot <- PreQCMetadata %>%
  ggplot(aes(x = orig.ident, fill = orig.ident)) + geom_bar() + theme_bw() + labs(x = "Sample ID", y = "Cell Count", title = "Cell Count by Sample") + theme(
    legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + scale_fill_manual(values = SampleIDPal)
PreQCSampleCountPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCSampleCountPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCSampleCountPlot
dev.off()
#Cell Proportion by Sample 
PreQCSamplePropPlot <- PreQCMetadata %>%
  ggplot(aes(x = "", fill = orig.ident)) + geom_bar(position = "fill") + theme_bw() + labs(x = "", y = "Cell Proportion", title = "Cell Proportions by Sample") + theme(
    plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + scale_fill_manual(values = SampleIDPal)
PreQCSamplePropPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCSamplePropPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCSamplePropPlot
dev.off()
#Cell Count by Donor 
PreQCDonorCountPlot <- PreQCMetadata %>%
  ggplot(aes(x = DonorID, fill = DonorID)) + geom_bar() + theme_bw() + labs(x = "Donor ID", y = "Cell Count", title = "Cell Count by Donor") + theme(
    legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + scale_fill_manual(values = MetBrewer::met.brewer("Signac", 8))
PreQCDonorCountPlot
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCDonorCountPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCSamplePropPlot
dev.off()

#-----Plotting Cutoffs for each metric-----#
#Plotting nCount_RNA Cutoffs
PreQCCountRidges <- RidgePlot(HHA, features = c("nCount_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  geom_vline(xintercept = 65000, linetype = "dashed", color = "black") + 
  labs(x = "", y = "", title = "HHA QC: nCount_RNA Cutoffs") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PreQCCountRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCCountRidgesPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCCountRidges
dev.off()
#Plotting nFeature_RNA Cutoffs
PreQCFeatureRidges <- RidgePlot(HHA, features = c("nFeature_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  geom_vline(xintercept = 800, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 8000, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "HHA QC: nFeature_RNA Cutoffs") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PreQCFeatureRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCFeatureRidgesPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCFeatureRidges
dev.off()
#Plotting percent.mt cutoffs
PreQCMitoRidges <- RidgePlot(HHA, features = c("percent.mt")) + scale_fill_manual(values = SampleIDPal) + 
  labs(x = "", y = "", title = "HHA QC: percent.mt Cutoffs") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PreQCMitoRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCMitoRidgesPlot_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCMitoRidges
dev.off()
#Plotting ribo
PreQCRiboRidges <- RidgePlot(HHA, features = c("percent.ribo")) + scale_fill_manual(values = MetBrewer::met.brewer("Signac", 11)) + 
  labs(x = "", y = "", title = "HHA QC: percent.ribo") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PreQCRiboRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PreFiltering_QC_Metrics/PreQCRiboRidges_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PreQCRiboRidges
dev.off()

#-----Filtering out low quality cells-----#
#Collecting low feature cells annotated as alive 
HHA.Rejects <- subset(HHA, subset = nFeature_RNA < 800 & Likely_Ruptured_075 == F)
HHA.Rejects #5,849 cells with less than 500 genes and annotated as likely viable. 
save(HHA.Rejects, file = paste0(SeuratDirectory, "HHA_Datasets_Lower_Bound_QC_Failures_10_10_24.RData"))
#Filtering to include cells passing quality control thresholds. 
#Performing Quality Control 
HHA <- subset(HHA, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & nCount_RNA < 65000 & Likely_Ruptured_075 == F) #using 75% posterior probability of being dead as cutoff. 
HHA #103,902 cells retained pre-doublet calling
PostQCVlns <- VlnPlot(HHA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 4, raster = F, cols = SampleIDPal)
PostQCVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PostFiltering_QC_Metrics/PostQCVlns_10_10_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostQCVlns
dev.off()

#-----Plotting PostQC RidgePlots-----#
#Plotting nFeature_RNA Cutoffs
PostQCCountRidges <- RidgePlot(HHA, features = c("nCount_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  geom_vline(xintercept = 65000, linetype = "dashed", color = "black") + 
  labs(x = "", y = "", title = "HHA QC: nCount_RNA Post-Filtering") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(0, 65000)
PostQCCountRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PostFiltering_QC_Metrics/PostQCCountRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostQCCountRidges
dev.off()
#Plotting nFeature_RNA Cutoffs
PostQCFeatureRidges <- RidgePlot(HHA, features = c("nFeature_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  geom_vline(xintercept = 800, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 8000, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "HHA QC: nFeature_RNA Post-Filtering") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(800, 8000)
PostQCFeatureRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PostFiltering_QC_Metrics/PostQCFeatureRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostQCFeatureRidges
dev.off()
#Plotting percent.mt cutoffs
PostQCMitoRidges <- RidgePlot(HHA, features = c("percent.mt")) + scale_fill_manual(values = SampleIDPal) + 
  labs(x = "", y = "", title = "HHA QC: percent.mt Post-Filtering") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(0, max(HHA$percent.mt))
PostQCMitoRidges 
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PostFiltering_QC_Metrics/PostQCMitoRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostQCMitoRidges
dev.off()
#Plotting percent.ribo
PostQCRiboRidges <- RidgePlot(HHA, features = c("percent.ribo")) + scale_fill_manual(values = SampleIDPal) + 
  labs(x = "", y = "", title = "HHA QC: percent.ribo Post-Filtering") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(0, max(HHA$percent.ribo))
PostQCRiboRidges 
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "PostFiltering_QC_Metrics/PostQCRiboRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostQCRiboRidges
dev.off()

#-----Saving Post-Filtering QC Metrics-----#
PostQCMetadata <- HHA[[]] %>% mutate(Filtering = "Post-Filtering")

#-----Running scDblFinder-----#
#Splitting into list of Seurat Objects
HHA <-  JoinLayers(HHA) %>% Seurat::SplitObject(split.by = 'orig.ident') #fixes weird "invalid rownames" error that doesn't appear when splitting in other contexts
HHA
#Running scDblFinder on each dataset, using default parameters
HHA <- map(HHA, ~ RunSCDblFinder(.x))
HHA <- merge(HHA[[1]], HHA[2:length(HHA)])
HHA
table(HHA$scDblFinder.class) #94,986 singlets, 8,789 doublets
#Adding factored metadata column with nice capitalization
HHA$DoubletStatus <- factor(str_to_title(HHA$scDblFinder.class), levels = c("Singlet", "Doublet"))

#-----Plotting QC Metrics by scDblFinder Prediction-----#
#Plotting Gene Count by Doublet Status
DoubletGenePlot <- VlnPlot(HHA, features = "nFeature_RNA", group.by = "SampleID", split.by = "DoubletStatus", split.plot = T, pt.size = 0.01, alpha = 0.025, 
                           cols = c("steelblue3", "firebrick"), raster = F) + 
  labs(x = "", y = "nFeature_RNA", title = "Detected Gene Count by Doublet Status", fill = "scDblFinder Prediction") + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold", size = 10), legend.key.size = unit(0.5, "cm"))
DoubletGenePlot
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "scDblFinder_Metrics/DoubletGenePlot_10_4_24.png"), res = 300, height = 6, width = 10, units = "in") 
DoubletGenePlot
dev.off()
#Plotting UMI Count by Doublet Status
DoubletUMIPlot <- VlnPlot(HHA, features = "nCount_RNA", group.by = "SampleID", split.by = "DoubletStatus", split.plot = T, pt.size = 0.01, alpha = 0.025, 
                          cols = c("steelblue3", "firebrick"), raster = F) + 
  labs(x = "", y = "nCount_RNA", title = "UMI Count by Doublet Status", fill = "scDblFinder Prediction")  + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold", size = 10), legend.key.size = unit(0.5, "cm"))
DoubletUMIPlot 
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "scDblFinder_Metrics/DoubletUMIPlot_10_4_24.png"), res = 300, height = 6, width = 10, units = "in") 
DoubletUMIPlot 
dev.off()
#Plotting Mitochondrial UMI Percentage by Doublet Status
DoubletMitoPlot <- VlnPlot(HHA, features = "percent.mt", group.by = "SampleID", split.by = "DoubletStatus", split.plot = T, pt.size = 0.01, alpha = 0.025, 
                           cols = c("steelblue3", "firebrick"), raster = F) + 
  labs(x = "", y = "percent.mt", title = "Mitochondrial UMI Percentage by Doublet Status", fill = "scDblFinder Prediction")  + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold", size = 10), legend.key.size = unit(0.5, "cm"))
DoubletMitoPlot
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "scDblFinder_Metrics/DoubletMitoPlot_10_4_24.png"), res = 300, height = 6, width = 10, units = "in") 
DoubletMitoPlot
dev.off()
#Plotting Ribosomal UMI Percentage by Doublet Status
DoubletRiboPlot <- VlnPlot(HHA, features = "percent.ribo", group.by = "SampleID", split.by = "DoubletStatus", split.plot = T, pt.size = 0.01, alpha = 0.025, 
                           cols = c("steelblue3", "firebrick"), raster = F) + 
  labs(x = "", y = "percent.ribo", title = "Ribosomal UMI Percentage by Doublet Status", fill = "scDblFinder Prediction")  + 
  theme(axis.title = element_text(face = "bold"), legend.title = element_text(face = "bold", size = 10), legend.key.size = unit(0.5, "cm"))
DoubletRiboPlot
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "scDblFinder_Metrics/DoubletRiboPlotPlot_10_4_24.png"), res = 300, height = 6, width = 10, units = "in") 
DoubletRiboPlot
dev.off()

#-----Saving Predicted Doublet/Singlet Ratios as csv-----#
DFMetadata <- HHA[[]]
DoubletsBySample <- DFMetadata %>% group_by(SampleID) %>% 
  summarize(TotalCells = n(), TotalSinglets = sum(DoubletStatus == "Singlet"), TotalDoublets = sum(DoubletStatus == "Doublet"), DoubletPercentage = sum(DoubletStatus == "Doublet") / n() * 100)
DoubletsBySample
write_csv(DoubletsBySample, file = "/projects/p31989/HHA_Spatial_SC_Atlas/Dataset_Metadata/HHA_Sample_Estimated_Doublet_Rates_10_3_24.csv")

#-----Plotting Filtering QC Metrics Comparisons-----#
#Removing doublets from data
HHA <- subset(HHA, DoubletStatus == "Singlet")
HHA #93,425 cells
PostDFMetadata <- HHA[[]] %>% mutate(Filtering = "scDblFinder Filtering")
QCMetadata <- PreQCMetadata %>% bind_rows(PostQCMetadata, PostDFMetadata)
QCMetadata$Filtering <- factor(QCMetadata$Filtering, levels = c("Pre-Filtering", "Post-Filtering", "scDblFinder Filtering"))
#Plotting Cells by Dataset
PostDFCellsBySamplePlot <-  ggplot(QCMetadata, aes(x = SampleID, fill = Filtering)) + geom_bar(position = "dodge") + theme_bw() + labs(x = "Sample", y = "Cell Count", title = "Cell Counts by Sample") + theme(
  plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + 
  scale_fill_manual(values = c("steelblue3", "firebrick", "springgreen4"))
PostDFCellsBySamplePlot
#saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCellsBySamplePlot_10_4_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCellsBySamplePlot
dev.off()
#Features per Sample
PostDFCompGeneVlns <- ggplot(QCMetadata, aes(x = orig.ident, fill = Filtering, y = nFeature_RNA)) + geom_violin(position = position_dodge(0.9)) + theme_bw() + labs(x = "Dataset", y = "Gene Count", title = "Gene Counts by Sample") + theme(
  plot.title = element_text(face = "bold", hjust = 0.5, size = 15), axis.title = element_text(face = "bold", size = 15)) + geom_boxplot(width=0.2, color="black", position = position_dodge(0.9), outliers = F) + 
  scale_fill_manual(values = c("steelblue3", "firebrick", "springgreen4"))
PostDFCompGeneVlns
#saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCompGeneVlns_10_4_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCompGeneVlns
dev.off()
#UMIs Per Sample
PostDFCompUMIVlns <- ggplot(QCMetadata, aes(x = orig.ident, fill = Filtering, y = nCount_RNA)) + geom_violin() + theme_bw() + labs(x = "Dataset", y = "UMI Count", title = "UMI Counts by Sample") + theme(
  plot.title = element_text(face = "bold", hjust = 0.5, size = 15), axis.title = element_text(face = "bold", size = 15)) + geom_boxplot(width=0.15, color="black", position = position_dodge(0.9), outliers = F) + 
  scale_fill_manual(values = c("steelblue3", "firebrick", "springgreen4"))
PostDFCompUMIVlns
#saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCompUMIVlnsPlot_10_4_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCompUMIVlns
dev.off()
#Mitochondrial Read Percentage
PostDFCompMitoVlns <- ggplot(QCMetadata, aes(x = orig.ident, fill = Filtering, y = percent.mt)) + geom_violin(position = position_dodge(0.9)) + theme_bw() + labs(x = "Dataset", y = "Mitochondrial Read Percentage", title = "Mitochondrial Read Percentage by Sample") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + geom_boxplot(width=0.05, color="black", position = position_dodge(0.9), outliers = F) + 
  scale_fill_manual(values = c("steelblue3", "firebrick", "springgreen4"))
PostDFCompMitoVlns
#saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCompMitoVlnsPlot_10_4_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCompMitoVlns
dev.off()
#Percent Ribo
PostDFRiboCompVlns <- ggplot(QCMetadata, aes(x = orig.ident, fill = Filtering, y = percent.ribo)) + geom_violin(position = position_dodge(0.9)) + theme_bw() + labs(x = "Dataset", y = "Ribosomal Read Percentage", title = "Ribosomal Read Percentage by Sample") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold")) + geom_boxplot(width=0.05, color="black", position = position_dodge(0.9), outliers = F) + 
  scale_fill_manual(values = c("steelblue3", "firebrick", "springgreen4"))
PostDFRiboCompVlns
#saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFRiboVlns_10_4_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFRiboCompVlns
dev.off()

#-----Final RidgePlots-----#
#Plotting nFeature_RNA Cutoffs
PostDFCountRidges <- RidgePlot(HHA, features = c("nCount_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  labs(x = "", y = "", title = "HHA Post-QC: nCount_RNA") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(0, 65000)
PostDFCountRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCountRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCountRidges
dev.off()
#Plotting nFeature_RNA Cutoffs
PostDFFeatureRidges <- RidgePlot(HHA, features = c("nFeature_RNA")) + scale_fill_manual(values = SampleIDPal) + 
  geom_vline(xintercept = 800, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = 8000, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "HHA QC: nFeature_RNA Post-QC") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(800, 8000)
PostDFFeatureRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFFeatureRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFFeatureRidges
dev.off()
#Plotting percent.mt cutoffs
PostDFMitoRidges <- RidgePlot(HHA, features = c("percent.mt")) + scale_fill_manual(values = SampleIDPal) + 
  labs(x = "", y = "", title = "HHA Post-QC: percent.mt") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))  + xlim(0, max(HHA$percent.mt))
PostDFMitoRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFMitoRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFMitoRidges
dev.off()
#Plotting percent.ribo
PostDFRiboRidges <- RidgePlot(HHA, features = c("percent.ribo")) + scale_fill_manual(values = MetBrewer::met.brewer("Signac", 12)) + 
  labs(x = "", y = "", title = "HHA Post-QC: percent.ribo") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + xlim(0, max(HHA$percent.ribo))
PostDFRiboRidges
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFRiboRidges_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFRiboRidges
dev.off()

#-----Final Violins-----#
#Plotting nCount_RNA Cutoffs
PostDFCountVlns <- VlnPlot(HHA, features = c("nCount_RNA"), cols = SampleIDPal, alpha = 0.05) + 
  labs(x = "", y = "nCount_RNA", title = "HHA Post-QC: UMI Counts by Dataset") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PostDFCountVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFCountVlns_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFCountVlns
dev.off()
#Plotting nFeature_RNA Cutoffs
PostDFFeatureVlns <- VlnPlot(HHA, features = c("nFeature_RNA"), cols = SampleIDPal, alpha = 0.05) + 
  labs(x = "", y = "nFeature_RNA", title = "HHA Post-QC: Unique Genes by Dataset") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PostDFFeatureVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFFeatureVlns_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFFeatureVlns
dev.off()
#Plotting Mitochondrial Read Percentage
PostDFMitoVlns <- VlnPlot(HHA, features = c("percent.mt"), cols = SampleIDPal, alpha = 0.05) + 
  labs(x = "", y = "percent.mt", title = "HHA Post-QC: Mitochondrial UMI Percentage by Dataset") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PostDFMitoVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFMitoVlns_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFMitoVlns
dev.off()
#Plotting Mitochondrial Read Percentage
PostDFRiboVlns <- VlnPlot(HHA, features = c("percent.ribo"), cols = SampleIDPal, alpha = 0.05) + 
  labs(x = "", y = "percent.ribo", title = "HHA Post-QC: Ribosomal UMI Percentage by Dataset") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
PostDFRiboVlns
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Post_Doublet_Filtering_QC_Metrics/PostDFRiboVlns_10_3_24.png"), res = 300, height = 6, width = 12, units = "in") 
PostDFRiboVlns
dev.off()

#-----Saving Seurat Object-----#
save(HHA, file = paste0(SeuratDirectory, "HHA_PreIntegration_10_10_24.RData"))

#-----Recording session information-----#
sessionInfo()
#writing to text file
writeLines(capture.output(sessionInfo()), "/home/cms2775/HHA_Spatial_scRNAseq_2024/sessionInfo/HHA_LMM_Initial_QC_10_3_24_sessionInfo.txt")



