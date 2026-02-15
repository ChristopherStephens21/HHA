#------Pre-Processing and Analysis of HHA Multiomics Data-----#
#This script contains our workflow for initial processing and integration of our single cell multiomics data. 
#We will generate multimodal Seurat objects, perform QC, run scDblFinder based on the RNA data and then integrate with the scRNA-seq data to get clusters. 
#-----Loading Packages-----#
library(Seurat)
library(scCustomize)
library(scDblFinder)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(presto)
library(AUCell)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Multiome_SeuratPlots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#----Helper Functions----#
#Creates an RNA-level Seurat Object from the Multiomic Data
GetMultiomeRNA <- function(input.file, project.name = "SeuratObj") {
  #Read in 10X object
  cat("Reading in Data \n")
  counts.combined <- scCustomize::Read_CellBender_h5_Mat(file_name = input.file)
  cat("Splitting into Gene Expression and Peak Matrices \n")
  feature.info <- rhdf5::h5read(input.file, name = "/matrix/features") #Reading feature information directly from h5 file 
  feature.df <- data.frame(name = make.unique(feature.info$name), feature_type = feature.info$feature_type) #Annotating each feature with its feature type 
  #Extract RNA and ATAC matrices
  counts.rna <- counts.combined[feature.df$feature_type == "Gene Expression",]
  counts.atac <- counts.combined[feature.df$feature_type == "Peaks",]
  cat("Creating Seurat Object \n")
  seuratobj <- CreateSeuratObject(counts = counts.rna, project = project.name)
  #Extract RNA and ATAC level data
  counts.rna <- counts.combined[feature.df$feature_type == "Gene Expression",]
  return(seuratobj)}

#Parses the metadata file and adds relevant dataset information, converting numerical variables to factors  
AddDatasetMetadata <- function(SeuratObj, Metadata = Dataset.Metadata, Sample = SeuratObj@project.name) {
  Sample.Metadata <- Dataset.Metadata %>% filter(SampleID == Sample)
  SeuratObj[["SampleID"]] <-  Sample.Metadata$SampleID
  SeuratObj[["DonorID"]] <- as.factor(Sample.Metadata$DonorID)
  SeuratObj[["Platform"]] <- factor(Sample.Metadata$Platform, levels = c("V3.1", "V4", "Multiome"))
  SeuratObj[["Isolation"]] <- factor(Sample.Metadata$Isolation, levels = c("4mm", "FUE", "Bulb"))
  SeuratObj[["Sex"]] <- factor(Sample.Metadata$Sex)
  SeuratObj[["Age"]] <- as.factor(Sample.Metadata$Age)
  SeuratObj[["Source"]] <- Sample.Metadata$Source
  return(SeuratObj)}

#Adds mitochondrial and ribosomal read percentages
AddMitoRibo <- function(SeuratObj) {
  SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT-", assay = 'RNA')
  SeuratObj[["percent.ribo"]] <- PercentageFeatureSet(SeuratObj, pattern = "^RP[SL]", assay = 'RNA')
  return(SeuratObj)}

#Reformats DefaultAssay() to tidy format for easy piping/use with purr::map()
DefaultAssayTidy <- function(SeuratObj, assay) {
  DefaultAssay(SeuratObj) <- assay
  return(SeuratObj)}

#Generates an SCE object from a Seurat count matrix and runs scDblFinder using the default parameters (clustering workflow)
RunSCDblFinder <- function(SeuratObj) {
  set.seed(1701)
  #Cluster specific information is not really worth keeping as it is not cell-type annotated and varies across datasets. 
  KeptCols <- c("scDblFinder.class", "scDblFinder.score", "scDblFinder.weighted", "scDblFinder.difficulty", "scDblFinder.cxds_score")
  #Running scDblFinder using default parameters
  SDF <-scDblFinder(GetAssayData(SeuratObj, layer ="counts", assay = "RNA"), clusters = T)
  SeuratObj@meta.data <- cbind(SeuratObj@meta.data, as.data.frame(colData(SDF)[, KeptCols])) #Adding output to Seurat Metadata
  return(SeuratObj)}

#Extracts raw count matrix and runs AUCell
RunAUCell <- function(SeuratObj, GeneLists) {
  #----Calculating AUCell Scores-----#
  #Pulling out raw count data
  cat("Pulling Count Matrix \n")
  SeuratObj.counts <- GetAssayData(object = SeuratObj, slot = "counts")
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
#-----Themes-----#
Metadata.Theme <- theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15), axis.title = element_text(face = "bold", size = 14), aspect.ratio = 1.4)

#-----Dataset Metadata-----#
#Relevant metadata for each (Multiome) sample included in atlas
Dataset.Metadata <- read_csv("/projects/p31989/HHA_Spatial_SC_Atlas/Dataset_Metadata/HHA_Multiome_Metadata_9_23_24.csv")

#-----Loading in Data-----#
#EL_C5
EL_C5 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/Multiome/EL_C5_Multi/cellbender_corr_matrix_filtered.h5"
#EL_C6
EL_C6 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/Multiome/EL_C6_Multi/cellbender_corr_matrix_filtered.h5"
#EL_E6
EL_E6 <-  "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/Multiome/EL_E6_Multi/cellbender_corr_matrix_filtered.h5"
#EL_E7
EL_E7 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/Multiome/EL_E7_Multi/cellbender_corr_matrix_filtered.h5"
#EL_E8
EL_E8 <- "/projects/p31989/HHA_Spatial_SC_Atlas/HHA_CellBender_Corrected_Datasets_10_1_24/Multiome/EL_E8_Multi/cellbender_corr_matrix_filtered.h5"

#-----Creating Seurat Objects and adding relevant metadata-----#
EL_C5 <- GetMultiomeRNA(input.file = EL_C5, project.name = "EL_C5") %>% AddDatasetMetadata() %>% AddMitoRibo()
EL_C6 <- GetMultiomeRNA(input.file = EL_C6, project.name = "EL_C6") %>% AddDatasetMetadata() %>% AddMitoRibo()
EL_E6 <- GetMultiomeRNA(input.file = EL_E6, project.name = "EL_E6") %>% AddDatasetMetadata() %>% AddMitoRibo()
EL_E7 <- GetMultiomeRNA(input.file = EL_E7, project.name = "EL_E7") %>% AddDatasetMetadata() %>% AddMitoRibo()
EL_E8 <- GetMultiomeRNA(input.file = EL_E8, project.name = "EL_E8") %>% AddDatasetMetadata() %>% AddMitoRibo()
#Merging objects
HHA.Multiome <- merge(EL_C5, list(EL_C6, EL_E6, EL_E7, EL_E8))
rm(EL_C5, EL_C6, EL_E6, EL_E7, EL_E8)

#-----Performing Quality Control By Individual Datasets-----#
#RNA Metrics
PreQCVlns <- VlnPlot(HHA.Multiome, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), cols = SampleIDPal, ncol = 4, alpha = 0)
PreQCVlns

#-----Merging Metadata Across Samples and Calculating Summary Statistics-----#
PreQCMetadata <- HHA.Multiome[[]]
PreQCMetadata %>% group_by(orig.ident) %>%
  summarize(across(.cols = c("nCount_RNA", "nFeature_RNA", "percent.mt"), .fn = median, .names = {"median_{.col}"}), 
            NCells = n(),
            nFeature20 = quantile(nFeature_RNA, prob = 0.20), nFeature95 = quantile(nFeature_RNA, prob = 0.95),
            nMito90 = quantile(percent.mt, prob = 0.90), nCount95 = quantile(nCount_RNA, prob = 0.95))

#-----Plotting RNA Summary Statistics-----#
#nCount_RNA
ggplot(PreQCMetadata, aes(x = orig.ident, fill = orig.ident, y = nCount_RNA)) + geom_violin() + 
  labs(title = "RNA-Level Data: nCount_RNA", x = "Sample", y = "nCount_RNA") + geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_hline(yintercept = 15000, linetype = "dashed") + Metadata.Theme + Metadata.Pal
#nFeature_RNA
ggplot(PreQCMetadata, aes(x = orig.ident, fill = orig.ident, y = nFeature_RNA)) + geom_violin() + 
  labs(title = "RNA-Level Data: nFeature_RNA", x = "Sample", y = "nFeature_RNA") + geom_hline(yintercept = 500, linetype = "dashed") +
  geom_hline(yintercept = 5000, linetype = "dashed") + Metadata.Theme + Metadata.Pal
#Percent MT
ggplot(PreQCMetadata, aes(x = orig.ident, fill = orig.ident, y = percent.mt)) + geom_violin() + 
  labs(title = "RNA-Level Data: percent.mt", x = "Sample", y = "percent.mt") + geom_hline(yintercept = 25, linetype = "dashed") +
  Metadata.Theme + Metadata.Pal
#Percent Ribo
ggplot(PreQCMetadata, aes(x = orig.ident, fill = orig.ident, y = percent.ribo)) + geom_violin() + 
  labs(title = "RNA-Level Data: percent.ribo", x = "Sample", y = "percent.ribo")  + Metadata.Theme + Metadata.Pal
#Plotting Counts vs Features
ggplot(PreQCMetadata, aes(x = nFeature_RNA, y = percent.mt, color = log(nCount_RNA))) + geom_point(alpha = 0.5, size = 0.25) + geom_hline(yintercept = 25, linetype = "dashed") +
  geom_vline(xintercept = 500, linetype = "dashed") +
  scale_color_viridis_c(option = "inferno") + facet_wrap(~orig.ident) & theme_browser()

#Computing proportion of cells retained following filtering from RNA Metrics alone. 
PreQCCellProps <- PreQCMetadata %>% group_by(orig.ident) %>% 
  summarize(PropRetained = mean(nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 15000 & percent.mt < 25),
            CellsRetained = PropRetained * n())
PreQCCellProps
PreQCCellProps %>% pull(CellsRetained) %>% sum() #40,306 cells pass RNA level filtering

#-----Filtering for RNA level Thresholds-----#
HHA.Multiome <- subset(HHA.Multiome, nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 15000 & percent.mt < 25)
HHA.Multiome

#-----Extracting PostQC Metadata-----#
PostQCMetadata <- HHA.Multiome[[]]
PostQCMetadata
#-----Plotting PostQC RNA Summary Statistics-----#
#nCount_RNA
ggplot(PostQCMetadata, aes(x = orig.ident, fill = orig.ident, y = nCount_RNA)) + geom_violin() + 
  labs(title = "Post-Filtering: nCount_RNA", x = "Sample", y = "nCount_RNA") + geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_hline(yintercept = 15000, linetype = "dashed") + Metadata.Theme + Metadata.Pal
#nFeature_RNA
ggplot(PostQCMetadata, aes(x = orig.ident, fill = orig.ident, y = nFeature_RNA)) + geom_violin() + 
  labs(title = "Post-Filtering: nFeature_RNA", x = "Sample", y = "nFeature_RNA") + geom_hline(yintercept = 500, linetype = "dashed") +
  geom_hline(yintercept = 5000, linetype = "dashed") + Metadata.Theme + Metadata.Pal
#Percent MT
ggplot(PostQCMetadata, aes(x = orig.ident, fill = orig.ident, y = percent.mt)) + geom_violin() + 
  labs(title = "Post-Filtering: percent.mt", x = "Sample", y = "percent.mt") + geom_hline(yintercept = 25, linetype = "dashed") +
  Metadata.Theme + Metadata.Pal
#Percent Ribo
ggplot(PostQCMetadata, aes(x = orig.ident, fill = orig.ident, y = percent.ribo)) + geom_violin() + 
  labs(title = "Post-Filtering: percent.ribo", x = "Sample", y = "percent.ribo")  + Metadata.Theme + Metadata.Pal

#-----Running scDblFinder on Multiome Data-----#
HHA.Multiome <- SplitObject(HHA.Multiome, split.by = "orig.ident")
HHA.Multiome <- HHA.Multiome%>% map(~RunSCDblFinder(.x))
#Function for adding scDblFinder output to SeuratObject
AddDFMetadata <- function(SeuratObj) {SeuratObj$DoubletStatus <- factor(str_to_title(SeuratObj$scDblFinder.class), levels = c("Singlet", "Doublet"))
return(SeuratObj)}
#Adding metadata and merging 
HHA.Multiome <- HHA.Multiome %>% map(~AddDFMetadata(.x))
HHA.Multiome <- merge(HHA.Multiome[[1]],HHA.Multiome[2:length(HHA.Multiome)])
HHA.Multiome$DoubletStatus <- factor(HHA.Multiome$DoubletStatus, levels = c("Singlet", "Doublet"))
table(HHA.Multiome$DoubletStatus) #35,285 singlets   5025 Doublets

#save.image(paste0(WorkspaceDirectory, "Multiome_RNA_Level_Workspace_11_15_24.Rdata"))

#-----Profiling scDblFinder Output-----#
PostDFMetadata <- HHA.Multiome[[]]
#nCount_RNA
ggplot(PostDFMetadata, aes(x = orig.ident, fill = DoubletStatus, y = nCount_RNA)) + geom_violin() + 
  labs(title = "Post- scDblFinder Filtering: nCount_RNA", x = "Sample", y = "nCount_RNA", fill = "scDblFinder Prediction") & 
  scale_fill_manual(values = c("steelblue3", "firebrick")) & Metadata.Theme
#nFeature_RNA
ggplot(PostDFMetadata, aes(x = orig.ident, fill = DoubletStatus, y = nFeature_RNA)) + geom_violin() + 
  labs(title = "Post- scDblFinder Filtering: nFeature_RNA", x = "Sample", y = "nFeature_RNA", fill = "scDblFinder Prediction") & 
  scale_fill_manual(values = c("steelblue3", "firebrick")) & Metadata.Theme
#percent mt
ggplot(PostDFMetadata, aes(x = orig.ident, fill = DoubletStatus, y = percent.mt)) + geom_violin() + 
  labs(title = "Post- scDblFinder Filtering: percent.mt", x = "Sample", y = "percent.mt", fill = "scDblFinder Prediction") & 
  scale_fill_manual(values = c("steelblue3", "firebrick")) & Metadata.Theme
#percent ribo
ggplot(PostDFMetadata, aes(x = orig.ident, fill = DoubletStatus, y = percent.ribo)) + geom_violin() + 
  labs(title = "Post- scDblFinder Filtering: percent.ribo", x = "Sample", y = "percent.ribo", fill = "scDblFinder Prediction") & 
  scale_fill_manual(values = c("steelblue3", "firebrick")) & Metadata.Theme

#----Writing Barcodes to Metadata for ArchR analysis-----#
#Reading in csv
FilteredBarcodes <- HHA.Multiome[[]] %>% mutate(OriginalBarcode = str_remove_all(rownames(HHA.Multiome[[]]), pattern = "_[0-9]")) %>%
  select(OriginalBarcode, SampleID)
FilteredBarcodes
write_csv(FilteredBarcodes, paste0(MetadataDirectory, "HHA_Multiome_PreDF_FilteredBarcodes_11_14_24.csv"))

#-----Reading in ATAC-Level Metadata-----#
PostQCATACMetadata <- read_csv(paste0(MetadataDirectory, "HHA_Multiome_PostfilteringATACMetadata_11_14_24.csv"))
PostQCATACMetadata$UnifiedBarcode
#Creating a column to easily unify barcodes
HHA.Multiome$UnifiedBarcode <- paste(HHA.Multiome$SampleID, str_remove_all(rownames(HHA.Multiome[[]]), pattern = "_[0-9]"), sep = "_")
#Creating a Unified Metadata Object
UnifiedMetadata <- HHA.Multiome[[]] %>% inner_join(PostQCATACMetadata, by = "UnifiedBarcode")

#-----Selecting for Singlets and cells passing ATAC QC-----#
HHA.Multiome <- subset(HHA.Multiome, DoubletStatus == "Singlet" & UnifiedBarcode %in% UnifiedMetadata$UnifiedBarcode)
HHA.Multiome #33,199 singlets passing QC thresholds

#-----Exporting Final Barcodes-----#
PostDFFilteringMetadata <-HHA.Multiome[[]]
PostDFFilteringMetadata
write_csv(PostDFFilteringMetadata, paste0(MetadataDirectory, "HHA_Multiome_PostDoubletFilteringMetadata_11_14_24.csv"))

#-------Scoring Dissociation Signature using AUCell-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation 
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Calculating stress response score with AUCell
Multiome.StressAUC <- HHA.Multiome %>% JoinLayers() %>% RunAUCell(list(DissocResponse = StressSig$Gene_Human))
#Adding to Seurat Object
HHA.Multiome$DissociationScore <- Multiome.StressAUC$DissocResponse
#Visualizing by Dataset with violin plot
VlnPlot(HHA.Multiome, features = "DissociationScore", group.by = "SampleID")

#-----Performing Log Normalization and Cell Cycle Scoring-----#
DefaultAssay(HHA.Multiome) <- "RNA"
HHA.Multiome <- NormalizeData(HHA.Multiome, normalization.method = "LogNormalize", scale.factor = 10000)
#Performing CC scoring on individual datasets
HHA.Multiome <- SplitObject(HHA.Multiome, split.by = "orig.ident")
HHA.Multiome <- map(HHA.Multiome, ~CellCycleScoring(.x, s.features= cc.genes.updated.2019$s.genes, g2m.features= cc.genes.updated.2019$g2m.genes,
                                                    set.ident=F))
HHA.Multiome <- merge(HHA.Multiome[[1]],HHA.Multiome[2:length(HHA.Multiome)])

#-----Unintegrated Clustering/DR-----#
DefaultAssay(HHA.Multiome) <- "RNA"
#Selecting 2000 variable features for dimensionality reduction 
HHA.Multiome <- FindVariableFeatures(HHA.Multiome, selection.method = "vst", nfeatures = 2000)
HHA.Multiome.TopFeaturePlot <- VariableFeaturePlot(HHA.Multiome) %>% LabelPoints(points = head(VariableFeatures(HHA.Multiome), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Multiome.TopFeaturePlot
HHA.Multiome <- ScaleData(HHA.Multiome)
HHA.Multiome <- RunPCA(HHA.Multiome)
DimHeatmap(HHA.Multiome, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Multiome <-FindNeighbors(HHA.Multiome, dims = 1:50, reduction = "pca")
HHA.Multiome <- RunUMAP(HHA.Multiome, dims = 1:50, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 1, cluster.name = "unintegrated_clusters")

#--------Plotting Unintegrated Results-----#
#Plotting unintegrated UMAP
UnintegratedUMAP <- DimPlot(HHA.Multiome, reduction = "UMAP_Unintegrated", raster = F) + 
  labs(title = "Unintegrated Datasets: Louvain Clusters", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15), plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAP
#Plotting unintegrated UMAP
UnintegratedUMAPBySample <-DimPlot(HHA.Multiome, reduction = "UMAP_Unintegrated", group.by = "orig.ident", raster = F, cols = SampleIDPal, shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAPBySample
#Plotting by Isolation Method 
UnintegratedUMAPByIsolation <- DimPlot(HHA.Multiome, reduction = "UMAP_Unintegrated", group.by = "Isolation", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen4")) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAPByIsolation
#Plotting QC Metrics
QCFeaturePlot <- FeaturePlot(HHA.Multiome, reduction = "UMAP_Unintegrated", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "DissociationScore")) & 
  scale_color_viridis_c(option = "inferno")
QCFeaturePlot

#----Preparing for Integration-----#
#-------Performing Scanpy Conversion for SCVI Integration-----#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#Subsetting Atlas for Variable Genes 
Diet.Multiome <- HHA.Multiome[VariableFeatures(HHA.Multiome)]

#----Writing Metadata----#
Multiome.metadata <- Diet.Multiome[[]]
Multiome.metadata$barcode <- rownames(Multiome.metadata)
write.csv(Multiome.metadata, file = paste0(AnnDataDirectory, "Multiome_Metadata_RNA_PreIntegration_3_2_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting count matrix
DefaultAssay(Diet.Multiome) <- "RNA"
DietMultiome.counts <- Diet.Multiome %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
#Writing to mtx file 
Matrix::writeMM(DietMultiome.counts , file = paste0(AnnDataDirectory, "Multiome_RNA_RawCounts_3_2_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietMultiome.counts)), file = paste0(AnnDataDirectory, "Multiome_RNA_GeneNames_3_2_25.csv"),
            quote = F, row.names = F, col.names = T)

#------Loading in SCVI Integration Output------#
Multi.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_Multiome_RNA_LatentRep_3_2_25.csv")) %>% column_to_rownames("Barcode")
head(Multi.scvi)
# put it back in our original Seurat object
HHA.Multiome[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Multi.scvi), key = "scvi_")
rm(Multi.scvi)

#-----Performing Final DR and Clustering-----#
HHA.Multiome <-FindNeighbors(HHA.Multiome, dims = 1:50, reduction = "scvi")
HHA.Multiome <- RunUMAP(HHA.Multiome, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Multiome <- FindClusters(HHA.Multiome, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#-----Plotting Clusters-----#
#Plotting Clusters
DimPlot(HHA.Multiome, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", raster = F, pt.size = 0.01) + 
  NoAxes() + theme(aspect.ratio = 1)
#Plotting by Sample ID
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, group.by = "SampleID", raster = F,
        cols = SampleIDPal, pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
#Plotting by Donor ID
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, group.by = "DonorID", raster = F,
        pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, split.by = "DonorID", raster = F, group.by = "scvi_clusters_0.75",
        pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
#Plotting by Sex
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, group.by = "Sex", raster = F,
        pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, split.by = "Sex", raster = F, group.by = "scvi_clusters_0.75",
        pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
#Plotting by Isolation Method
HHA.Multiome$Isolation <- factor(HHA.Multiome$Isolation, levels = c("4mm", "FUE", "Bulb"))
DimPlot(HHA.Multiome, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", split.by = "Isolation", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()
DimPlot(HHA.Multiome, reduction = "UMAP", label = F, shuffle = T, group.by = "Isolation", raster = F, pt.size = 0.1) + 
  theme(aspect.ratio = 1) + NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Multiome, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "DissociationScore")) & 
  scale_color_viridis_c(option = "inferno")

#------Plotting Marker Genes-------#
#Major Groups
FeaturePlot(HHA.Multiome, features = c("KRT5", "VIM", "COL1A1", "PECAM1", "PTPRC", "RGS5", "DCT", "MPZ", 'KRT7'), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Matrix
FeaturePlot(HHA.Multiome, features = c("KRT35", "KRT85", "GATA3", "FABP4", "ALDH3A1", "FOXQ1", "KRT23", "KRT75", "TCHH"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#ORS
FeaturePlot(HHA.Multiome, features = c("KRT15", "KRT17", "KRT6A", "DIO2"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Cuticle and Late Cortex 
FeaturePlot(HHA.Multiome, features = c("SOX21", "KRT32", "KRT82", "KRT36"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#IFE and UpHF
FeaturePlot(HHA.Multiome, features = c("S100A8", "IGFBP3", "KRT1", "SOX9"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#T Cells
FeaturePlot(HHA.Multiome, features = c("PTPRC", "CD74", "CD3G", "CD68", "CD8A", "NKG7"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Immune Cells
FeaturePlot(HHA.Multiome, features = c("CD207", "MZB1", "CD79A", "MS4A1"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Sebaceous Gland
FeaturePlot(HHA.Multiome, features = c("AWAT2", "SCD", "PPARG", "MGST1"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Proliferating
FeaturePlot(HHA.Multiome, features = c("MKI67"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)

#-----Finding Markers for Each Cluster-----#
#Getting markers for each cluster 
MultiomeMarkers.TFIDF <- HHA.Multiome %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Multiome$scvi_clusters_0.75, N = 20, FDR = 0.01)
MultiomeMarkers.TFIDF
MultiomeMarkers.TFIDFList <- MultiomeMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(MultiomeMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Multiome$scvi_clusters_0.75))-1))
MultiomeMarkers.TFIDFList


#------Annotating Clusters-----#
MultiomeClusters.RNA <- c("Cortex I", #0
                          "Suprabasal ORS", #1
                          "Fibroblasts", #2
                          "Cortex II", #3
                          "Bulge", #4
                          "Basal Isthmus/Infundibulum", #5
                          "Inner Root Sheath", #6
                          "Late Cortex", #7
                          "Suprabasal IFE", #8
                          "Endothelial", #9
                          "Proximal ORS", #10
                          "Basal ORS", #11
                          "Smooth Muscle/Dermal Sheath", #12
                          "Cuticle", #13
                          "Antigen Presenting Cells I", #14
                          "Basal IFE", #15
                          "T Cells I", #16
                          "Suprabasal Isthmus", #17
                          "Upper Companion", #18
                          "Melanocytes I", #19
                          "B Cells", #20
                          "Melanocytes II", #21
                          "Proximal Companion", #22
                          "Eccrine Epithelial", #23
                          "Eccrine Myoepithelial", #24
                          "Suprabasal Upper HF", #25
                          "T Cells II", #26
                          "Lymph Vasculature", #27
                          "Schwann", #28
                          "Antigen Presenting Cells II", #29
                          "Dermal Papilla") #30


#Renaming Idents
Idents(HHA.Multiome) <- "scvi_clusters_0.75"
names(MultiomeClusters.RNA) <- levels(HHA.Multiome)
HHA.Multiome <- RenameIdents(HHA.Multiome, MultiomeClusters.RNA)
HHA.Multiome$MultiomeAnnotationRNA <- HHA.Multiome@active.ident
DimPlot(HHA.Multiome, reduction = "UMAP", label = T, repel = T)  & NoLegend() & NoAxes() 
#Factoring Idents
HHA.Multiome$MultiomeAnnotationOrdered <- factor(HHA.Multiome$MultiomeAnnotationRNA, levels = c("Cortex I", "Cortex II", "Late Cortex", "Cuticle", "Inner Root Sheath", "Proximal Companion", "Proximal ORS", "Basal ORS", "Bulge",
                                                                                                "Suprabasal ORS", "Suprabasal Isthmus", "Basal Isthmus/Infundibulum", "Suprabasal Upper HF", "Basal IFE", "Suprabasal IFE", "Upper Companion",
                                                                                                "Eccrine Epithelial", "Eccrine Myoepithelial", "Smooth Muscle/Dermal Sheath", 'Endothelial', "Lymph Vasculature", "Fibroblasts", "Dermal Papilla",
                                                                                                "Melanocytes I", "Melanocytes II", "Schwann", "Antigen Presenting Cells I", "Antigen Presenting Cells II", "T Cells I", "T Cells II", "B Cells"))
DimPlot(HHA.Multiome, reduction = "UMAP", group.by = "MultiomeAnnotationOrdered", label = T, repel = T)  & NoLegend() & NoAxes() & theme(aspect.ratio = 1)


DotPlot(HHA.Multiome, group.by = "MultiomeAnnotationOrdered", features = c("MSX2", "KRT35", "KRT85", "KRT31", "KRT32", "KRT25", "GATA3", "KRT75", "TBX1", 'KRT5', "KRT15", "DIO2",  "LGR5", "KRT17", "CTNND2", "BARX2", "SOX9", "KRT6A", "S100A8", "IGFBP3", "KRT1", 'KRT20', "KRT18", "KRT7",
                                                                           "ACTA2", "TAGLN", "CHRM3", "RGS5", "ITGA8", "PECAM1", "CDH5", "PROX1", "CCL21", "PDGFRA", "DCN", "RSPO3", "DCT", "PMEL", "CDH19", "MPZ",
                                                                           "NRXN1", "PTPRC", "FCER1G", "CD14", 'CD68', "CD207", "TPSB2", "CD3D", "CD3G", "NKG7", "MS4A1", "CD79A", "MZB1", "IGHG1", 'VIM')) + RotatedAxis() + 
  scale_color_gradientn(colors = c("black","#2488F0","#7F3F98","#E22929", "#FCB31A")) + 
  theme(aspect.ratio = 1, panel.background = element_rect(color = "black"))
table(HHA.Multiome$MultiomeAnnotationRNAOrdered)

#----Transferring Cluster Information to csv for ArchR-----#
write_csv(HHA.Multiome[[]], paste0(MetadataDirectory, "HHA_Multiome_MultiomeAnnotatedMetadata_11_15_24.csv"))

#------Saving Object and Workspace-----#
save(HHA.Multiome, file = paste0(SeuratDirectory, "HHA_Multiome_Integrated_Annotated_3_2_25.rds"))
save.image(paste0(WorkspaceDirectory, "HHA_Multiome_Integration_Workspace_3_2_25.Rdata"))

#-----Recording session information-----#
sessionInfo()
#writing to text file
writeLines(capture.output(sessionInfo()), "/home/cms2775/HHA_Spatial_scRNAseq_2024/sessionInfo/HHA_Multiome_RNA_Level_Clustering_3_2_25_sessionInfo.txt")
