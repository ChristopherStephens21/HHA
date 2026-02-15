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
library(clusterProfiler)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Bulb_Seurat_Plots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#----Helper Functions----#
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

#---------Loading Checkpoint-------#
load(paste0(WorkspaceDirectory, "HHA_Bulb_Recluster_5_22_25.RData"))

#---------Loading Full Object-------#
load(paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))
HHA.Full

#GlobalAnnotation
DimPlot(HHA.Full, reduction = "UMAP", group.by = "GlobalAnnotation", label = T, repel = T, raster = F) +
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()
#Broad Annotation
DimPlot(HHA.Full, reduction = "UMAP", group.by = "BroadAnnotation", label = T, repel = T, raster = F) +
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()

#-------Subsetting Out Bulb Cells-----#
HHA.Bulb <- subset(HHA.Full, BroadAnnotation == "Matrix_Derived")
HHA.Bulb #24,028 cells, 12344 Multiome cells
#-----Finding Variable Features-----#
#Selecting 2000 variable features for dimensionality reduction
HHA.Bulb <- FindVariableFeatures(HHA.Bulb, selection.method = "vst", nfeatures = 2000)
HHA.Bulb.TopFeaturePlot <- VariableFeaturePlot(HHA.Bulb) %>% LabelPoints(points = head(VariableFeatures(HHA.Bulb), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Bulb.TopFeaturePlot
#Scaling Data
HHA.Bulb <- ScaleData(HHA.Bulb)
#Running PCA
HHA.Bulb <- RunPCA(HHA.Bulb)
DimHeatmap(HHA.Bulb, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Bulb <-FindNeighbors(HHA.Bulb, dims = 1:30, reduction = "pca")
HHA.Bulb <- RunUMAP(HHA.Bulb, dims = 1:30, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Bulb <- FindClusters(HHA.Bulb, resolution = 1, cluster.name = "unintegrated_clusters")

#------Plotting UMAP-----#
#Plotting unintegrated UMAP
DimPlot(HHA.Bulb, reduction = "UMAP_Unintegrated", raster = F) +
  labs(title = "Unintegrated Datasets: Louvain Clusters", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15), plot.title = element_text(size = 20, hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Sample
DimPlot(HHA.Bulb, reduction = "UMAP_Unintegrated", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T) +
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Donor
DimPlot(HHA.Bulb, reduction = "UMAP_Unintegrated", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T) +
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting By Platform
HHA.Bulb$Platform <- factor(HHA.Bulb$Platform, levels = c("V3.1", "V4", "Multiome"))
DimPlot(HHA.Bulb, reduction = "UMAP_Unintegrated", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) +
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting By Annotation
DimPlot(HHA.Bulb, reduction = "UMAP_Unintegrated", group.by = "GlobalAnnotation", raster = F, shuffle = T) +
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP_Unintegrated", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()
#Plotting Dissociation Score
FeaturePlot(HHA.Bulb, reduction = "UMAP_Unintegrated", features = c("DissociationScore")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()

#------Removing Cell Cycle genes from Variable Features------#
#-----Removing Cell Cycle Genes and Stress Genes from HVGs-----#
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
VariableFeatures(HHA.Bulb) <- VariableFeatures(HHA.Bulb)[!VariableFeatures(HHA.Bulb) %in% CellCycleGenes]
VariableFeatures(HHA.Bulb) %>% length() #1955 genes
#------Pulling Cell Cycle Associated Genes-----#
#Tirosh et al Cell Cycle Genes
CanonicalCellCycle <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Subsetting Seurat Object for genes in stress signature
Bulb.CellCycle <- HHA.Bulb[CanonicalCellCycle,]

#-------Performing Scanpy Conversion for SCVI Integration: No Dissociation Signature Gene Removal-----#
AnnDataDirectory <- "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/"
#Subsetting Atlas for Variable Genes
Diet.Bulb <- HHA.Bulb[VariableFeatures(HHA.Bulb)]
#----Writing Metadata----#
Bulb.metadata <- Diet.Bulb[[]]
Bulb.metadata$barcode <- rownames(Bulb.metadata)
write.csv(Bulb.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_Metadata_5_22_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.Bulb) <- "RNA"
#Writing log normalized counts to mtx
DietBulb.counts <- Diet.Bulb %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietBulb.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietBulb.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.Bulb, DietBulb.counts, Bulb.metadata)
#------Writing Normalized Cell Cycle Gene Expression------#
#Writing as csv
CellCycle.GeneMat <- Bulb.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_CellCycle_LogNorm_5_22_25.csv"))

#-----Getting SCVI Integration Results-----#
Bulb.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Bulb_LatentRep_5_22_25.csv")) %>% column_to_rownames("Barcode")
head(Bulb.scvi)
# put it back in our original Seurat object
HHA.Bulb[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Bulb.scvi), key = "scvi_")
rm(Bulb.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Bulb <-FindNeighbors(HHA.Bulb, dims = 1:50, reduction = "scvi")
HHA.Bulb <- RunUMAP(HHA.Bulb, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Bulb, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 0.5) 
#-----Plotting by Sample-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T)
#-----Plotting by Donor-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T, pt.size = 0.5) 
#------Plotting By Platform-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) 
#-----Plotting By Annotation-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T)
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("scDblFinder.score"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-------Plotting Markers--------#
#Inner layer and ORS
FeaturePlot(HHA.Bulb, features = c("KRT5", "COL17A1", "MKI67", "LGR5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Inner Root Sheath
FeaturePlot(HHA.Bulb, features = c("GATA3", "TCHH", "KRT72", "KRT73"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cortex
FeaturePlot(HHA.Bulb, features = c("ODC1", "KRT5", "KRT35", "KRT85", "LEF1", "KRT31"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cuticle
FeaturePlot(HHA.Bulb, features = c("KRT32", "KRT82", "CRAT", "SOX21"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proximal Companion and Medulla
FeaturePlot(HHA.Bulb, features = c("KRT75", "TCHH", "FOXQ1", "ALDH3A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#ORS and Companion Beyond Matrix
FeaturePlot(HHA.Bulb, features = c("KRT75", "FGF5", "KRT6A", "BARX2"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proliferation
FeaturePlot(HHA.Bulb, features = c("MKI67"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(HHA.Bulb, features = c("COL17A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Problem Genes
FeaturePlot(HHA.Bulb, features = c("FOS", "JUN", "KRT35", "KRT5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Violin Plot of Problematic Genes and Dissociation Score
VlnPlot(HHA.Bulb, features = c("FOS", "JUN", "DissociationScore"), alpha = 0.1)

#-----Removing Dissociation-Heavy Cluster-----#
HHA.Bulb <- subset(HHA.Bulb, scvi_clusters_0.4 != 0)
HHA.Bulb #16,749 cells
#-----Plotting Clusters-----#
DimPlot(HHA.Bulb, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 0.5) 
table(HHA.Bulb$Platform) #Multiome cells are mostly ok

#-----Reclustering after cleaning -----#
#Selecting 2000 variable features for dimensionality reduction
HHA.Bulb <- FindVariableFeatures(HHA.Bulb, selection.method = "vst", nfeatures = 2000)
HHA.Bulb.TopFeaturePlot <- VariableFeaturePlot(HHA.Bulb) %>% LabelPoints(points = head(VariableFeatures(HHA.Bulb), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Bulb.TopFeaturePlot
#Scaling Data
HHA.Bulb <- ScaleData(HHA.Bulb)
#Running PCA
HHA.Bulb <- RunPCA(HHA.Bulb)
DimHeatmap(HHA.Bulb, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Bulb <-FindNeighbors(HHA.Bulb, dims = 1:30, reduction = "pca")
HHA.Bulb <- RunUMAP(HHA.Bulb, dims = 1:30, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Bulb <- FindClusters(HHA.Bulb, resolution = 1, cluster.name = "unintegrated_clusters")

#------Getting Stress Signature Associated Genes-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Subsetting Seurat Object for genes in stress signature
Bulb.StressSig <- HHA.Bulb[StressSig$Gene_Human,]
#------Pulling Cell Cycle Associated Genes-----#
#Tirosh et al Cell Cycle Genes
CanonicalCellCycle <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Subsetting Seurat Object for genes in stress signature
Bulb.CellCycle <- HHA.Bulb[CanonicalCellCycle,]

#------Removing Dissociation Signature and Cell Cycle genes from Variable Features------#
#-----Removing Cell Cycle Genes and Stress Genes from HVGs-----#
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
VariableFeatures(HHA.Bulb) <- VariableFeatures(HHA.Bulb)[!VariableFeatures(HHA.Bulb) %in% CellCycleGenes]
VariableFeatures(HHA.Bulb) %>% length() #1959 genes
#Removing Stress Signature Genes
VariableFeatures(HHA.Bulb) <- VariableFeatures(HHA.Bulb)[!VariableFeatures(HHA.Bulb) %in% StressSig$Gene_Human]
VariableFeatures(HHA.Bulb) %>% length() #1938 genes genes

#-------Performing Scanpy Conversion for SCVI Integration: Cell Cycle and Stress Removed-----#
#Subsetting Atlas for Variable Genes
Diet.Bulb <- HHA.Bulb[VariableFeatures(HHA.Bulb)]
#----Writing Metadata----#
Bulb.metadata <- Diet.Bulb[[]]
Bulb.metadata$barcode <- rownames(Bulb.metadata)
write.csv(Bulb.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_Metadata_Cleaned_5_22_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.Bulb) <- "RNA"
#Writing log normalized counts to mtx
DietBulb.counts <- Diet.Bulb %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietBulb.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_RawCountsCleaned_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietBulb.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_GeneNames_Cleaned_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.Bulb, DietBulb.counts, Bulb.metadata)
#------Writing Normalized Cell Cycle Gene Expression------#
#Writing as csv
CellCycle.GeneMat <- Bulb.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_CellCycle_LogNorm_Cleaned_5_22_25.csv"))
#------Writing Normalized Stress Gene Expression------#
StressSig.GeneMat <- Bulb.StressSig %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
StressSig.GeneMat$Barcode <- rownames(StressSig.GeneMat)
#Writing to csv
write_csv(StressSig.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Bulb_StressSig_LogNorm_Cleaned_5_22_25.csv"))

#-----Getting SCVI Integration Results-----#
Bulb.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Bulb_PostCleaning_LatentRep_5_22_25.csv")) %>% column_to_rownames("Barcode")
head(Bulb.scvi)
# put it back in our original Seurat object
HHA.Bulb[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Bulb.scvi), key = "scvi_")
rm(Bulb.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Bulb <-FindNeighbors(HHA.Bulb, dims = 1:50, reduction = "scvi")
HHA.Bulb <- RunUMAP(HHA.Bulb, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Bulb, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", raster = F, pt.size = 0.5) 
#-----Plotting by Sample-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T)
#-----Plotting by Donor-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T, pt.size = 0.5) 
#------Plotting By Platform-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) 
#-----Plotting By Annotation-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T)
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("scDblFinder.score"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-------Plotting Markers--------#
#Inner layer and ORS
FeaturePlot(HHA.Bulb, features = c("KRT5", "COL17A1", "MKI67", "LGR5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Inner Root Sheath
FeaturePlot(HHA.Bulb, features = c("GATA3", "TCHH", "KRT72", "KRT74"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cortex
FeaturePlot(HHA.Bulb, features = c("ODC1", "KRT5", "KRT35", "KRT85", "LEF1", "KRT31"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cuticle
FeaturePlot(HHA.Bulb, features = c("KRT32", "KRT82", "CRAT", "SOX21"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proximal Companion and Medulla
FeaturePlot(HHA.Bulb, features = c("KRT75", "TCHH", "FOXQ1", "ALDH3A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#ORS and Companion Beyond Matrix
FeaturePlot(HHA.Bulb, features = c("KRT75", "FGF5", "KRT6A", "BARX2"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proliferation
FeaturePlot(HHA.Bulb, features = c("MKI67"), reduction = "UMAP", pt.size = 1) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(HHA.Bulb, features = c("COL17A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Problem Genes
FeaturePlot(HHA.Bulb, features = c("FOS", "JUN", "KRT35", "KRT5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Violin Plot of Problematic Genes and Dissociation Score
VlnPlot(HHA.Bulb, group.by = "scvi_clusters_0.75", features = c("FOS", "JUN", "DissociationScore", "nCount_RNA", "nFeature_RNA", "percent.mt"), alpha = 0.1)

#--------Removing Problem Clusters-------#
#Clusters 12 and 10 are holdouts from the original bad cluster removal: 10 has abnormally high counts, and is probably holdouts that failed to reintegrate well, 12 shows clear elevation of FOS/Jun/Dissociation Score
#Cluster 11 is BARX2 and KRT6A positive, which maps well above the LPC and proximal companion, it represents the companion layer with some early suprabasal ORS and is too far up for the analysis we want to do here on fate specification. 
#There is a really clear set of poor quality cells within the COL17A1 + layer that is a clear outlier for dissociation scores
#-------Clustering and removing Outliers from  Cluster 2-------#
DissociationCleanup <- subset(HHA.Bulb, scvi_clusters_0.75 == 2)
DissociationCleanup
DissociationCleanup<- FindClusters(DissociationCleanup, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
DissociationCleanup<- FindClusters(DissociationCleanup, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
DissociationCleanup<- FindClusters(DissociationCleanup, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
#Cluster 1 represents the outliers and can be removed
DimPlot(DissociationCleanup, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 0.5) 
DissociationCleanup <- subset(DissociationCleanup, scvi_clusters_0.4 == 1)
#--------Filtering for Problem Clusters------#
HHA.Bulb$Barcode <- colnames(HHA.Bulb)
#Removing outlier FOS/JUN cluster in cluster 2
HHA.Bulb <- subset(HHA.Bulb, Barcode %in% colnames(DissociationCleanup) == F)
#Removing other problematic clusters
HHA.Bulb <- subset(HHA.Bulb, scvi_clusters_0.75 %in% c(10, 11, 12) == F)
HHA.Bulb #15,531 cells are left
DimPlot(HHA.Bulb, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", raster = F, pt.size = 0.5) 
#-----Saving list of retained cells to csv-----#
#Not going to rerun variable feature selection here, as changes should be pretty minimal 
Bulb.FilteredBarcodes <- data.frame(Barcode = colnames(HHA.Bulb))
write_csv(Bulb.FilteredBarcodes, file = paste0(AnnDataDirectory, "HHA_Bulb_Doublet_Filtered_Barcodes_5_22_25.csv"))

#---------------One Last Ride-----------#
#-----Getting SCVI Integration Results-----#
Bulb.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Bulb_Filtered_LatentRep_5_22_25.csv")) %>% column_to_rownames("Barcode")
head(Bulb.scvi)
# put it back in our original Seurat object
HHA.Bulb[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Bulb.scvi), key = "scvi_")
rm(Bulb.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Bulb <-FindNeighbors(HHA.Bulb, dims = 1:50, reduction = "scvi", return.neighbor = T)
HHA.Bulb <- RunUMAP(HHA.Bulb, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Bulb<- FindClusters(HHA.Bulb, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Bulb, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.6", raster = F, pt.size = 0.5) 
#-----Plotting by Sample-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T, pt.size = 0.5)
#-----Plotting by Donor-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T, pt.size = 0.5) 
#------Plotting By Platform-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T,
        cols = c("steelblue3", "firebrick", "darkseagreen"), pt.size = 0.5) 
#-----Plotting By Annotation-----#
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T, pt.size = 0.5)
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Bulb, reduction = "UMAP", features = c("scDblFinder.score"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-------Plotting Markers--------#
#Inner layer and ORS
FeaturePlot(HHA.Bulb, features = c("KRT5", "COL17A1", "MKI67", "LGR5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Inner Root Sheath
FeaturePlot(HHA.Bulb, features = c("GATA3", "TCHH", "KRT72", "KRT74"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cortex
FeaturePlot(HHA.Bulb, features = c("ODC1", "KRT5", "KRT35", "KRT85", "LEF1", "KRT31"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cuticle
FeaturePlot(HHA.Bulb, features = c("KRT32", "KRT82", "CRAT", "SOX21"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proximal Companion and Medulla
FeaturePlot(HHA.Bulb, features = c("KRT75", "TCHH", "FOXQ1", "ALDH3A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#ORS and Companion Beyond Matrix
FeaturePlot(HHA.Bulb, features = c("KRT75", "LGR5", "KRT6A", "BARX2"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proliferation
FeaturePlot(HHA.Bulb, features = c("MKI67"), reduction = "UMAP", pt.size = 1) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(HHA.Bulb, features = c("COL17A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Problem Genes
FeaturePlot(HHA.Bulb, features = c("FOS", "JUN", "KRT35", "KRT5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Violin Plot of Problematic Genes and Dissociation Score
VlnPlot(HHA.Bulb, group.by = "scvi_clusters_0.75", features = c("FOS", "JUN", "DissociationScore", "nCount_RNA", "nFeature_RNA", "percent.mt"), alpha = 0.1)

#-------Subdividing Matrix into Zones-------#
#Each zone will be subclustered in order to map individual layers. 
BulbAnnotationsBroad <- c(
  "IRS", #0
  "Cortex", #1
  "Cortex", #2
  "Cortex", #3
  "Cuticle", #4
  "Inner Layer", #5
  "Cortex", #6
  "LPC", #6
  "Proximal Companion") #8

#-----Renaming Idents-----#
Idents(HHA.Bulb) <- "scvi_clusters_0.6"
names(BulbAnnotationsBroad) <- levels(HHA.Bulb)
HHA.Bulb <- RenameIdents(HHA.Bulb, BulbAnnotationsBroad)
HHA.Bulb$BulbAnnotationBroad <- HHA.Bulb@active.ident
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "BulbAnnotationBroad", label = T, repel = T, raster = F, pt.size = 1.0)
table(HHA.Bulb$BulbAnnotationBroad)
table(HHA.Bulb$BulbAnnotationBroad, HHA.Bulb$Platform)

#------Subclustering Analysis of Matrix-------#
#-------Inner Layer----------#
#This stands for germinative layer, and maps to the COL17A1+ cells along the DP
HHA.GL <- subset(HHA.Bulb, BulbAnnotationBroad == "Inner Layer")
HHA.GL #1574 cells
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.GL <-FindNeighbors(HHA.GL, dims = 1:50, reduction = "scvi")
HHA.GL <- RunUMAP(HHA.GL, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.GL<- FindClusters(HHA.GL, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.GL<- FindClusters(HHA.GL, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.GL<- FindClusters(HHA.GL, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.GL<- FindClusters(HHA.GL, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.GL<- FindClusters(HHA.GL, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(HHA.GL, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.6", raster = F, pt.size = 3) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(HHA.GL, reduction = "UMAP", label = F, group.by = "Platform", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(HHA.GL, reduction = "UMAP", features = c("DissociationScore"), pt.size = 3) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-----Finding Markers for Each Cluster-----#
GLMarkers.TFIDF <- HHA.GL %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.GL$scvi_clusters_0.35, N = 20, FDR = 0.01)
GLMarkers.TFIDF
GLMarkers.TFIDFList <- GLMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(GLMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.GL$scvi_clusters_0.35))-1))
GLMarkers.TFIDFList

#-------Wilcoxon Test-------#
Idents(HHA.GL) <- "scvi_clusters_0.35"
GLMarkers.Wilcox <- HHA.GL %>% JoinLayers() %>% FindAllMarkers(assay = "RNA", test.use = "wilcox", only_pos = T, verbose = T)
GLMarkers.Wilcox %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% group_split()

#-------Plotting Markers------#
FeaturePlot(HHA.GL, features = c("NFIB", "COL17A1", "FOXQ1", "KRT14", "KRT75", "CPA6", "GATA3", "FRY", "SHH"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Lower Cluster
FeaturePlot(HHA.GL, features = c("FRY", "SHH", "GATA3", "SAMD5"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
FeaturePlot(HHA.GL, features = c("ODC1", "CDH13", "TNC", "LAMB3"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
FeaturePlot(HHA.GL, features = c("MKI67"), reduction = "UMAP")  &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-------Annotating Clusters------#
GLAnnotations <- c("Delaminating COl17", #0
                   "Lower COL17", #1
                   "Upper COL17", #2
                   "Medulla", #3
                   "Middle COL17", #4
                   "Medulla") #5
#-----Renaming Idents-----#
Idents(HHA.GL) <- "scvi_clusters_0.6"
names(GLAnnotations) <- levels(HHA.GL)
HHA.GL <- RenameIdents(HHA.GL, GLAnnotations)
HHA.GL$GLAnnotation <- HHA.GL@active.ident
DimPlot(HHA.GL, reduction = "UMAP", group.by = "GLAnnotation", label = T, repel = T, raster = F, pt.size = 3.0)
table(HHA.GL$GLAnnotation)

#--------Exporting to Python for cNMF Analysis------#
#----Writing Metadata----#
GL.metadata <- HHA.GL[[]] %>% mutate(UMAP_1 = HHA.GL@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.GL@reductions$UMAP@cell.embeddings[,2])
GL.metadata$barcode <- rownames(GL.metadata)
GL.metadata <- GL.metadata
write.csv(GL.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_GL_FullConversion_Metadata_5_22_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(HHA.GL@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_GL_FullConversion_scVI_Embeddings_5_22_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Pulling Raw Counts Data
GL.rawcounts <- HHA.GL %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(GL.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_GL_FullConversion_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(GL.rawcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_GL_FullConversion_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)

#------Adding Idents to Seurat Object------$
GL.Metadata <- data.frame(Barcode = colnames(HHA.GL), GLAnnotation = HHA.GL$GLAnnotation)
GL.Metadata
GL.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Bulb)), data.frame(Barcode = colnames(HHA.GL), GLAnnotation = HHA.GL$GLAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Bulb$GLAnnotation <- GL.Metadata$GLAnnotation
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "GLAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#------Subclustering Analysis of Matrix-------#
#-------Inner Layer----------#
#This stands for germinative layer, and maps to the COL17A1+ cells along the DP
HHA.IRS <- subset(HHA.Bulb, BulbAnnotationBroad == "IRS")
HHA.IRS #1574 cells
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.IRS <-FindNeighbors(HHA.IRS, dims = 1:50, reduction = "scvi")
HHA.IRS <- RunUMAP(HHA.IRS, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.IRS<- FindClusters(HHA.IRS, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(HHA.IRS, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.5", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(HHA.IRS, reduction = "UMAP", label = F, group.by = "Platform", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(HHA.IRS, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-----Finding Markers for Each Cluster-----#
IRSMarkers.TFIDF <- HHA.IRS %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.IRS$scvi_clusters_0.5, N = 20, FDR = 0.01)
IRSMarkers.TFIDF
IRSMarkers.TFIDFList <- IRSMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(IRSMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.IRS$scvi_clusters_0.5))-1))
IRSMarkers.TFIDFList

#-------Plotting Markers------#
FeaturePlot(HHA.IRS, features = c("KRT73", "KRT74", "TCHH", "GATA3","SHH", "SOX21"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-------Annotating Clusters------#
IRSAnnotations <- c("IRS Cuticle Progenitor", #0
                    "Huxley Progenitor", #1
                    "IRS Cuticle", #2
                    "IRS Henle", #3
                    "IRS Huxley") #4
#-----Renaming Idents-----#
Idents(HHA.IRS) <- "scvi_clusters_0.5"
names(IRSAnnotations) <- levels(HHA.IRS)
HHA.IRS <- RenameIdents(HHA.IRS, IRSAnnotations)
HHA.IRS$IRSAnnotation <- HHA.IRS@active.ident
DimPlot(HHA.IRS, reduction = "UMAP", group.by = "IRSAnnotation", label = T, repel = T, raster = F, pt.size = 3.0)
table(HHA.IRS$IRSAnnotation)

#------Adding Idents to Seurat Object------$
IRS.Metadata <- data.frame(Barcode = colnames(HHA.IRS), IRSAnnotation = HHA.IRS$IRSAnnotation)
IRS.Metadata
IRS.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Bulb)), data.frame(Barcode = colnames(HHA.IRS), IRSAnnotation = HHA.IRS$IRSAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Bulb$IRSAnnotation <- IRS.Metadata$IRSAnnotation
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "IRSAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#------Subclustering Analysis of Matrix-------#
#-------Inner Layer----------#
#This stands for germinative layer, and maps to the COL17A1+ cells along the DP
HHA.Cortex <- subset(HHA.Bulb, BulbAnnotationBroad == "Cortex")
HHA.Cortex #1574 cells
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Cortex <-FindNeighbors(HHA.Cortex, dims = 1:50, reduction = "scvi")
HHA.Cortex <- RunUMAP(HHA.Cortex, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.Cortex<- FindClusters(HHA.Cortex, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(HHA.Cortex, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(HHA.Cortex, reduction = "UMAP", label = F, group.by = "Platform", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(HHA.Cortex, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
DimPlot(HHA.Cortex, reduction = "UMAP", group.by = "DonorID", pt.size = 1) 

#-----Finding Markers for Each Cluster-----#
CortexMarkers.TFIDF <- HHA.Cortex %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Cortex$scvi_clusters_0.3, N = 20, FDR = 0.01)
CortexMarkers.TFIDF
CortexMarkers.TFIDFList <- CortexMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(CortexMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Cortex$scvi_clusters_0.3))-1))
CortexMarkers.TFIDFList

#-------Plotting Markers------#
FeaturePlot(HHA.Cortex, features = c("KRT5", "MKI67", "KRT35", "LEF1", "KRT31", "KRT75"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-------Annotating Clusters------#
CortexAnnotations <- c("Middle Cortex",
                       "Early Cortex",
                       "Late Cortex") #4
#-----Renaming Idents-----#
Idents(HHA.Cortex) <- "scvi_clusters_0.25"
names(CortexAnnotations) <- levels(HHA.Cortex)
HHA.Cortex <- RenameIdents(HHA.Cortex, CortexAnnotations)
HHA.Cortex$CortexAnnotation <- HHA.Cortex@active.ident
DimPlot(HHA.Cortex, reduction = "UMAP", group.by = "CortexAnnotation", label = T, repel = T, raster = F, pt.size = 3.0)
table(HHA.Cortex$CortexAnnotation)

#------Adding Idents to Seurat Object------$
Cortex.Metadata <- data.frame(Barcode = colnames(HHA.Cortex), CortexAnnotation = HHA.Cortex$CortexAnnotation)
Cortex.Metadata
Cortex.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Bulb)), data.frame(Barcode = colnames(HHA.Cortex), CortexAnnotation = HHA.Cortex$CortexAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Bulb$CortexAnnotation <- Cortex.Metadata$CortexAnnotation
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "CortexAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#------Subclustering Analysis of Matrix-------#
#-------Inner Layer----------#
#This stands for germinative layer, and maps to the COL17A1+ cells along the DP
HHA.Cuticle <- subset(HHA.Bulb, BulbAnnotationBroad == "Cuticle")
HHA.Cuticle #1574 cells
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Cuticle <-FindNeighbors(HHA.Cuticle, dims = 1:50, reduction = "scvi")
HHA.Cuticle <- RunUMAP(HHA.Cuticle, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.Cuticle<- FindClusters(HHA.Cuticle, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(HHA.Cuticle, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(HHA.Cuticle, reduction = "UMAP", label = F, group.by = "Platform", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(HHA.Cuticle, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
DimPlot(HHA.Cuticle, reduction = "UMAP", group.by = "DonorID", pt.size = 1) 

#-----Finding Markers for Each Cluster-----#
CuticleMarkers.TFIDF <- HHA.Cuticle %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Cuticle$scvi_clusters_0.3, N = 20, FDR = 0.01)
CuticleMarkers.TFIDF
CuticleMarkers.TFIDFList <- CuticleMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(CuticleMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Cuticle$scvi_clusters_0.3))-1))
CuticleMarkers.TFIDFList

#-------Plotting Markers------#
FeaturePlot(HHA.Cuticle, features = c("KRT5", "MKI67", "KRT35", "LEF1", "KRT31", "KRT75", "GATA3", "SOX21", "SHH"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-------Annotating Clusters------#
CuticleAnnotations <- c("Early Cuticle",
                        "Late Cuticle")
#-----Renaming Idents-----#
Idents(HHA.Cuticle) <- "scvi_clusters_0.25"
names(CuticleAnnotations) <- levels(HHA.Cuticle)
HHA.Cuticle <- RenameIdents(HHA.Cuticle, CuticleAnnotations)
HHA.Cuticle$CuticleAnnotation <- HHA.Cuticle@active.ident
DimPlot(HHA.Cuticle, reduction = "UMAP", group.by = "CuticleAnnotation", label = T, repel = T, raster = F, pt.size = 3.0)
table(HHA.Cuticle$CuticleAnnotation)

#------Adding Idents to Seurat Object------$
Cuticle.Metadata <- data.frame(Barcode = colnames(HHA.Cuticle), CuticleAnnotation = HHA.Cuticle$CuticleAnnotation)
Cuticle.Metadata
Cuticle.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Bulb)), data.frame(Barcode = colnames(HHA.Cuticle), CuticleAnnotation = HHA.Cuticle$CuticleAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Bulb$CuticleAnnotation <- Cuticle.Metadata$CuticleAnnotation
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "CuticleAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#------Subclustering Analysis of Matrix-------#
#-------Inner Layer----------#
#This stands for germinative layer, and maps to the COL17A1+ cells along the DP
HHA.LPCComp <- subset(HHA.Bulb, BulbAnnotationBroad %in% c("LPC", "Proximal Companion"))
HHA.LPCComp #1574 cells
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.LPCComp <-FindNeighbors(HHA.LPCComp, dims = 1:50, reduction = "scvi")
HHA.LPCComp <- RunUMAP(HHA.LPCComp, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.LPCComp<- FindClusters(HHA.LPCComp, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(HHA.LPCComp, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.2", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(HHA.LPCComp, reduction = "UMAP", label = F, group.by = "Platform", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(HHA.LPCComp, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
DimPlot(HHA.LPCComp, reduction = "UMAP", group.by = "DonorID", pt.size = 1) 

#-----Finding Markers for Each Cluster-----#
LPCCompMarkers.TFIDF <- HHA.LPCComp %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.LPCComp$scvi_clusters_0.3, N = 20, FDR = 0.01)
LPCCompMarkers.TFIDF
LPCCompMarkers.TFIDFList <- LPCCompMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(LPCCompMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.LPCComp$scvi_clusters_0.3))-1))
LPCCompMarkers.TFIDFList

#-------Plotting Markers------#
FeaturePlot(HHA.LPCComp, features = c("KRT5", "KRT75", "LGR5", "BARX2", "KRT6A", "SORCS2"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-------Annotating Clusters------#
LPCCompAnnotations <- c("Early LPC", #0
                        "Outer Progenitor", #1
                        "Early Proximal Companion", #2
                        "Late LPC", #3
                        "Late Proximal Companion") #4
LPCCompAnnotations <- c("Outer Progenitor", #1
                        "Proximal Companion", #2
                        "LPC") #3
#-----Renaming Idents-----#
Idents(HHA.LPCComp) <- "scvi_clusters_0.2"
names(LPCCompAnnotations) <- levels(HHA.LPCComp)
HHA.LPCComp <- RenameIdents(HHA.LPCComp, LPCCompAnnotations)
HHA.LPCComp$LPCCompAnnotation <- HHA.LPCComp@active.ident
DimPlot(HHA.LPCComp, reduction = "UMAP", group.by = "LPCCompAnnotation", label = T, repel = T, raster = F, pt.size = 3.0)
table(HHA.LPCComp$LPCCompAnnotation)

#------Adding Idents to Seurat Object------$
LPCComp.Metadata <- data.frame(Barcode = colnames(HHA.LPCComp), LPCCompAnnotation = HHA.LPCComp$LPCCompAnnotation)
LPCComp.Metadata
LPCComp.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Bulb)), data.frame(Barcode = colnames(HHA.LPCComp), LPCCompAnnotation = HHA.LPCComp$LPCCompAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Bulb$LPCCompAnnotation <- LPCComp.Metadata$LPCCompAnnotation
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "LPCCompAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#--------Adding Column with All Clusters--------#
HHA.Bulb$BulbAnnotationFine <- case_when(HHA.Bulb$BulbAnnotationBroad == "IRS" ~ HHA.Bulb$IRSAnnotation,
                                         HHA.Bulb$BulbAnnotationBroad == "Inner Layer" ~ HHA.Bulb$GLAnnotation,
                                         HHA.Bulb$BulbAnnotationBroad %in% c("LPC", "Proximal Companion") ~ HHA.Bulb$LPCCompAnnotation,
                                         HHA.Bulb$BulbAnnotationBroad == "Cortex" ~ HHA.Bulb$CortexAnnotation,
                                         HHA.Bulb$BulbAnnotationBroad == "Cuticle" ~ HHA.Bulb$CuticleAnnotation,)
DimPlot(HHA.Bulb, reduction = "UMAP", group.by = "BulbAnnotationFine", label = T, repel = T, raster = F, pt.size = 1.0, cols = DiscretePalette(22)) 

#---------Saving Checkpoint-------#
save.image(paste0(WorkspaceDirectory, "HHA_Bulb_Recluster_5_22_25.RData"))

#--------Exporting to Python for Downstream Analysis------#
#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
Bulb.metadata <- HHA.Bulb[[]] %>% mutate(UMAP_1 = HHA.Bulb@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.Bulb@reductions$UMAP@cell.embeddings[,2])
Bulb.metadata$barcode <- rownames(Bulb.metadata)
Bulb.metadata <- Bulb.metadata
write.csv(Bulb.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Bulb_FullConversion_Metadata_5_22_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(HHA.Bulb@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Bulb_FullConversion_scVI_Embeddings_5_22_25.csv"), quote = F, row.names = F)

#------Exporting SNN Graph Connectivities------#
#Getting long format DF containing all possible connections and their presence/absence
Bulb.SNN <- as.data.frame(summary(HHA.Bulb@graphs$RNA_nn)) %>% #Matching Pythonic Indexing (starting at 0)
  mutate(i = i-1, j = j-1)
write_csv(Bulb.SNN, file = paste0(AnnDataDirectory, "HHA_Multiome_Bulb_SNN_Connectivities_5_22_25.csv"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Bulb.logcounts <- HHA.Bulb %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(Bulb.logcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Bulb_FullConversion_Expression_LogCounts_5_22_25.mtx"))
#Pulling Raw Counts Data
Bulb.rawcounts <- HHA.Bulb %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(Bulb.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Bulb_FullConversion_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Bulb.logcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Bulb_FullConversion_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)

#------Pulling in Force Atlas 2 Graph from Scanpy-------#
Bulb.FA2 <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Bulb_ForceAtlas_4_1_25.csv")) %>% column_to_rownames("Barcode")
# put it back in our original Seurat object
HHA.Bulb[["FA2"]] <- CreateDimReducObject(embeddings = as.matrix(Bulb.FA2), key = "FA_")
rm(Bulb.FA2)
DimPlot(HHA.Bulb, reduction = "FA2", label = T, group.by = "BulbAnnotation", raster = F, pt.size = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2")

#---------Saving Barcodes for ATAC analysis-------#
#Subsetting for Multiome Cells
Bulb.Multiome <- HHA.Bulb %>% subset(Platform == "Multiome")
#Writing to csv for ArchR Analysis
Bulb.ArchR <- Bulb.Multiome[[]] %>% dplyr::select(UnifiedBarcode, BulbAnnotationBroad, BulbAnnotationFine)
write_csv(Bulb.ArchR, paste0(SeuratDirectory, 'HHA_Multiome_Bulb_Barcodes.csv'), col_names = T)

#-------Subsetting out Proximal Companion and ORS----------#
HHA.Matrix <- subset(HHA.Bulb, BulbAnnotationFine %in% c("LPC", "Proximal Companion") == F)
HHA.Matrix
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "BulbAnnotationFine",
        label = T, repel = T, raster = F, pt.size = 1.0, cols = DiscretePalette(22)) 

#-----Reclustering after cleaning -----#
#Selecting 2000 variable features for dimensionality reduction
HHA.Matrix <- FindVariableFeatures(HHA.Matrix, selection.method = "vst", nfeatures = 2000)
HHA.Matrix.TopFeaturePlot <- VariableFeaturePlot(HHA.Matrix) %>% LabelPoints(points = head(VariableFeatures(HHA.Matrix), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Matrix.TopFeaturePlot
#Scaling Data
HHA.Matrix <- ScaleData(HHA.Matrix)
#Running PCA
HHA.Matrix <- RunPCA(HHA.Matrix)
DimHeatmap(HHA.Matrix, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Matrix <-FindNeighbors(HHA.Matrix, dims = 1:30, reduction = "pca")
HHA.Matrix <- RunUMAP(HHA.Matrix, dims = 1:30, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Matrix <- FindClusters(HHA.Matrix, resolution = 1, cluster.name = "unintegrated_clusters")

#------Getting Stress Signature Associated Genes-----#
#Subsetting Seurat Object for genes in stress signature
Matrix.StressSig <- HHA.Matrix[StressSig$Gene_Human,]
#------Pulling Cell Cycle Associated Genes-----#
#Subsetting Seurat Object for genes in stress signature
Matrix.CellCycle <- HHA.Matrix[CanonicalCellCycle,]

#------Removing Dissociation Signature and Cell Cycle genes from Variable Features------#
#Removing Stress Signature Genes
VariableFeatures(HHA.Matrix) <- VariableFeatures(HHA.Matrix)[!VariableFeatures(HHA.Matrix) %in% StressSig$Gene_Human]
VariableFeatures(HHA.Matrix) %>% length() #1979 genes
#Removing CellCycle Signature Genes
VariableFeatures.CCIncluded <- VariableFeatures(HHA.Matrix)
VariableFeatures(HHA.Matrix) <- VariableFeatures(HHA.Matrix)[!VariableFeatures(HHA.Matrix) %in% CanonicalCellCycle]
VariableFeatures(HHA.Matrix) %>% length() #1932 genes

#-------Performing Scanpy Conversion for SCVI Integration: Cell Cycle and Stress Removed-----#
#Subsetting Atlas for Variable Genes
Diet.Matrix <- HHA.Matrix[VariableFeatures(HHA.Matrix)]
#----Writing Metadata----#
Matrix.metadata <- Diet.Matrix[[]]
Matrix.metadata$barcode <- rownames(Matrix.metadata)
write.csv(Matrix.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_Metadata_5_22_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.Matrix) <- "RNA"
#Writing raw counts to mtx
DietMatrix.counts <- Diet.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietMatrix.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietMatrix.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.Matrix, DietMatrix.counts, Matrix.metadata)
#----Writing Count Matrix (Cell Cycle Included)------#
Diet.CC <- HHA.Matrix[VariableFeatures.CCIncluded]
Diet.CC
DefaultAssay(Diet.CC) <- "RNA"
#Writing raw counts to mtx
DietCC.counts <- Diet.CC %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietCC.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_CC_Included_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietCC.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_CC_Included_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)
#------Writing Normalized Cell Cycle Gene Expression------#
#Writing as csv
CellCycle.GeneMat <- Matrix.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_CellCycle_LogNorm_5_22_25.csv"))
#------Writing Normalized Stress Gene Expression------#
StressSig.GeneMat <- Matrix.StressSig %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
StressSig.GeneMat$Barcode <- rownames(StressSig.GeneMat)
#Writing to csv
write_csv(StressSig.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Matrix_StressSig_LogNorm_5_22_25.csv"))

#-----Getting SCVI Integration Results-----#
Matrix.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Matrix_LatentRep_5_22_25.csv")) %>% column_to_rownames("Barcode")
head(Matrix.scvi)
# put it back in our original Seurat object
HHA.Matrix[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Matrix.scvi), key = "scvi_")
rm(Matrix.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Matrix <-FindNeighbors(HHA.Matrix, dims = 1:50, reduction = "scvi")
HHA.Matrix <- RunUMAP(HHA.Matrix, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Matrix, reduction = "UMAP", label = T, group.by = "BulbAnnotationBroad", raster = F, pt.size = 0.5) 
#-----Plotting by Sample-----#
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T)
#-----Plotting by Donor-----#
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T, pt.size = 0.5) 
#------Plotting By Platform-----#
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) 
#-----Plotting By Annotation-----#
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "BulbAnnotationBroad", raster = F, shuffle = T)
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "BulbAnnotationFine", raster = F, shuffle = T)
#Plotting QC Metrics
FeaturePlot(HHA.Matrix, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Matrix, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Matrix, reduction = "UMAP", features = c("scDblFinder.score"), pt.size = 0.5) &
  scale_color_viridis_c(option = "turbo") & NoAxes()

#-------Plotting Markers--------#
#Inner layer and ORS
FeaturePlot(HHA.Matrix, features = c("COL17A1", "ITGA6", "NFIB", "RUNX1"), reduction = "UMAP", pt.size = 0.6) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Inner layer and ORS
FeaturePlot(HHA.Matrix, features = c("KRT5", "COL17A1", "MKI67", "LGR5"), reduction = "UMAP", pt.size = 0.6) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Inner Root Sheath
FeaturePlot(HHA.Matrix, features = c("GATA3", "TCHH", "KRT72", "KRT74"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cortex
FeaturePlot(HHA.Matrix, features = c("ODC1", "KRT5", "KRT35", "KRT85", "LEF1", "KRT31"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Cuticle
FeaturePlot(HHA.Matrix, features = c("KRT32", "KRT82", "CRAT", "SOX21"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proximal Companion and Medulla
FeaturePlot(HHA.Matrix, features = c("KRT75", "TCHH", "FOXQ1", "ALDH3A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#ORS and Companion Beyond Matrix
FeaturePlot(HHA.Matrix, features = c("KRT75", "FGF5", "KRT6A", "BARX2"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Proliferation
FeaturePlot(HHA.Matrix, features = c("MKI67"), reduction = "UMAP", pt.size = 1) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(HHA.Matrix, features = c("COL17A1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Differentiation Drivers
FeaturePlot(HHA.Matrix, features = c("GATA3", "SHH", "FRY", "HOXC13", "LEF1", "SOX21"), reduction = "UMAP", ncol = 3) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#--------Exporting to Python for Downstream Analysis------#
#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
Matrix.metadata <- HHA.Matrix[[]] %>% mutate(UMAP_1 = HHA.Matrix@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.Matrix@reductions$UMAP@cell.embeddings[,2])
Matrix.metadata$barcode <- rownames(Matrix.metadata)
Matrix.metadata <- Matrix.metadata
write.csv(Matrix.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_Metadata_5_22_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(HHA.Matrix@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_scVI_Embeddings_5_22_25.csv"), quote = F, row.names = F)

#------Exporting SNN Graph Connectivities------#
#Getting long format DF containing all possible connections and their presence/absence
Matrix.SNN <- as.data.frame(summary(HHA.Matrix@graphs$RNA_nn)) %>% #Matching Pythonic Indexing (starting at 0)
  mutate(i = i-1, j = j-1)
write_csv(Matrix.SNN, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_SNN_Connectivities_5_22_25.csv"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Matrix.logcounts <- HHA.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(Matrix.logcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_Expression_LogCounts_5_22_25.mtx"))
#Pulling Raw Counts Data
Matrix.rawcounts <- HHA.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(Matrix.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_RawCounts_5_22_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Matrix.logcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_GeneNames_5_22_25.csv"),
            quote = F, row.names = F, col.names = T)

#------Pulling in Force Atlas 2 Graph from Scanpy-------#
Matrix.FA2 <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Matrix_FA_Projection_5_22_25.csv")) %>% column_to_rownames("Barcode")
# put it back in our original Seurat object
HHA.Matrix[["FA2"]] <- CreateDimReducObject(embeddings = as.matrix(Matrix.FA2), key = "FA_")
rm(Matrix.FA2)
DimPlot(HHA.Matrix, reduction = "FA2", label = T, group.by = "BulbAnnotationBroad", raster = F, pt.size = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Matrix, reduction = "FA2", label = T, group.by = "scvi_clusters_0.6", raster = F, pt.size = 0.5) 
DimPlot(HHA.Matrix, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", raster = F, pt.size = 0.5) 
#Differentiation Drivers
FeaturePlot(HHA.Matrix, features = c("GATA3", "SHH", "FRY", "HOXC13", "LEF1", "SOX21"), reduction = "FA2", ncol = 3) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#------Grouping Cells into Broad Lineages------#
MatrixAnnotationsBroad <- c("Cortex", #0
                           "Cortex", #1
                           "IRS", #2
                           "Cuticle", #3
                           "Cortex", #4
                           "Cortex", #5
                           "COL17", #6
                           "Cortex", #7
                           "IRS", #8
                           "COL17", #9
                           "LPC", #10
                           "IRS", #11
                           "Cuticle", #12
                           "IRS") #13
#-----Renaming Idents-----#
Idents(HHA.Matrix) <- "scvi_clusters_0.75"
names(MatrixAnnotationsBroad) <- levels(HHA.Matrix)
HHA.Matrix <- RenameIdents(HHA.Matrix, MatrixAnnotationsBroad)
HHA.Matrix$MatrixAnnotationBroad <- HHA.Matrix@active.ident
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "MatrixAnnotationBroad", label = T, repel = T, raster = F, pt.size = 1.0)
DimPlot(HHA.Matrix, reduction = "FA2", group.by = "MatrixAnnotationBroad", label = T, repel = T, raster = F, pt.size = 1.0)

#-------Resolving IRS and Curticle Lineages-------#
Matrix.GL <- subset(HHA.Matrix, MatrixAnnotationBroad %in% c("COL17", "LPC"))
Matrix.GL
#-----Performing Final DR and Clustering in scVI Space-----#
Matrix.GL <-FindNeighbors(Matrix.GL, dims = 1:50, reduction = "scvi")
#Matrix.GL <- RunUMAP(Matrix.GL, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 1, cluster.name = "scvi_clusters_1.0")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
Matrix.GL<- FindClusters(Matrix.GL, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(Matrix.GL, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.6", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(Matrix.GL, reduction = "UMAP", label = F, group.by = "BulbAnnotationFine", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(Matrix.GL, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting Markers
FeaturePlot(Matrix.GL, features = c("GATA3", "SHH", "FRY", "HOXC13", "LEF1", "COL17A1", "MKI67", "SOX21", "KRT35"), reduction = "UMAP", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(Matrix.GL, features = c("GATA3", "SHH", "FRY", "HOXC13", "LEF1", "COL17A1", "MKI67", "SOX21", "KRT35"), reduction = "FA2", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-------Annotating Clusters------#
GLAnnotations <- c("Lower COL17", #0
                   "Upper COL17", #1
                   "LPC", #2,
                   "Upper COL17", #3
                   "LPC", #4
                   "Early_Cortex",
                   "Medulla") #5
#-----Renaming Idents-----#
Idents(Matrix.GL) <- "scvi_clusters_0.6"
names(GLAnnotations) <- levels(Matrix.GL)
Matrix.GL <- RenameIdents(Matrix.GL, GLAnnotations)
Matrix.GL$GLAnnotation <- Matrix.GL@active.ident
Matrix.GL$GLAnnotationOrdered <- factor(Matrix.GL$GLAnnotation, levels = c("LPC", "Lower COL17", "Upper COL17", "Medulla", "Early_Cortex"))
DimPlot(Matrix.GL, reduction = "UMAP", group.by = "GLAnnotationOrdered", label = T, repel = T, raster = F)
table(Matrix.GL$GLAnnotation)

#------Adding Idents to Seurat Object------$
GL.Metadata <- data.frame(Barcode = colnames(Matrix.GL), GLAnnotation = Matrix.GL$GLAnnotation)
GL.Metadata
GL.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Matrix)), data.frame(Barcode = colnames(Matrix.GL), GLAnnotation = Matrix.GL$GLAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Matrix$GLAnnotation <- GL.Metadata$GLAnnotation
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "GLAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#-------Resolving IRS and Curticle Lineages-------#
Matrix.IRS <- subset(HHA.Matrix, MatrixAnnotationBroad %in% c("IRS"))
Matrix.IRS
#-----Performing Final DR and Clustering in scVI Space-----#
Matrix.IRS <-FindNeighbors(Matrix.IRS, dims = 1:50, reduction = "scvi")
#Matrix.IRS <- RunUMAP(Matrix.IRS, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 1, cluster.name = "scvi_clusters_1.0")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
Matrix.IRS<- FindClusters(Matrix.IRS, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(Matrix.IRS, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(Matrix.IRS, reduction = "FA2", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 1)
DimPlot(Matrix.IRS, reduction = "UMAP", label = F, group.by = "BulbAnnotationFine", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(Matrix.IRS, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting Markers
FeaturePlot(Matrix.IRS, features = c("GATA3", "SHH", "FRY", "TCHH", "SOX21", "KRT73", "KRT74"), reduction = "UMAP", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(Matrix.IRS, features = c("GATA3", "SHH", "FRY", "HOXC13", "LEF1", "COL17A1", "MKI67", "SOX21", "KRT35"), reduction = "FA2", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-------Annotating Clusters------#
IRSAnnotations <- c("IRS_Henle", #0
                   "Early_IRS_II", #1
                   "IRS_Cuticle", #2
                   "Early_IRS_I", #3
                   "IRS_Huxley", #4
                   "IRS_Huxley") #5 Late Huxley, but grouping in with the other Huxley cells
#-----Renaming Idents-----#
Idents(Matrix.IRS) <- "scvi_clusters_0.4"
names(IRSAnnotations) <- levels(Matrix.IRS)
Matrix.IRS <- RenameIdents(Matrix.IRS, IRSAnnotations)
Matrix.IRS$IRSAnnotation <- Matrix.IRS@active.ident
Matrix.IRS$IRSAnnotationOrdered <- factor(Matrix.IRS$IRSAnnotation, levels = c("Early_IRS_I", "Early_IRS_II", "IRS_Cuticle", "IRS_Henle", "IRS_Huxley"))
DimPlot(Matrix.IRS, reduction = "UMAP", group.by = "IRSAnnotationOrdered", label = T, repel = T, raster = F)
table(HHA.IRS$IRSAnnotation)

#------Adding Idents to Seurat Object------$
IRS.Metadata <- data.frame(Barcode = colnames(Matrix.IRS), IRSAnnotation = Matrix.IRS$IRSAnnotation)
IRS.Metadata
IRS.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Matrix)), data.frame(Barcode = colnames(Matrix.IRS), IRSAnnotation = Matrix.IRS$IRSAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Matrix$IRSAnnotation <- IRS.Metadata$IRSAnnotation
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "IRSAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#-------Resolving Cuticle Lineages-------#
Matrix.Cuticle <- subset(HHA.Matrix, MatrixAnnotationBroad %in% c("Cuticle"))
Matrix.Cuticle
#-----Performing Final DR and Clustering in scVI Space-----#
Matrix.Cuticle <-FindNeighbors(Matrix.Cuticle, dims = 1:50, reduction = "scvi")
#Matrix.Cuticle <- RunUMAP(Matrix.Cuticle, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 1, cluster.name = "scvi_clusters_1.0")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
Matrix.Cuticle<- FindClusters(Matrix.Cuticle, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(Matrix.Cuticle, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(Matrix.Cuticle, reduction = "FA2", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1)
DimPlot(Matrix.Cuticle, reduction = "UMAP", label = F, group.by = "BulbAnnotationFine", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(Matrix.Cuticle, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting Markers
FeaturePlot(Matrix.Cuticle, features = c("GATA3", "HOXC13", "SOX21", "KRT32", "KRT82", "CRAT"), reduction = "UMAP", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(Matrix.Cuticle, features = c("GATA3", "HOXC13", "SOX21", "KRT32", "KRT82", "CRAT"), reduction = "FA2", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-------Annotating Clusters------#
CuticleAnnotations <- c("Early_Cuticle", #0
                        "Middle_Cuticle", #1
                        "Late_Cuticle") #2
#-----Renaming Idents-----#
Idents(Matrix.Cuticle) <- "scvi_clusters_0.25"
names(CuticleAnnotations) <- levels(Matrix.Cuticle)
Matrix.Cuticle <- RenameIdents(Matrix.Cuticle, CuticleAnnotations)
Matrix.Cuticle$CuticleAnnotation <- Matrix.Cuticle@active.ident
Matrix.Cuticle$CuticleAnnotationOrdered <- factor(Matrix.Cuticle$CuticleAnnotation, levels = c("Early_Cuticle", "Middle_Cuticle", "Late_Cuticle"))
DimPlot(Matrix.Cuticle, reduction = "UMAP", group.by = "CuticleAnnotationOrdered", label = T, repel = T, raster = F)

#------Adding Idents to Seurat Object------$
Cuticle.Metadata <- data.frame(Barcode = colnames(Matrix.Cuticle), CuticleAnnotation = Matrix.Cuticle$CuticleAnnotation)
Cuticle.Metadata
Cuticle.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Matrix)), data.frame(Barcode = colnames(Matrix.Cuticle), CuticleAnnotation = Matrix.Cuticle$CuticleAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Matrix$CuticleAnnotation <- Cuticle.Metadata$CuticleAnnotation
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "CuticleAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)
DimPlot(HHA.Matrix, reduction = "FA2", group.by = "CuticleAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#-------Resolving Cortex Lineages-------#
Matrix.Cortex <- subset(HHA.Matrix, MatrixAnnotationBroad %in% c("Cortex"))
Matrix.Cortex
#-----Performing Final DR and Clustering in scVI Space-----#
Matrix.Cortex <-FindNeighbors(Matrix.Cortex, dims = 1:50, reduction = "scvi")
#Matrix.Cortex <- RunUMAP(Matrix.Cortex, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 1, cluster.name = "scvi_clusters_1.0")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
Matrix.Cortex<- FindClusters(Matrix.Cortex, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
#Plotting Clusters
DimPlot(Matrix.Cortex, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
DimPlot(Matrix.Cortex, reduction = "FA2", label = T, group.by = "scvi_clusters_0.25", raster = F, pt.size = 1)
DimPlot(Matrix.Cortex, reduction = "UMAP", label = F, group.by = "BulbAnnotationFine", raster = F, pt.size = 1) +
  labs(x = "UMAP 1", y = "UMAP 2")
#Plotting QC Metrics
FeaturePlot(Matrix.Cortex, reduction = "UMAP", features = c("DissociationScore"), pt.size = 1) &
  scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting Markers
FeaturePlot(Matrix.Cortex, features = c("ODC1", "HOXC13", "LEF1", "KRT35", "KRT85", "KRT31"), reduction = "UMAP", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
FeaturePlot(Matrix.Cortex, features = c("ODC1", "HOXC13", "LEF1", "KRT35", "KRT85", "KRT31"), reduction = "FA2", ncol = 3, pt.size = 0.25) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-------Annotating Clusters------#
CortexAnnotations <- c("Middle_Cortex", #0
                       "Early_Cortex", #1
                       "Late_Cortex") #2
#-----Renaming Idents-----#
Idents(Matrix.Cortex) <- "scvi_clusters_0.25"
names(CortexAnnotations) <- levels(Matrix.Cortex)
Matrix.Cortex <- RenameIdents(Matrix.Cortex, CortexAnnotations)
Matrix.Cortex$CortexAnnotation <- Matrix.Cortex@active.ident
Matrix.Cortex$CortexAnnotationOrdered <- factor(Matrix.Cortex$CortexAnnotation, levels = c("Early_Cortex", "Middle_Cortex", "Late_Cortex"))
DimPlot(Matrix.Cortex, reduction = "UMAP", group.by = "CortexAnnotationOrdered", label = T, repel = T, raster = F)
#------Adding Idents to Seurat Object------$
Cortex.Metadata <- data.frame(Barcode = colnames(Matrix.Cortex), CortexAnnotation = Matrix.Cortex$CortexAnnotation)
Cortex.Metadata
Cortex.Metadata <- full_join(data.frame(Barcode = colnames(HHA.Matrix)), data.frame(Barcode = colnames(Matrix.Cortex), CortexAnnotation = Matrix.Cortex$CortexAnnotation), by = "Barcode")
#Mapping back to full object
HHA.Matrix$CortexAnnotation <- Cortex.Metadata$CortexAnnotation
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "CortexAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)
DimPlot(HHA.Matrix, reduction = "FA2", group.by = "CortexAnnotation", label = T, repel = T, raster = F, pt.size = 1.0)

#--------Adding Column with All Clusters--------#
stallion = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC",
            "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8",  "#6E4B9E", "#0C727C", "#7E1416", "#D8A767")
HHA.Matrix$MatrixAnnotationFine <- case_when(HHA.Matrix$MatrixAnnotationBroad == "IRS" ~ HHA.Matrix$IRSAnnotation,
                                             HHA.Matrix$MatrixAnnotationBroad %in% c("COL17", "LPC") ~ HHA.Matrix$GLAnnotation,
                                             HHA.Matrix$MatrixAnnotationBroad == "Cortex" ~ HHA.Matrix$CortexAnnotation,
                                             HHA.Matrix$MatrixAnnotationBroad == "Cuticle" ~ HHA.Matrix$CuticleAnnotation,)
#Plotting UMAP
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "MatrixAnnotationFine", label = T, repel = T, raster = F, pt.size = 1.0, cols = stallion) 

#-------Getting Identities for Removed Cells-----------#
#Getting csv with idents for removed cells
Full.Mask <- HHA.Full[[]] %>% filter(GeneralAnnotation == "Bulb")
Full.Mask$Barcode <- rownames(Full.Mask)
#Adding column for cells that were removed. 
Full.Mask$QCFiltered <- Full.Mask$Barcode %in% HHA.Bulb$Barcode == F
table(Full.Mask$QCFiltered)
#Adding column for cells that exist but are in the matrix. 
Full.Mask$InMatrix <- Full.Mask$Barcode %in% HHA.Matrix$Barcode
table(Full.Mask$InMatrix)
write_csv(Full.Mask, paste0(SeuratDirectory, "HHA_Full_Matrix_Mask_5_30_25.csv"))

#--------Mapping Final Annotation and General Annotation onto Matrix Object--------#
Annotation.Metadata <- HHA.Matrix[[]] %>% select(Barcode, MatrixAnnotationFine, MatrixAnnotationBroad)
Annotation.Metadata
#Getting Annotation Columns from original object
Annotation.Metadata <- full_join(Annotation.Metadata, Full.Mask[Full.Mask$InMatrix == T,], by = "Barcode")
rownames(Annotation.Metadata) <- Annotation.Metadata$Barcode
Annotation.Metadata
#Checking order
identical(colnames(HHA.Matrix), rownames(Annotation.Metadata))
HHA.Matrix$FinalAnnotation <- Annotation.Metadata$FinalAnnotation
HHA.Matrix$GeneralAnnotation <- Annotation.Metadata$GeneralAnnotation
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "FinalAnnotation")
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "GeneralAnnotation")

#--------Exporting to Python for Downstream Analysis------#
#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
Matrix.metadata <- HHA.Matrix[[]] %>% mutate(UMAP_1 = HHA.Matrix@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.Matrix@reductions$UMAP@cell.embeddings[,2])
Matrix.metadata$barcode <- rownames(Matrix.metadata)
Matrix.metadata <- Matrix.metadata
write.csv(Matrix.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_Metadata_5_30_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(HHA.Matrix@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_scVI_Embeddings_5_30_25.csv"), quote = F, row.names = F)

#------Exporting SNN Graph Connectivities------#
#Getting long format DF containing all possible connections and their presence/absence
Matrix.SNN <- as.data.frame(summary(HHA.Matrix@graphs$RNA_nn)) %>% #Matching Pythonic Indexing (starting at 0)
  mutate(i = i-1, j = j-1)
write_csv(Matrix.SNN, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_SNN_Connectivities_5_22_25.csv"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Matrix.logcounts <- HHA.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(Matrix.logcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_Expression_LogCounts_5_30_25.mtx"))
#Pulling Raw Counts Data
Matrix.rawcounts <- HHA.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(Matrix.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_RawCounts_5_30_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Matrix.logcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Matrix_FullConversion_GeneNames_5_30_25.csv"),
            quote = F, row.names = F, col.names = T)

#-----Reclustering without cell cycle removal-----#
#-----Getting SCVI Integration Results-----#
Matrix.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Matrix_CC_Included_LatentRep_5_22_25.csv")) %>% column_to_rownames("Barcode")
head(Matrix.scvi)
# put it back in our original Seurat object
HHA.Matrix[["scvi_cc"]] <- CreateDimReducObject(embeddings = as.matrix(Matrix.scvi), key = "scvi_")
rm(Matrix.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Matrix <-FindNeighbors(HHA.Matrix, dims = 1:50, reduction = "scvi_cc")
HHA.Matrix <- RunUMAP(HHA.Matrix, dims = 1:50, reduction = "scvi_cc", reduction.name = "UMAP_CC")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.85, cluster.name = "cc_clusters_0.85")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.8, cluster.name = "cc_clusters_0.8")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.75, cluster.name = "cc_clusters_0.75")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.6, cluster.name = "cc_clusters_0.6")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.5, cluster.name = "cc_clusters_0.5")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.4, cluster.name = "cc_clusters_0.4")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.35, cluster.name = "cc_clusters_0.35")
HHA.Matrix<- FindClusters(HHA.Matrix, resolution = 0.3, cluster.name = "cc_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Matrix, reduction = "UMAP_CC", label = T, group.by = "MatrixAnnotationFine", raster = F, pt.size = 0.5) 
DimPlot(HHA.Matrix, reduction = "UMAP_CC", label = T, group.by = "cc_clusters_0.5", raster = F, pt.size = 0.5) 
DimPlot(HHA.Matrix, reduction = "UMAP_CC", label = T, group.by = "Phase", raster = F, pt.size = 0.5) 

#------Annotating Lineages with Proliferation------#
MatrixAnnotationsCC <- c("Early Cortex", #0
                            "Middle Cortex", #1
                            "IRS", #2
                            "COL17", #3
                            "Cuticle", #4
                            "Proliferating Cortex/Cuticle", #5
                            "Late Cortex", #6
                            "Proliferating IRS", #7
                            "LPC") #8
#-----Renaming Idents-----#
Idents(HHA.Matrix) <- "cc_clusters_0.5"
names(MatrixAnnotationsCC) <- levels(HHA.Matrix)
HHA.Matrix <- RenameIdents(HHA.Matrix, MatrixAnnotationsCC)
HHA.Matrix$MatrixAnnotationCC<- HHA.Matrix@active.ident
DimPlot(HHA.Matrix, reduction = "UMAP_CC", group.by = "MatrixAnnotationCC", label = T, repel = T, raster = F, pt.size = 1.0)

#--------Plotting Proliferative Clusters on Gray Background---------#
HHA.Matrix$ProliferativeGrouping <- case_when(HHA.Matrix$MatrixAnnotationCC == "Proliferating Cortex/Cuticle" ~ "Proliferating Cortex/Cuticle",
                                              HHA.Matrix$MatrixAnnotationCC == "Proliferating IRS" ~ "Proliferating IRS",
                                              .default = "Other")
DimPlot(HHA.Matrix, reduction = "UMAP_CC", group.by = "ProliferativeGrouping", label = T, repel = T, raster = F, pt.size = 1.0, cols = c("gray", "firebrick",  "cornflowerblue"))
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "ProliferativeGrouping", label = T, repel = T, raster = F, pt.size = 1.0, cols = c("gray", "firebrick",  "cornflowerblue"))

#---------Plotting UMAPs as TIFFs-----#
#-------Bulb Proliferative Clusters: CC UMAP, Gray Background Clusters-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_CC_Clusters_Highlighted_CC_UMAP.png"), res = 300, height = 6, width = 8, units = "in") 
DimPlot(HHA.Matrix, reduction = "UMAP_CC", group.by = "ProliferativeGrouping", raster = F, pt.size = 0.5, cols = c("gray", "firebrick",  "cornflowerblue")) &
  NoAxes() & theme(aspect.ratio = 1) & labs(title = "")
dev.off()
#-------Bulb Proliferative Clusters: CC Regressed UMAP, Gray Background Clusters-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_CC_Clusters_Highlighted_UMAP.png"), res = 300, height = 6, width = 8, units = "in") 
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "ProliferativeGrouping", raster = F, pt.size = 0.5, cols = c("gray", "firebrick",  "cornflowerblue")) &
  NoAxes() & theme(aspect.ratio = 1) & labs(title = "")
dev.off()

#----------Plotting Clusters-------#
MatrixPal <- c("#89C75F", "#3BBCA8",  "#208A42", "#0C727C", "#9ECAE1", "#4292C6", "#08306B", "#E6C2DC", "#C06CAB",  "#89288F", "#D8A767", "#F47D2B", "#F37B7D",  "#7E1416", "#D24B27")
HHA.Matrix$MatrixAnnotationOrdered <- factor(HHA.Matrix$MatrixAnnotationFine, levels = c("Lower COL17", "Upper COL17", "LPC", "Medulla", 
                                                                                         "Early_Cortex", "Middle_Cortex", "Late_Cortex", 
                                                                                         "Early_Cuticle", "Middle_Cuticle", "Late_Cuticle",
                                                                                         "Early_IRS_I", "Early_IRS_II", "IRS_Henle", "IRS_Huxley", "IRS_Cuticle"))
#Generating Plot
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_Clusters_UMAP.png"), res = 300, height = 6, width = 8, units = "in") 
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "MatrixAnnotationOrdered", raster = F, pt.size = 0.25, cols = MatrixPal) &
  NoAxes() & theme(aspect.ratio = 1) & labs(title = "")
dev.off()

#-------Plotting COL17A1 Expression on UMAP-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_COL17A1_UMAP.png"), res = 300, height = 6, width = 8, units = "in") 
FeaturePlot(HHA.Matrix,  features = 'COL17A1', reduction = "UMAP", raster = F, pt.size = 0.5) &
  NoAxes() & theme(aspect.ratio = 1) & scale_color_viridis_c(option = "magma")
dev.off()

#---------Saving Object and Workspace--------#
#Saving Object
save(HHA.Matrix, file = paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))
#Saving Workspace
save.image(paste0(WorkspaceDirectory, "HHA_Bulb_Recluster_5_22_25.RData"))
