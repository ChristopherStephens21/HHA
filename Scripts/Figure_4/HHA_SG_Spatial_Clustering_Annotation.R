#-------Clustering and Integration of HHA scRNA-seq Data-----#
#This is our workflow for analysis of our Visium HD spatial transcriptomics data
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)

#-----Loading in Packages----#
library(Seurat) #scRNA-seq
library(EnsDb.Hsapiens.v86) #Genome for ClusterProfiler
library(clusterProfiler) #Gene set overrepresentation analysis 
library(tidyverse) #Data handling
library(scCustomize) #Reading in CellBender output
library(Signac) #trying out 
#library(Banksy) #Spatial Clustering 

#----Directories-----#
SeuratDirectory <- "/projects/b1217/HHA/Spatial/Seurat_Data/"
PlotsDirectory <- "/projects/b1217/HHA/Spatial/Plots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
AnnDataDirectory <- "/projects/b1217/HHA/SG_Trajectory/AnnData/"

#------Palettes-----#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
#From MetBrewer::met.brewer("Signac", 10) "#92c051" "#2b9b81" "#1f6e9c" "#633372" "#9f5691" "#e87b89" "#de597c" "#d8443c" "#fe9b00" "#f4c40f"
SampleIDPal <- c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f")
Metadata.Pal <- scale_fill_manual(values = SampleIDPal)
solarExtra = c('#3361A5','#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D') 

#-------Loading Checkpoint-----#
load(paste0(WorkspaceDirectory, "HHA_Spatial_UpperHF_SG_Subclustering_4_30_25.RData"))

#------Loading Full Spatial Object-------#
load(paste0(SeuratDirectory, "HHA_Integrated_Spatial_Initial_Annotations_12_11_24.rds"))

#------Plotting UMAP-----#
#Initial Annotation
DimPlot(HHA.Merged, reduction = "full.umap.sketch", label = T, repel = T, raster = T, group.by = "InitialAnnotationFine") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()
#Broad Annotation
DimPlot(HHA.Merged, reduction = "full.umap.sketch", label = T, repel = T, raster = T, group.by = "InitialAnnotationBroad") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()

#------Subsetting for Populations within the Upper HF and SG------#
Spatial.Upper <- subset(HHA.Merged, InitialAnnotationFine %in% c("Sebaceous I", "Sebaceous II", "Sebaceous III", "Sebaceous IV", "Sebaceous Duct", 
                                                                 "Basal Infundibulum", "Suprabasal Infundibulum"))
Spatial.Upper #86,523 spots
#Sanity checking selection
DimPlot(Spatial.Upper, reduction = "full.umap.sketch", label = T, repel = T, raster = T, group.by = "InitialAnnotationFine") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()

#-----Performing clustering and DR of 8 um data-----#
#Finding variable features and scaling data
DefaultAssay(Spatial.Upper) <- "Spatial.008um"
Spatial.Upper <- FindVariableFeatures(Spatial.Upper, selection.method = "vst", nfeatures = 2000)
Spatial.Upper <- ScaleData(Spatial.Upper)
Spatial.Upper <- RunPCA(Spatial.Upper, assay =  "Spatial.008um", reduction.name = "pca")
Spatial.Upper <- FindNeighbors(Spatial.Upper, assay =  "Spatial.008um", reduction = "pca", dims = 1:10)
Spatial.Upper <- RunUMAP(Spatial.Upper, reduction = "pca", reduction.name = "UMAP_Unintegrated", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "unintegrated_clusters", resolution = 0.6)

#--------Plotting Initial Clustering------#
#Visualizing Results: Clusters
DimPlot(Spatial.Upper, reduction = "UMAP_Unintegrated", label = T, repel = T, raster = F, group.by = "unintegrated_clusters") + 
  labs(x = "UMAP 1", y = "UMAP 2")
#Visualizing Results: Samples
DimPlot(Spatial.Upper, reduction = "UMAP_Unintegrated", raster = F, group.by = "orig.ident") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)")
#Visualizing Results: Initial Annotation
DimPlot(Spatial.Upper, reduction = "UMAP_Unintegrated", label = T, repel = T, raster = F, group.by = "InitialAnnotationFine") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()

#-----Performing Sketch integration with harmony------#
set.seed(1701)
Spatial.Upper <- IntegrateLayers(Spatial.Upper, method = HarmonyIntegration, assay = "Spatial.008um", orig = "pca", new.reduction = "harmony",
                                 dims = 1:10)
#Reclustering and DR in Harmony Space
Spatial.Upper <- FindNeighbors(Spatial.Upper, assay = "Spatial.008um", reduction = "harmony", dims = 1:10)
Spatial.Upper <- RunUMAP(Spatial.Upper, reduction = "harmony", reduction.name = "UMAP", return.model = T, dims = 1:10)
#Clustering cells using Louvain algorithm
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.3", resolution = 0.3)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.5", resolution = 0.5)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.6", resolution = 0.5)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.75", resolution = 0.75)

#--------Plotting Initial Clustering------#
#Visualizing Results: Clusters
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, group.by = "harmony_clusters_0.5", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")
#Visualizing Results: Samples
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, repel = T, raster = F, group.by = "orig.ident") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()
#Visualizing Results: Initial Annotaton
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, repel = T, raster = F, group.by = "InitialAnnotationFine") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()

#------Plotting Markers-----#
#Sebaceous Gland
FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("KRT5", "AWAT1", "MGST1", "KRT79", "PPARG", "KRT6A")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Dermal Sheath
FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("VIM", "COL1A1", "COL6A1", "SPARC")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#-----Looking at QC Metrics by Cluster-----#
VlnPlot(Spatial.Upper, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"),
        group.by = "harmony_clusters_0.5", alpha = 0.1)
#-----Looking at Counts by Cluster-----#
VlnPlot(Spatial.Upper, features = c("nCount_Spatial.008um", "nFeature_Spatial.008um"),
        group.by = "harmony_clusters_0.5", split.by = "orig.ident", alpha = 0.01)
#Examining Sample Ratios by Cluster
ggplot(Spatial.Upper[[]], aes(x = harmony_clusters_0.5, fill = orig.ident)) + geom_bar(position = "fill")

#------Adding Log-Transformed Count and Feature Columns-------#
Spatial.Upper$LogFeature <- log1p(Spatial.Upper$nFeature_Spatial.008um)
Spatial.Upper$LogCount <- log1p(Spatial.Upper$nCount_Spatial.008um)
FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("LogFeature", "LogCount")) &
  scale_color_viridis_c() & NoAxes()
#-----Looking at QC Metrics by Cluster-----#
VlnPlot(Spatial.Upper, features = c("LogFeature", "LogCount"),
        group.by = "harmony_clusters_0.5", alpha = 0.1)
#-----Looking at QC Metrics by Cluster-----#
VlnPlot(Spatial.Upper, features = c("LogFeature", "LogCount"),
        group.by = "harmony_clusters_0.5", split.by = "orig.ident", alpha = 0.01)

#-----Performing Initial QC and Reclustering-----#
Spatial.Upper$QCThreshold <- case_when(Spatial.Upper$nFeature_Spatial.008um > 150 ~ "Above Threshold",
                                       Spatial.Upper$nFeature_Spatial.008um <= 150 ~ "Below Threshold")
#Plotting Thresholds: By Cluster
ggplot(Spatial.Upper[[]], aes(x = harmony_clusters_0.5, y = nFeature_Spatial.008um, fill = harmony_clusters_0.5)) + 
  geom_violin() + theme_bw() + geom_hline(yintercept = 150, linetype = "dashed") + theme_bw()
#Plotting Thresholds: By Cluster and Sample
ggplot(Spatial.Upper[[]], aes(x = harmony_clusters_0.5, y = nFeature_Spatial.008um, fill = orig.ident)) + 
  geom_violin() + geom_hline(yintercept = 150) + theme_bw()
#Plotting Thresholds: On UMAP
DimPlot(Spatial.Upper, reduction = "UMAP", raster = F, group.by = "QCThreshold") 

#------Removing Low Feature Cells and Reclustering-----#
#Saving Object
#save(Spatial.Upper, file = paste0(SeuratDirectory, "Spatial_UpperHF_SG_Pre_Filterng_4_13_25.rds"))
#Filtering at Threshold
Spatial.Upper <- subset(Spatial.Upper, QCThreshold == "Above Threshold")
Spatial.Upper #73,861 spots
#Visualizing Results: Clusters
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, group.by = "harmony_clusters_0.5", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")

#-----Reclustering after Removal of Poor Quality Cells------#
#Finding variable features and scaling data
DefaultAssay(Spatial.Upper) <- "Spatial.008um"
Spatial.Upper <- ScaleData(Spatial.Upper)
Spatial.Upper <- RunPCA(Spatial.Upper, assay =  "Spatial.008um", reduction.name = "pca")
Spatial.Upper <- FindNeighbors(Spatial.Upper, assay =  "Spatial.008um", reduction = "pca", dims = 1:10)
Spatial.Upper <- RunUMAP(Spatial.Upper, reduction = "pca", reduction.name = "UMAP_Unintegrated", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "unintegrated_clusters", resolution = 0.6)

#-----Performing Sketch integration with harmony------#
set.seed(1701)
Spatial.Upper <- IntegrateLayers(Spatial.Upper, method = HarmonyIntegration, assay = "Spatial.008um", orig = "pca", new.reduction = "harmony",
                                 dims = 1:10)
#Reclustering and DR in Harmony Space
Spatial.Upper <- FindNeighbors(Spatial.Upper, assay = "Spatial.008um", reduction = "harmony", dims = 1:10)
Spatial.Upper <- RunUMAP(Spatial.Upper, reduction = "harmony", reduction.name = "UMAP", return.model = T, dims = 1:10)
#Clustering cells using Louvain algorithm
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.35", resolution = 0.35)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.3", resolution = 0.3)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.4", resolution = 0.4)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.5", resolution = 0.5)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.6", resolution = 0.6)
Spatial.Upper <- FindClusters(Spatial.Upper, cluster.name = "harmony_clusters_0.75", resolution = 0.75)

#--------Plotting Initial Clustering------#
#Visualizing Results: Clusters
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, group.by = "harmony_clusters_0.3", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")

#------Plotting Markers-----#
#Sebaceous Gland
FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("KRT5", "AWAT1", "MGST1", "KRT79", "PPARG", "KRT6A", "KRT1", "VIM")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Dermal Sheath
FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("VIM", "COL1A1", "COL6A1", "SPARC")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#------Plotting Spatial DimPlot of Annotated Clusters-----#
set.seed(1701)
UpperPal <- DiscretePalette(12)

ragg::agg_tiff(filename = paste0(PlotsDirectory, "Upper_SG_SpatialDimPlot_4_28_25.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(Spatial.Upper, group.by = "harmony_clusters_0.3", images =  "Donor1_1.008um", label = F) &
  scale_fill_manual(values = UpperPal)
dev.off()

#-------Annotating Clusters-----#
UpperAnnotation <- c(
  "Basal Infundibulum", #0
  "Sebaceous III", #1
  "Sebaceous II", #2
  "Sebaceous I", #3
  "Transitional Basal", #4
  "Suprabasal Infundibulum I", #5
  "Suprabasal Duct", #6
  "Suprabasal Infundibulum II", #7,
  "Basal Infundibulum") #8
#Renaming Idents
Idents(Spatial.Upper) <- Spatial.Upper$harmony_clusters_0.4
names(UpperAnnotation) <- sort(as.numeric(levels(Spatial.Upper)))
Spatial.Upper <- RenameIdents(Spatial.Upper, UpperAnnotation)
Spatial.Upper$UpperAnnotation <- Spatial.Upper@active.ident

#-------Plotting UMAP------#
DimPlot(Spatial.Upper, reduction = "UMAP", label = T, repel = T, group.by = "UpperAnnotation", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2") & theme(aspect.ratio = 1)

#------Plotting Spatial DimPlot of Annotated Clusters-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Upper_SG_SpatialDimPlot_4_28_25.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(Spatial.UpperSG, group.by = "UpperSGAnnotation", images =  "Donor3.008um", label = F) &
  scale_fill_manual(values = c("#30123BFF", "#4777EFFF", "#1BD0D5FF", "#62FC6BFF", 
                               "#D2E935FF", "#FE9B2DFF", "#DB3A07FF", "#7A0403FF"))
dev.off()

FeaturePlot(Spatial.Upper, reduction = 'UMAP', features = c("KRT1", "MGST1", "KRT79", "KRT6A")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#------Full Conversion to Scanpy------#
#----Writing Metadata----#
Upper.metadata <- Spatial.Upper[[]] %>% mutate(UMAP_1 = Spatial.Upper@reductions$UMAP@cell.embeddings[,1], UMAP_2 = Spatial.Upper@reductions$UMAP@cell.embeddings[,2])
Upper.metadata$barcode <- rownames(Upper.metadata)
Upper.metadata <- Upper.metadata
write.csv(Upper.metadata, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_Metadata_4_28_25.csv"), quote = F, row.names = F)
#----Writing Harmony Embeddings----#
write.csv(Spatial.Upper@reductions$harmony@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_Harmony_Embeddings_4_28_25.csv"), quote = F, row.names = F)
#----Writing Adjacency Matrices (Graphs)---#
#Nearest Neighbor Graph
RNAnn.adjmat <- as(object = Spatial.Upper@graphs$Spatial.008um_nn, Class = "Matrix") #Converting to Matrix
RNAnn.adjmat <- as(object = RNAnn.adjmat, Class = "dgCMatrix") #Converting to Sparse
str(RNAnn.adjmat)
Matrix::writeMM(RNAnn.adjmat, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_NNAdjMat_4_28_25.mtx"))
#Shared Nearest Neighbor Graph
RNAsnn.adjmat <- as(object = Spatial.Upper@graphs$Spatial.008um_snn, Class = "Matrix") #Conveerting to matrix
RNAsnn.adjmat <- as(object = RNAsnn.adjmat, Class = "dgCMatrix") #Converting to Matrix
str(RNAsnn.adjmat)
Matrix::writeMM(RNAsnn.adjmat, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_SNNAdjMat_4_28_25.mtx"))
#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Upper.logcounts <- Spatial.Upper %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'Spatial.008um', layer = 'data')
Matrix::writeMM(Upper.logcounts, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_Expression_LogCounts_4_28_25.mtx"))
#Pulling Log Transformed, Unnormalized Counts
Upper.logcounts_unnorm <- Spatial.Upper %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'Log_008um', layer = 'data')
Matrix::writeMM(Upper.logcounts_unnorm, file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_Expression_LogCounts_Unnormalized_4_28_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Upper.logcounts)), file = paste0(AnnDataDirectory, "HHA_Spatial_UpperHF_SG_Conversion_GeneNames_4_28_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing intermediate objects
rm(Upper.metadata, RNAnn.adjmat, RNAsnn.adjmat, Upper.logcounts, Upper.logcounts_unnorm)

#------Subsetting for Sebaceous Gland and Basal Cells-----#
Spatial.SG <- subset(Spatial.Upper, UpperAnnotation %in% c("Basal Infundibulum", "Transitional Basal", "Sebaceous I", "Sebaceous II",
                                                   "Sebaceous III", "Suprabasal Duct"))
Spatial.SG #61,900 spots
#-------Plotting UMAP------#
DimPlot(Spatial.SG, reduction = "UMAP", label = T, repel = T, group.by = "UpperAnnotation", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2") & theme(aspect.ratio = 1)

#-----Reclustering Cells-----#
#Finding variable features and scaling data
DefaultAssay(Spatial.SG) <- "Spatial.008um"
Spatial.SG <- FindVariableFeatures(Spatial.SG, selection.method = "vst", nfeatures = 2000)
Spatial.SG <- ScaleData(Spatial.SG)
Spatial.SG <- RunPCA(Spatial.SG, assay =  "Spatial.008um", reduction.name = "pca")
Spatial.SG <- FindNeighbors(Spatial.SG, assay =  "Spatial.008um", reduction = "pca", dims = 1:10)
Spatial.SG <- RunUMAP(Spatial.SG, reduction = "pca", reduction.name = "UMAP_Unintegrated", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "unintegrated_clusters", resolution = 0.6)
#-----Performing Sketch integration with harmony------#
set.seed(1701)
Spatial.SG <- IntegrateLayers(Spatial.SG, method = HarmonyIntegration, assay = "Spatial.008um", orig = "pca", new.reduction = "harmony",
                              dims = 1:10)
#Reclustering and DR in Harmony Space
Spatial.SG <- FindNeighbors(Spatial.SG, assay = "Spatial.008um", reduction = "harmony", dims = 1:10)
Spatial.SG <- RunUMAP(Spatial.SG, reduction = "harmony", reduction.name = "UMAP", return.model = T, dims = 1:10)
#Clustering cells using Louvain algorithm
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.3", resolution = 0.3)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.5", resolution = 0.5)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.6", resolution = 0.5)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.75", resolution = 0.75)

#--------Plotting Initial Clustering------#
#Visualizing Results: Clusters
DimPlot(Spatial.SG, reduction = "UMAP", label = T, group.by = "harmony_clusters_0.3", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")
#Visualizing Results: Clusters
DimPlot(Spatial.SG, reduction = "UMAP", label = T, group.by = "UpperAnnotation", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")

#------Plotting Markers-----#
#Sebaceous Gland
FeaturePlot(Spatial.SG, reduction = 'UMAP', features = c("KRT5", "AWAT1", "MGST1", "KRT79", "PPARG", "KRT6A", "KRT1", "VIM", "IVL")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Dermal Sheath
FeaturePlot(Spatial.SG, reduction = 'UMAP', features = c("VIM", "COL1A1", "COL6A1", "SPARC")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Cluster 7 is residual infundibulum 
ragg::agg_tiff(filename = paste0(PlotsDirectory, "SG_SpatialDimPlot_4_28_25.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(Spatial.SG, group.by = "harmony_clusters_0.75", images =  "Donor1_1.008um", label = F) 
dev.off()

#------Subsetting for Sebaceous Gland and Basal Cells-----#
Spatial.SG <- subset(Spatial.SG, harmony_clusters_0.3 != 7)
Spatial.SG #61,900 spots

#-----Reclustering Cells-----#
#Finding variable features and scaling data
DefaultAssay(Spatial.SG) <- "Spatial.008um"
Spatial.SG <- FindVariableFeatures(Spatial.SG, selection.method = "vst", nfeatures = 2000)
Spatial.SG <- ScaleData(Spatial.SG)
Spatial.SG <- RunPCA(Spatial.SG, assay =  "Spatial.008um", reduction.name = "pca")
Spatial.SG <- FindNeighbors(Spatial.SG, assay =  "Spatial.008um", reduction = "pca", dims = 1:10)
Spatial.SG <- RunUMAP(Spatial.SG, reduction = "pca", reduction.name = "UMAP_Unintegrated", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "unintegrated_clusters", resolution = 0.6)
#-----Performing Sketch integration with harmony------#
set.seed(1701)
Spatial.SG <- IntegrateLayers(Spatial.SG, method = HarmonyIntegration, assay = "Spatial.008um", orig = "pca", new.reduction = "harmony",
                              dims = 1:10)
#Reclustering and DR in Harmony Space
Spatial.SG <- FindNeighbors(Spatial.SG, assay = "Spatial.008um", reduction = "harmony", dims = 1:10)
Spatial.SG <- RunUMAP(Spatial.SG, reduction = "harmony", reduction.name = "UMAP", return.model = T, dims = 1:10)
#Clustering cells using Louvain algorithm
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.25", resolution = 0.25)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.3", resolution = 0.3)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.35", resolution = 0.35)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.4", resolution = 0.4)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.5", resolution = 0.5)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.6", resolution = 0.5)
Spatial.SG <- FindClusters(Spatial.SG, cluster.name = "harmony_clusters_0.75", resolution = 0.75)

#--------Plotting Initial Clustering------#
#Visualizing Results: Clusters
DimPlot(Spatial.SG, reduction = "UMAP", label = T, group.by = "harmony_clusters_0.35", raster = F) + 
  labs(x = "UMAP 1", y = "UMAP 2")

#------Plotting Markers-----#
#Sebaceous Gland
FeaturePlot(Spatial.SG, reduction = 'UMAP', features = c("KRT5", "AWAT1", "MGST1", "KRT79", "PPARG", "KRT6A", "KRT1", "VIM", "IVL")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Dermal Sheath
FeaturePlot(Spatial.SG, reduction = 'UMAP', features = c("VIM", "COL1A1", "COL6A1", "SPARC")) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Cluster 7 is residual infundibulum 
ragg::agg_tiff(filename = paste0(PlotsDirectory, "SG_SpatialDimPlot_4_28_25.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(Spatial.SG, group.by = "harmony_clusters_0.3", images =  "Donor1_1.008um", label = F) 
dev.off()

#-------Annotating Clusters-----#
SGAnnotation <- c(
  "Maturation I", #0
  "Maturation III", #1
  "Basal Duct", #2
  "Peripheral Zone", #3
  "Maturation II", #4
  "Sebaceous Duct", #5
  "Differentiated Zone") #6
#Renaming Idents
Idents(Spatial.SG) <- Spatial.SG$harmony_clusters_0.35
names(SGAnnotation) <- sort(as.numeric(levels(Spatial.SG)))
Spatial.SG <- RenameIdents(Spatial.SG, SGAnnotation)
Spatial.SG$SGAnnotation <- Spatial.SG@active.ident
#Factoring cluster labels
Spatial.SG$SGAnnotationOrdered <- factor(Spatial.SG$SGAnnotation, levels = c("Basal Duct", "Peripheral Zone",
                                                                      "Maturation I", "Maturation II", "Maturation III", "Differentiated Zone",
                                                                      "Sebaceous Duct"))
#--------Plotting Clustering------#
#ArchR Summer Night
SGPal <- c("#9cdff0", "#2a7185", "#a64027", "#fbdf72","#60824f", "#022336", "#725ca5")
#Visualizing Results: Clusters
DimPlot(Spatial.SG, reduction = "UMAP", label = F, group.by = "SGAnnotationOrdered", cols = SGPal) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "")  + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)

#------Loading in Pseudotime Results-----#
#Reading palantir results
SG.Pseudotime <- read_csv(paste0(AnnDataDirectory, "HHA_SG_palantir_pseudotime_results_4_30_25.csv"))
#Adding to SeuratObject
Spatial.SG$PalantirPseudotime <- SG.Pseudotime$palantir_pseudotime

#--------Plotting Pseudotime-------#
#Plotting on UMAP
FeaturePlot(Spatial.SG, reduction = "UMAP", features = "PalantirPseudotime") & scale_color_viridis_c(option = "magma") & 
  NoAxes() & theme(aspect.ratio = 1)
#Plotting by cluster
Spatial.SG[[]] %>% filter(SGAnnotation != "Basal Duct") %>% 
ggplot(aes(x = PalantirPseudotime, y = SGAnnotationOrdered, color = SGAnnotationOrdered)) + ggbeeswarm::geom_quasirandom(size = 0.01) + 
  theme_bw() + scale_color_manual(values = SGPal) +
  labs(y = "Annotation", color = "Annotation", x = "Palantir Pseudotime") + theme(aspect.ratio = 1)
#Plotting by cluster
Spatial.SG[[]] %>% filter(SGAnnotation != "Basal Duct") %>% 
  ggplot(aes(x = PalantirPseudotime, y = SGAnnotationOrdered, color = PalantirPseudotime)) + ggbeeswarm::geom_quasirandom(size = 0.01) + 
  theme_bw() + scale_color_viridis_c(option = "magma") +
  labs(y = "Annotation", x = "Palantir Pseudotime") + theme(aspect.ratio = 1)

#------Full Conversion to Scanpy------#
#----Writing Metadata----#
SG.metadata <- Spatial.SG[[]] %>% mutate(UMAP_1 = Spatial.SG@reductions$UMAP@cell.embeddings[,1], UMAP_2 = Spatial.SG@reductions$UMAP@cell.embeddings[,2])
SG.metadata$barcode <- rownames(SG.metadata)
SG.metadata <- SG.metadata
write.csv(SG.metadata, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_Metadata_4_30_25.csv"), quote = F, row.names = F)
#----Writing Harmony Embeddings----#
write.csv(Spatial.SG@reductions$harmony@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_Harmony_Embeddings_4_30_25.csv"), quote = F, row.names = F)
#----Writing Adjacency Matrices (Graphs)---#
#Nearest Neighbor Graph
RNAnn.adjmat <- as(object = Spatial.SG@graphs$Spatial.008um_nn, Class = "Matrix") #Converting to Matrix
RNAnn.adjmat <- as(object = RNAnn.adjmat, Class = "dgCMatrix") #Converting to Sparse
str(RNAnn.adjmat)
Matrix::writeMM(RNAnn.adjmat, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_NNAdjMat_4_30_25.mtx"))
#Shared Nearest Neighbor Graph
RNAsnn.adjmat <- as(object = Spatial.SG@graphs$Spatial.008um_snn, Class = "Matrix") #Conveerting to matrix
RNAsnn.adjmat <- as(object = RNAsnn.adjmat, Class = "dgCMatrix") #Converting to Matrix
str(RNAsnn.adjmat)
Matrix::writeMM(RNAsnn.adjmat, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_SNNAdjMat_4_30_25.mtx"))
#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
SG.logcounts <- Spatial.SG %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'Spatial.008um', layer = 'data')
Matrix::writeMM(SG.logcounts, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_Expression_LogCounts_4_30_25.mtx"))
#Pulling Log Transformed, Unnormalized Counts
SG.logcounts_unnorm <- Spatial.SG %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'Log_008um', layer = 'data')
Matrix::writeMM(SG.logcounts_unnorm, file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_Expression_LogCounts_Unnormalized_4_30_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(SG.logcounts)), file = paste0(AnnDataDirectory, "HHA_Spatial_SG_Conversion_GeneNames_4_30_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing intermediate objects
rm(SG.metadata, RNAnn.adjmat, RNAsnn.adjmat, SG.logcounts, SG.logcounts_unnorm)

#--------Saving Workspace Image-------#
save.image(paste0(WorkspaceDirectory, "HHA_Spatial_UpperHF_SG_Subclustering_4_30_25.RData"))
#Saving Seurat Objects
save(Spatial.Upper, file = paste0(SeuratDirectory, "HHA_Spatial_UpperHF_SG_4_30_25.rds"))
#Saving Seurat Objects
save(Spatial.SG, file = paste0(SeuratDirectory, "HHA_Spatial_SG_4_30_25.rds"))


