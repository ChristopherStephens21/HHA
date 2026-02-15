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

#-----Loading Final Object-----#
load(paste0(SeuratDirectory, "HHA_Integrated_Spatial_Initial_Annotations_12_11_24.rds"))

#------Palettes-------#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
solarExtra = c('#3361A5','#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D') 

#------Helper Functions----------#
#Reformats DefaultAssay() to tidy format for easy piping/use with purr::map()
DefaultAssayTidy <- function(SeuratObj, assay) {
  DefaultAssay(SeuratObj) <- assay
  return(SeuratObj)}

#-----Loading in data-----#
#Donor 1 Slice 1
EL_S1 <- Load10X_Spatial("/projects/b1217/Chris/Spatial/20240612EL_HuScalp/outs/",
                         filename = "filtered_feature_bc_matrix.h5",
                         slice = "Donor1_1",
                         bin.size = 8)
EL_S1$orig.ident <- "EL_S1"
#Donor 1 Slice 2
EL_S2 <- Load10X_Spatial("/projects/b1217/Chris/Spatial/20240612EL_HuScalp_V2_SpatialData/outs/",
                         filename = "filtered_feature_bc_matrix.h5",
                         slice = "Donor1_2",
                         bin.size = 8)
EL_S2$orig.ident <- "EL_S2"
#Donor 2
EL_S3 <- Load10X_Spatial("/projects/b1217/Chris/Spatial/20241002MK_HuScalp/outs/",
                         filename = "filtered_feature_bc_matrix.h5",
                         slice = "Donor2",
                         bin.size = 8)
EL_S3$orig.ident <- "EL_S3"
#Donor 3
EL_S4 <- Load10X_Spatial("/projects/b1217/Chris/Spatial/20240610LG_HuScalp/outs/",
                         filename = "filtered_feature_bc_matrix.h5",
                         slice = "Donor3",
                         bin.size = 8)
EL_S4$orig.ident <- "EL_S4"

#-----Log Normalizing Data------#
#---EL_S1---#
EL_S1 <- EL_S1 %>% DefaultAssayTidy("Spatial.008um") %>% NormalizeData()
#Adding Log Transformed (but otherwise unnormalized) Assay
EL_S1[["Log_008um"]] <- CreateAssayObject(counts = log1p(GetAssayData(EL_S1, assay = "Spatial.008um", slot = "counts")))
#---EL_S2---#
EL_S2 <- EL_S2 %>% DefaultAssayTidy("Spatial.008um") %>% NormalizeData()
#Adding Log Transformed (but otherwise unnormalized) Assay
EL_S2[["Log_008um"]] <- CreateAssayObject(counts = log1p(GetAssayData(EL_S2, assay = "Spatial.008um", slot = "counts")))
#---EL_S3---#
EL_S3 <- EL_S3 %>% DefaultAssayTidy("Spatial.008um") %>% NormalizeData()
#Adding Log Transformed (but otherwise unnormalized) Assay
EL_S3[["Log_008um"]] <- CreateAssayObject(counts = log1p(GetAssayData(EL_S3, assay = "Spatial.008um", slot = "counts")))
#---EL_S2---#
EL_S4 <- EL_S4 %>% DefaultAssayTidy("Spatial.008um") %>% NormalizeData()
#Adding Log Transformed (but otherwise unnormalized) Assay
EL_S4[["Log_008um"]] <- CreateAssayObject(counts = log1p(GetAssayData(EL_S4, assay = "Spatial.008um", slot = "counts")))

#------Plotting Counts by Spot--------#
#EL_S1
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S1_SpatialCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S1 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nCount_Spatial.008um", shape = 22, images = c("Donor1_1.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#EL_S2
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S2_SpatialCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S2 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nCount_Spatial.008um", shape = 22, images = c("Donor1_2.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#EL_S3
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S3_SpatialCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S3 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nCount_Spatial.008um", shape = 22, images = c("Donor2.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#EL_S4
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S4_SpatialCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S4 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nCount_Spatial.008um", shape = 22, images = c("Donor3.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()

#------Plotting Gene Counts by Spot--------#
#EL_S1
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S1_SpatialFeaturePlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S1 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nFeature_Spatial.008um", shape = 22, images = c("Donor1_1.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#_EL_S2
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S2_SpatialFeaturePlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S2 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nFeature_Spatial.008um", shape = 22, images = c("Donor1_2.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#EL_S3
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S3_SpatialFeaturePlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S3 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nFeature_Spatial.008um", shape = 22, images = c("Donor2.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()
#_EL_S4
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S4_SpatialFeaturePlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
EL_S4 %>% DefaultAssayTidy("Spatial.008um") %>% 
  SpatialFeaturePlot(features = "nFeature_Spatial.008um", shape = 22, images = c("Donor3.008um")) & 
  scale_color_viridis_c(option = "turbo") & theme(legend.position = "right")
dev.off()

#-----Plotting Regions with Low Counts-----#
#EL_S1
EL_S1$LowCounts <- EL_S1$nCount_Spatial.008um < 10
table(EL_S1$LowCounts)
#EL_S2
EL_S2$LowCounts <- EL_S2$nCount_Spatial.008um < 10
table(EL_S2$LowCounts)
#EL_S1
EL_S3$LowCounts <- EL_S3$nCount_Spatial.008um < 10
table(EL_S3$LowCounts)
#EL_S2
EL_S4$LowCounts <- EL_S4$nCount_Spatial.008um < 10
table(EL_S4$LowCounts)

#-------Plotting Low versus Non-Low Count Regions-----#
#EL_S1
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S1_LowCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(EL_S1, group.by = "LowCounts", shape = 22, alpha = 0.25)
dev.off()
#EL_S2
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S2_LowCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(EL_S2, group.by = "LowCounts", shape = 22, alpha = 0.25)
dev.off()
#EL_S3
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S3_LowCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(EL_S3, group.by = "LowCounts", shape = 22, alpha = 0.25)
dev.off()
#EL_S4
ragg::agg_tiff(filename = paste0(PlotsDirectory, "EL_S4_LowCountPlot_8um_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(EL_S4, group.by = "LowCounts", shape = 22, alpha = 0.25)
dev.off()

#-----Removing Low Count Regions for Better Clustering-----#
#Saving original aspect ratios
EL_S1_bbox.original <- EL_S1@images$Donor1_1.008um$centroids@bbox
EL_S1_ratio.original <- (EL_S1_bbox.original[1,2] - EL_S1_bbox.original[1,1]) / (EL_S1_bbox.original[2,2] - EL_S1_bbox.original[2,1])
EL_S2_bbox.original <- EL_S2@images$Donor1_2.008um$centroids@bbox
EL_S2_ratio.original <- (EL_S2_bbox.original[1,2] - EL_S2_bbox.original[1,1]) / (EL_S2_bbox.original[2,2] - EL_S2_bbox.original[2,1])
EL_S3_bbox.original <- EL_S3@images$Donor2.008um$centroids@bbox
EL_S3_ratio.original <- (EL_S3_bbox.original[1,2] - EL_S3_bbox.original[1,1]) / (EL_S3_bbox.original[2,2] - EL_S3_bbox.original[2,1])
EL_S4_bbox.original <- EL_S4@images$Donor3.008um$centroids@bbox
EL_S4_ratio.original <- (EL_S4_bbox.original[1,2] - EL_S4_bbox.original[1,1]) / (EL_S4_bbox.original[2,2] - EL_S4_bbox.original[2,1])

#------Subsetting Objects for spots above 10 reads------#
EL_S1 <- subset(EL_S1, LowCounts == F)
EL_S1 #430,054 spots retained
EL_S2 <- subset(EL_S2, LowCounts == F)
EL_S2 #393,721 spots retained
EL_S3 <- subset(EL_S3, LowCounts == F)
EL_S3 #248,145 spots retained
EL_S4 <- subset(EL_S4, LowCounts == F)
EL_S4 #283,382 spots retained

#------Merging Datasets-----#
HHA.Merged <- merge(EL_S1, list(EL_S2, EL_S3, EL_S4))
HHA.Merged
rm(EL_S1, EL_S2, EL_S3, EL_S4)

#-----Performing sketch-based clustering and DR of 8 um data-----#
#Finding variable features and scaling data
DefaultAssay(HHA.Merged) <- "Spatial.008um"
HHA.Merged <- FindVariableFeatures(HHA.Merged, selection.method = "vst", nfeatures = 2000)
HHA.Merged <- ScaleData(HHA.Merged)
#Performing Sketch-based clustering
HHA.Merged <- SketchData(object = HHA.Merged, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
HHA.Merged

#-----Using Sketched Cells for DR/Clustering-----#
DefaultAssay(HHA.Merged) <- "sketch"
#Performing Dimensional Reduction on sketched data
HHA.Merged <- FindVariableFeatures(HHA.Merged)
HHA.Merged <- ScaleData(HHA.Merged)
HHA.Merged <- RunPCA(HHA.Merged, assay = "sketch", reduction.name = "pca.sketch")
HHA.Merged <- FindNeighbors(HHA.Merged, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
HHA.Merged <- RunUMAP(HHA.Merged, reduction = "pca.sketch", reduction.name = "umap.sketch.unintegrated", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "unintegrated.sketched_1.0", resolution = 1)
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "unintegrated.sketched_2.0", resolution = 2)
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "unintegrated.sketched_3.0", resolution = 3)

#------Plotting UMAP of Sketched Data-----#
DefaultAssay(HHA.Merged) <- "sketch"
Idents(HHA.Merged) <- "unintegrated.sketched_3.0"
UnintegratedSketchUMAP <- DimPlot(HHA.Merged, reduction = "umap.sketch.unintegrated", label = T, group.by = "unintegrated.sketched_3.0", raster = F) + 
  labs(x = "UMAP 1 (Sketched)", y = "UMAP 2 (Sketched)", title = "Sketched Clustering (50,000 Spots)")
#Saving as png
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_Merged_Unintegrated_UMAP_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
UnintegratedSketchUMAP
dev.off()
#Plotting by dataset
#Saving as png
UnintegratedSketchUMAPByDataset <- DimPlot(HHA.Merged, reduction = "umap.sketch.unintegrated", label = T, group.by = "orig.ident", raster = F) + 
  labs(x = "UMAP 1 (Sketched)", y = "UMAP 2 (Sketched)", title = "Sketched Clustering (100,000 Spots)")
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_Merged_Unintegrated_UMAP_By_Dataset_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
UnintegratedSketchUMAPByDataset #Minor batch effect
dev.off()

#-----Performing Sketch integration with harmony------#
set.seed(1701)
HHA.Merged <- IntegrateLayers(HHA.Merged, method = HarmonyIntegration, assay = "sketch", orig = "pca.sketch", new.reduction = "harmony.sketch")
#Reclustering and DR in Harmony Space
HHA.Merged <- FindNeighbors(HHA.Merged, assay = "sketch", reduction = "harmony.sketch", dims = 1:50)
HHA.Merged <- RunUMAP(HHA.Merged, reduction = "harmony.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
#Clustering cells using Louvain algorithm
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "harmony.sketched_1.0", resolution = 1)
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "harmony.sketched_2.0", resolution = 2)
HHA.Merged <- FindClusters(HHA.Merged, cluster.name = "harmony.sketched_3.0", resolution = 3)

#------Plotting UMAP-----#
SketchUMAP <- DimPlot(HHA.Merged, reduction = "umap.sketch", label = T, group.by = "harmony.sketched_3.0", raster = F) + 
  labs(x = "UMAP 1 (Sketched)", y = "UMAP 2 (Sketched)", title = "Sketched Clustering (200,000 Spots)")
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_Sketch_UMAP_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
SketchUMAP #Minor batch effect
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_Sketch_UMAP_by_Sample_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
DimPlot(HHA.Merged, reduction = "umap.sketch", label = F, group.by = "orig.ident", raster = F, shuffle = T) + 
  labs(x = "UMAP 1 (Sketched)", y = "UMAP 2 (Sketched)", title = "Sketched Clustering (200,000 Spots)")
dev.off()

#Saving 
save.image(paste0(WorkspaceDirectory, "EL_S1_S2_Integration_12_10_24.Rdata"))

#-----Projecting Integration to Full Dataset-----#
HHA.Merged <- ProjectIntegration(object = HHA.Merged, sketched.assay = "sketch", assay = "Spatial.008um", reduction = "harmony.sketch")

#-----Projecting full dataset onto sketched data-----#
HHA.Merged <- ProjectData(object = HHA.Merged, assay = "Spatial.008um", full.reduction = "harmony.sketch.full",
                         sketched.assay = "sketch", sketched.reduction = "harmony.sketch.full", umap.model = "umap.sketch", dims = 1:50,
                         refdata = list(harmony.projected_1.0 = "harmony.sketched_1.0", harmony.projected_2.0 = "harmony.sketched_2.0",
                                        harmony.projected_3.0 = "harmony.sketched_3.0"))

#Saving 
save.image(paste0(WorkspaceDirectory, "EL_S1_S2_Integration_12_10_24.Rdata"))



#------Plotting UMAP of Sketched Data-----#
FullUMAP <- DimPlot(HHA.Merged, reduction = "full.umap.sketch", label = T, raster = F, group.by = "harmony.projected_2.0") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)")

ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_Full_UMAP_12_10_24.png"), res = 600, height = 10, width = 12, units = "in")
FullUMAP
dev.off()

DefaultAssay(HHA.Merged) <- "Spatial.008um"
FeaturePlot(HHA.Merged, reduction = "full.umap.sketch", features = c("KRT5", "KRT1", "IVL", "FLG"), raster = F) & scale_color_gradientn(colors = GrayMagma)

#-----Getting Markers from Sketched Data-----#
SketchedMarkers.TFIDF <- HHA.Merged %>% JoinLayers() %>% GetAssayData(assay = "sketch", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Merged$harmony.sketched_2.0[!is.na(HHA.Merged$harmony.sketched_2.0)])
#Converting to List and Naming 
SketchedMarkersList.TFIDF <- SketchedMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(SketchedMarkersList.TFIDF) <- paste("Cluster", map(SketchedMarkersList.TFIDF, ~ .x$cluster[1]), "Markers")
SketchedMarkersList.TFIDF[22:41]

#-----Plotting Clusters-----#
Idents(HHA.Merged) <- "harmony.projected_2.0"
PlotClusters <- function(SeuratObj, Cluster, image) {
  cat("Plotting Cluster", Cluster, "\n")
  cells <- CellsByIdentities(SeuratObj, idents = c(Cluster))
  HighlightPlot <- SpatialDimPlot(SeuratObj, shape = 22,
                                  cells.highlight = cells[setdiff(names(cells), "NA")],
                                  cols.highlight = c("blue", "grey50"), facet.highlight = T, combine = T,
                                  alpha = c(0.25, 1), images = image) + NoLegend()
  return(HighlightPlot)
}
#Saving to pdf
pdf(paste0(PlotsDirectory, "EL_S2_Merged_Spatial_Clustering_Res2_12_10_24.pdf"))
EL_S2.Clusters <- HHA.Merged[[]] %>% filter(orig.ident == "EL_S2") %>% pull(harmony.projected_2.0)
for (i in 0:(length(levels(HHA.Merged))-1)) {
  if (i %in% EL_S2.Clusters == T) { #Avoids issues with low quality sample specific clusters throwing errors. 
  print(PlotClusters(HHA.Merged, i, image = "Donor1_2.008um"))}}
dev.off()
#Saving to pdf
pdf(paste0(PlotsDirectory, "EL_S1_Merged_Spatial_Clustering_Res2_12_10_24.pdf"))
EL_S1.Clusters <- HHA.Merged[[]] %>% filter(orig.ident == "EL_S1") %>% pull(harmony.projected_2.0)
for (i in 0:(length(levels(HHA.Merged))-1)) {
  if (i %in% EL_S1.Clusters == T) { #Avoids issues with low quality sample specific clusters throwing errors. 
    print(PlotClusters(HHA.Merged, i, image = "Donor1_1.008um"))}}
dev.off()

ClustersGrouped_3.0 <- list(
  c("Suprabasal ORS I", "Suprabasal ORS"), #0
  c("Dermis", "Dermis"), #1
  c("Dermis", "Dermis"), #2
  c("Dermal Sheath", "Dermal Sheath"), #3
  c("Granular Layer", "Upper HF/IFE"), #4
  c("Suprabasal ORS II", "Suprabasal ORS"), #5
  c("Suprabasal Infundibulum", "Upper HF/IFE"), #6
  c("Basal/Spinous IFE", "Upper HF/IFE"), #7
  c("Sebaceous IV", "Sebaceous"), #8
  c("Sebaceous Duct", "Sebaceous"), #9
  c("Endothelial", "Endothelial"), #10
  c("Upper Companion I", "Upper Companion"), #11
  c("Peri-Eccrine", "Eccrine"), #12
  c("Upper Companion II", "Upper Companion"), #13
  c("Sebaceous III", "Sebaceous"), #14
  c("Stratum Corneum", "Upper HF/IFE"), #15
  c("Dermis", "Dermis"), #16
  c("Dermis", "Dermis"), #17
  c("Proximal ORS", "Basal ORS"), #18
  c("Cortex II", "Cortex and Cuticle"), #19
  c("Eccrine I", "Eccrine"), #20
  c("Dermis", "Dermis"), #21
  c("Dermis", "Dermis"), #22
  c("Basal Infundibulum", "Upper HF/IFE"), #23
  c("Upper IRS", "IRS"), #24
  c("Bulge Suprabasal", "Suprabasal ORS"), #25
  c("Sebaceous II", "Sebaceous"), #26
  c("Upper IRS", "IRS"), #27
  c("Arrector Pili/VSM I", "Muscle"), #28
  c("Dermis", "Dermis"), #29
  c('Eccrine II', 'Eccrine'), #30
  c("Lower IRS", "IRS"), #31
  c("Dermis", "Dermis"), #32
  c("Matrix", "Matrix"), #33
  c("Outer Bulge", "Outer Bulge"), #34
  c("Cuticle", "Cuticle"), #35
  c("Basal ORS I", "Basal ORS"), #36
  c("Arrector Pili/VSM II", "Muscle"), #37
  c("Sebaceous I", "Sebaceous"), #38
  c("Cortex III", "Cortex"), #39
  c("Eccrine III", "Eccrine"), #40
  c("Sebaceous Sheath", "Sebaceous"), #41
  c("Bulge Dermal Sheath", "Dermal Sheath"), #42
  c("Catagen Club Hair", "Catagen"), #43
  c("Cortex I", "Cortex"), #44
  c("Basal ORS II", "Basal ORS"), #45
  c("Proximal Companion", "Companion"), #46
  c("Dermal Papilla", "Dermal Papilla"), #47
  c("Dermis", "Dermis"), #48
  c("Dermis", "Dermis"), #49
  c("Basal ORS III", "Basal ORS"), #50
  c("Dermis", "Dermis"), #51
  c("Dermis", "Dermis"), #52
  c("Dermis", "Dermis"), #53
  c("Dermis", "Dermis"), #54
  c("Dermis", "Dermis"), #55
  c("Poor Quality", "Poor Quality"), #56
  c("Poor Quality", "Poor Quality"), #57
  c("Dermis", "Dermis"), #58
  c("Dermis", "Dermis"), #59
  c("Upper Cuticle", "Cuticle"), #60
  c("Dermis", "Dermis"), #61
  c("Cortex IV", "Cortex"), #62
  c("Dermis", "Dermis"), #63
  c("IFE Melanocyte", "Upper HF/IFE"), #64
  c("Dermis", "Dermis"), #65
  c("Dermis", "Dermis"), #66
  c("Catagen II", "Catagen"), #67
  c("Dermis", "Dermis"), #68
  c("Dermis", "Dermis"), #69
  c("Dermis", "Dermis"), #70
  c("Dermis", "Dermis"), #71
  c("Dermis", "Dermis"), #72
  c("Dermis", "Dermis"), #73
  c("Dermis", "Dermis"), #77
  c("Dermis", "Dermis") #78
) %>% map(~ data.frame(InitialAnnotationFine = .x[1], InitialAnnotationBroad = .x[2])) %>% bind_rows()

#--------Renaming Idents--------#
InitialGlobalAnnotation <- ClustersGrouped_3.0$InitialAnnotationFine
InitialBroadAnnotation <- ClustersGrouped_3.0$InitialAnnotationBroad
#Broad Annotation
Idents(HHA.Merged) <- HHA.Merged$harmony.projected_2.0
names(InitialBroadAnnotation) <- sort(as.numeric(levels(HHA.Merged)))
HHA.Merged <- RenameIdents(HHA.Merged, InitialBroadAnnotation)
HHA.Merged$InitialAnnotationBroad <- HHA.Merged@active.ident
#Fine Annotation
Idents(HHA.Merged) <- HHA.Merged$harmony.projected_2.0
names(InitialGlobalAnnotation) <- sort(as.numeric(levels(HHA.Merged)))
HHA.Merged <- RenameIdents(HHA.Merged, InitialGlobalAnnotation)
HHA.Merged$InitialAnnotationFine <- HHA.Merged@active.ident

#------Plotting UMAP of Sketched Data-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_InitialFineAnnotation_UMAP_12_11_24.png"), res = 600, height = 10, width = 12, units = "in")
DimPlot(HHA.Merged, reduction = "full.umap.sketch", label = T, repel = T, raster = F, group.by = "InitialAnnotationFine") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)") + NoLegend()
dev.off()

#------Plotting UMAP of Sketched Data-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_InitialBroadAnnotation_UMAP_12_11_24.png"), res = 600, height = 10, width = 12, units = "in")
DimPlot(HHA.Merged, reduction = "full.umap.sketch", label = T, repel = T, raster = F, group.by = "InitialAnnotationBroad") + 
  labs(x = "UMAP 1 (Projected)", y = "UMAP 2 (Projected)", title = "Projected Clustering (1,355,302 Spots)")
dev.off()

#------Plotting Spatial DimPlot of Annotated Clusters-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_InitialBroadAnnotation_SpaialDimPlot_12_11_24.png"), res = 600, height = 10, width = 12, units = "in")
SpatialDimPlot(HHA.Merged, group.by = "InitialAnnotationBroad", images =  "Donor1_2.008um", label = T)
dev.off()
                          
SpatialDimPlot(HHA.Merged, group.by = "InitialAnnotationBroad", images =  "Donor1_2.008um")
SpatialDimPlot(HHA.Merged, group.by = "harmony.projected_2.0", images =  "Donor1_2.008um")


PlotClusters(HHA.Merged, 67, image = "Donor3.008um")

rm(EL_S1.Clusters)

#------Saving Seurat Object as Rds-----#
table(HHA.Merged$InitialAnnotationFine)
save(HHA.Merged, file = paste0(SeuratDirectory, "HHA_Integrated_Spatial_Initial_Annotations_12_11_24.rds"))

#---------Saving to AnnData------#
#----Writing Metadata----#
Full.metadata <- HHA.Merged[[]] %>% mutate(UMAP_1 = HHA.Merged@reductions$full.umap.sketch@cell.embeddings[,1], UMAP_2 = HHA.Merged@reductions$full.umap.sketch@cell.embeddings[,2])
Full.metadata$barcode <- rownames(Full.metadata)
write.csv(Full.metadata, file = paste0(AnnDataDirectory, "HHA_Spatial_Full_Conversion_Metadata_4_30_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
HHA.Merged <- JoinLayers(HHA.Merged)
#Pulling Normalized Counts Data
Full.logcounts <- HHA.Merged %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'Spatial.008um', layer = 'data')
Matrix::writeMM(Full.logcounts, file = paste0(AnnDataDirectory, "HHA_Spatial_Full_Conversion_Expression_LogCounts_4_30_25.mtx"))
#Pulling Log Transformed, Unnormalized Counts
Full.logcounts_unnorm <- HHA.Merged %>% SeuratObject::LayerData(assay = 'Log_008um', layer = 'data')
Matrix::writeMM(Full.logcounts_unnorm, file = paste0(AnnDataDirectory, "HHA_Spatial_Full_Conversion_Expression_LogCounts_Unnormalized_4_30_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Full.logcounts)), file = paste0(AnnDataDirectory, "HHA_Spatial_Full_Conversion_GeneNames_4_30_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing intermediate objects
rm(Full.metadata, Full.logcounts, Full.logcounts_unnorm)



