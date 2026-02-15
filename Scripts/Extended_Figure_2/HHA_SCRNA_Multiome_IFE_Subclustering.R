#------Subclustering of Upper HF and IFE-----#
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
PlotsDirectory <- "/projects/b1217/HHA/Multiome_SeuratPlots/"
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

#------Loading Checkpoint-----#
load(paste0(WorkspaceDirectory, "HHA_UpperHF_IFE_Subclustering_4_1_25.RData"))

#-----Loading Global Object-----#
load(paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))

#--------------Analysis of Upper HF and IFE Populations------------------#
#Checking for Bulb cells in Upper HF
table(HHA.Upper$BroadAnnotation, HHA.Upper$Isolation) # 1 Bulb cell present in Upper HF/IFE (from sample EL_B8)
#Selecting cells, removing EL_B8 Bulb singlet
HHA.Upper <- subset(HHA.Full, BroadAnnotation %in% c("Upper_HF_IFE") & Isolation != "Bulb")
HHA.Upper #31,208 cells
#Sanity checking selection
DimPlot(HHA.Upper, reduction = "UMAP", label = T, group.by = "GlobalAnnotation", raster = F, pt.size = 0.01)
#rm(HHA.Full)
#Looking at cell counts by sample and platform
table(HHA.Upper$SampleID)
table(HHA.Upper$Platform)

#-----Finding Variable Features-----#
#Selecting 2000 variable features for dimensionality reduction 
HHA.Upper <- FindVariableFeatures(HHA.Upper, selection.method = "vst", nfeatures = 2000)
HHA.Upper.TopFeaturePlot <- VariableFeaturePlot(HHA.Upper) %>% LabelPoints(points = head(VariableFeatures(HHA.Upper), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Upper.TopFeaturePlot
#Scaling Data
HHA.Upper <- ScaleData(HHA.Upper)
#Running PCA
HHA.Upper <- RunPCA(HHA.Upper)
DimHeatmap(HHA.Upper, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Upper <-FindNeighbors(HHA.Upper, dims = 1:30, reduction = "pca")
HHA.Upper <- RunUMAP(HHA.Upper, dims = 1:30, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Upper <- FindClusters(HHA.Upper, resolution = 1, cluster.name = "unintegrated_clusters")

#------Plotting UMAP-----#
#Plotting unintegrated UMAP
DimPlot(HHA.Upper, reduction = "UMAP_Unintegrated", raster = F) + 
  labs(title = "Unintegrated Datasets: Louvain Clusters", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15), plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Sample
DimPlot(HHA.Upper, reduction = "UMAP_Unintegrated", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Donor
DimPlot(HHA.Upper, reduction = "UMAP_Unintegrated", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting By Platform 
HHA.Upper$Platform <- factor(HHA.Upper$Platform, levels = c("V3.1", "V4", "Multiome"))
DimPlot(HHA.Upper, reduction = "UMAP_Unintegrated", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting By Annotation
DimPlot(HHA.Upper, reduction = "UMAP_Unintegrated", group.by = "GlobalAnnotation", raster = F, shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting QC Metrics
FeaturePlot(HHA.Upper, reduction = "UMAP_Unintegrated", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()
#Plotting Dissociation Score
FeaturePlot(HHA.Upper, reduction = "UMAP_Unintegrated", features = c("DissociationScore")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()

#------Getting Stress Signature Associated Genes-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation 
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Subsetting Seurat Object for genes in stress signature 
Upper.StressSig <- HHA.Upper[StressSig$Gene_Human,]

#------Pulling Cell Cycle Associated Genes-----#
#Tirosh et al Cell Cycle Genes
CanonicalCellCycle <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Subsetting Seurat Object for genes in stress signature 
Upper.CellCycle <- HHA.Upper[CanonicalCellCycle,]

#------Removing Dissociation Signature and Cell Cycle genes from Variable Features------#
#-----Removing Cell Cycle Genes and Stress Genes from HVGs-----#
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
VariableFeatures(HHA.Upper) <- VariableFeatures(HHA.Upper)[!VariableFeatures(HHA.Upper) %in% CellCycleGenes] 
VariableFeatures(HHA.Upper) %>% length() #1953 genes
#Removing Stress Signature Genes
VariableFeatures(HHA.Upper) <- VariableFeatures(HHA.Upper)[!VariableFeatures(HHA.Upper) %in% StressSig$Gene_Human] 
VariableFeatures(HHA.Upper) %>% length() #1932 genes genes


#-------Performing Scanpy Conversion for SCVI Integration: Cell Cycle Removed-----#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#Subsetting Atlas for Variable Genes 
Diet.Upper <- HHA.Upper[VariableFeatures(HHA.Upper)]
#----Writing Metadata----#
Upper.metadata <- Diet.Upper[[]]
Upper.metadata$barcode <- rownames(Upper.metadata)
write.csv(Upper.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Upper_Metadata_4_1_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.Upper) <- "RNA"
#Writing log normalized counts to mtx
DietUpper.counts <- Diet.Upper %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietUpper.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Upper_RawCounts_4_1_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietUpper.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_Upper_GeneNames_4_1_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.Upper, DietUpper.counts, Upper.metadata)

#------Writing Normalized Cell Cycle Gene Expression------#
#Writing as csv
CellCycle.GeneMat <- Upper.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Upper_CellCycle_LogNorm_4_1_25.csv"))
#------Writing Normalized Stress Gene Expression------#
StressSig.GeneMat <- Upper.StressSig %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
StressSig.GeneMat$Barcode <- rownames(StressSig.GeneMat)
#Weriting to csv
write_csv(StressSig.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Upper_StressSig_LogNorm_4_1_25.csv"))

#-----Getting SCVI Integration Results-----#
Upper.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_Upper_LatentRep_4_1_25.csv")) %>% column_to_rownames("Barcode")
head(Upper.scvi)
# put it back in our original Seurat object
HHA.Upper[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Upper.scvi), key = "scvi_")
rm(Upper.scvi)

#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Upper <-FindNeighbors(HHA.Upper, dims = 1:50, reduction = "scvi")
HHA.Upper <- RunUMAP(HHA.Upper, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.Upper<- FindClusters(HHA.Upper, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Upper, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.5", raster = F, pt.size = 0.5) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Sample-----#
DimPlot(HHA.Upper, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 11), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Donor-----#
DimPlot(HHA.Upper, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Donor") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#------Plotting By Platform-----#
DimPlot(HHA.Upper, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting By Annotation-----#
DimPlot(HHA.Upper, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "GlobalAnnotation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting QC Metrics
FeaturePlot(HHA.Upper, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()
#Plotting QC Metrics
FeaturePlot(HHA.Upper, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "inferno") & NoAxes()

#-----Plotting Marker Genes----#
#IFE Differentiation
FeaturePlot(HHA.Upper, features = c("KRT5", "KRT1", "KRTDAP", "IVL"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Upper HF versus IFE
FeaturePlot(HHA.Upper, features = c("WNT3", "KRT15", "GATA3", "S100A8", "S100A9", "SOX9", "MGST1", "DEFB1", "PHLDB1"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

FeaturePlot(HHA.Upper, features = c("SOX9"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-----Finding Markers for Each Cluster-----#
UpperMarkers.TFIDF <- HHA.Upper %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Upper$scvi_clusters_0.5, N = 20, FDR = 0.01)
UpperMarkers.TFIDF
UpperMarkers.TFIDFList <- UpperMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(UpperMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Upper$scvi_clusters_0.5))-1))
UpperMarkers.TFIDFList

#------Annotating Cells------#
UpperAnnotations <- c("Suprabasal_IFE_I", #0
                      "Basal_IFE", #1
                      "Basal_IFN", #2
                      "Suprabasal_IFN", #3
                      "Early_Seb", #4
                      "Suprabasal_IFE_II", #5
                      "Upper_Bulge", #6
                      "Sebocyte") #7)

#-----Renaming Idents-----#
Idents(HHA.Upper) <- "scvi_clusters_0.5"
names(UpperAnnotations) <- levels(HHA.Upper)
HHA.Upper <- RenameIdents(HHA.Upper, UpperAnnotations)
HHA.Upper$UpperAnnotation <- HHA.Upper@active.ident
DimPlot(HHA.Upper, reduction = "UMAP", group.by = "UpperAnnotation", label = T, repel = T, raster = F, pt.size = 0.5) &
  theme(aspect.ratio = 1)

DimPlot(HHA.Upper, reduction = "UMAP", group.by = "GlobalAnnotation", label = T, repel = T, raster = F, pt.size = 0.5) &
  theme(aspect.ratio = 1)
VlnPlot(HHA.Upper, features = "FOS", group.by = "UpperAnnotation")

#--------------Subclustering of IFE Populations---------#
#-------Subsetting for IFE Populations-------#
#Subsetting for IFE
HHA.IFE <- subset(HHA.Upper, UpperAnnotation %in% c("Basal_IFE", "Suprabasal_IFE_I", "Suprabasal_IFE_II"))
HHA.IFE #15,884 cells
#Sanity checking selection
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "UpperAnnotation", label = T, repel = T, raster = F, pt.size = 0.5) &
  theme(aspect.ratio = 1)

#------Getting Normalized Expression for Stress Signature Associated Genes-----#
IFE.StressSig <- HHA.IFE %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "data")
#Subsetting Seurat Object for genes in stress signature 
IFE.StressSig <- IFE.StressSig[rownames(IFE.StressSig) %in% StressSig$Gene_Human,] %>% t() %>% as.data.frame()
IFE.StressSig$Barcode <- rownames(IFE.StressSig)

#------Getting Normalized Expression for Cell Cycle Genes-----#
IFE.CellCycle <- HHA.IFE %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "data")
#Subsetting Seurat Object for genes in stress signature 
IFE.CellCycle <- IFE.CellCycle[rownames(IFE.CellCycle) %in% c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes),] %>% t() %>% as.data.frame()
IFE.CellCycle$Barcode <- rownames(IFE.CellCycle)

#-----Finding Variable Features-----#
#Selecting 2000 variable features for dimensionality reduction 
HHA.IFE <- FindVariableFeatures(HHA.IFE, selection.method = "vst", nfeatures = 2000)
#Removing Stress Signature genes from Variable Features
VariableFeatures(HHA.IFE) <- VariableFeatures(HHA.IFE)[VariableFeatures(HHA.IFE) %in% StressSig$Gene_Human ==F]
VariableFeatures(HHA.IFE) %>% length() #1975 genes
#Removing Cell Cycle Genes from Variable Features
VariableFeatures(HHA.IFE) <- VariableFeatures(HHA.IFE)[VariableFeatures(HHA.IFE) %in% CanonicalCellCycle ==F]
VariableFeatures(HHA.IFE) %>% length() #1923 genes

#-------Performing Scanpy Conversion for SCVI Integration: Cell Cycle Removed-----#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#Subsetting Atlas for Variable Genes 
Diet.IFE <- HHA.IFE[VariableFeatures(HHA.IFE)]
#----Writing Metadata----#
IFE.metadata <- Diet.IFE[[]]
IFE.metadata$barcode <- rownames(IFE.metadata)
write.csv(IFE.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_IFE_Metadata_4_1_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.IFE) <- "RNA"
#Writing log normalized counts to mtx
DietIFE.counts <- Diet.IFE %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietIFE.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_IFE_RawCounts_4_1_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietIFE.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_IFE_GeneNames_4_1_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.IFE, DietIFE.lognorm, IFE.metadata)
#------Writing Normalized Cell Cycle Gene Expression------#
#Writing to csv
write_csv(IFE.CellCycle, file = paste0(AnnDataDirectory, "SCRNA_Multiome_IFE_CellCycle_LogNorm_4_1_25.csv"))
#------Writing Normalized Stress Gene Expression------#
write_csv(IFE.StressSig, file = paste0(AnnDataDirectory, "SCRNA_Multiome_IFE_StressSig_LogNorm_4_1_25.csv"))

#-----Getting SCVI Integration Results-----#
IFE.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_IFE_LatentRep_4_1_25.csv")) %>% column_to_rownames("Barcode")
head(IFE.scvi)
# put it back in our original Seurat object
HHA.IFE[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(IFE.scvi), key = "scvi_")
rm(IFE.scvi)

#-----Performing Final DR and Clustering in scVI Space-----#
HHA.IFE <-FindNeighbors(HHA.IFE, dims = 1:50, reduction = "scvi")
HHA.IFE <- RunUMAP(HHA.IFE, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.IFE, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", raster = F, pt.size = 0.5) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Sample-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 11), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Donor-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Donor") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#------Plotting By Platform-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting By Annotation-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "GlobalAnnotation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#------Plotting QC Metrics------#
#Standard QC
FeaturePlot(HHA.IFE, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()
#Dissociation Signature
FeaturePlot(HHA.IFE, reduction = "UMAP", features = c("DissociationScore", "FOS"), pt.size = 0.5) &
  scale_color_viridis_c(option = "plasma") & NoAxes()
#Problem Genes
FeaturePlot(HHA.IFE, reduction = "UMAP", features = c("FOS", "SOX9", "S100A8", "WNT3"), pt.size = 0.5) &
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#------Removing Stress Gene Outlier Clusters------#
VlnPlot(HHA.IFE, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "DissociationScore", 'FOS', "SOX9"), group.by = "scvi_clusters_0.75",
        alpha = 0.1)
HHA.IFE <- subset(HHA.IFE, scvi_clusters_0.75 %in% c(7,9) == F)
HHA.IFE

#-----Saving list of retained cells to csv-----#
IFE.FilteredBarcodes <- data.frame(Barcode = colnames(HHA.IFE))
write_csv(IFE.FilteredBarcodes, file = paste0(AnnDataDirectory, "HHA_IFE_Stress_Filtered_Barcodes_4_1_25.csv"))

#-----Getting SCVI Integration Results Post-Filtering-----#
IFE.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_RNA_IFE_PostFiltering_LatentRep_4_1_25.csv")) %>% column_to_rownames("Barcode")
head(IFE.scvi)
# put it back in our original Seurat object
HHA.IFE[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(IFE.scvi), key = "scvi_")
rm(IFE.scvi)

#-----Performing Final DR and Clustering in scVI Space-----#
HHA.IFE <-FindNeighbors(HHA.IFE, dims = 1:50, reduction = "scvi")
HHA.IFE <- RunUMAP(HHA.IFE, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
HHA.IFE<- FindClusters(HHA.IFE, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.IFE, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.4", raster = F, pt.size = 0.5) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Sample-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "SampleID", raster = F, cols = MetBrewer::met.brewer("Signac", 11), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Donor-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "DonorID", raster = F, cols = MetBrewer::met.brewer("Signac", 9), shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Donor") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#------Plotting By Platform-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting By Annotation-----#
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "GlobalAnnotation", raster = F, shuffle = T) + 
  labs(x = "UMAP 1", y = "UMAP 2", color = "GlobalAnnotation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#------Plotting QC Metrics------#
#Standard QC
FeaturePlot(HHA.IFE, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) &
  scale_color_viridis_c(option = "inferno") & NoAxes()
#Dissociation Signature
FeaturePlot(HHA.IFE, reduction = "UMAP", features = c("DissociationScore"), pt.size = 0.5) &
  scale_color_viridis_c(option = "plasma") & NoAxes()

#-----Plotting Markers----#
#Differentiation
FeaturePlot(HHA.IFE, features = c("COL17A1", "KRT5", "MKI67", "TOP2A", "KRT10", "KRT1", "KRTDAP", "IVL","FLG"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Transcription Factors
FeaturePlot(HHA.IFE, features = c("GATA3", "SOX6", "POU3F1", "GRHL3"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)
#Stress Response
FeaturePlot(HHA.IFE, features = c("FOS", "JUN",  "GADD45B"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#-----Finding Markers for Each Cluster-----#
IFEMarkers.TFIDF <- HHA.IFE %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.IFE$scvi_clusters_0.4, N = 20, FDR = 0.01)
IFEMarkers.TFIDF-
  IFEMarkers.TFIDFList <- IFEMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(IFEMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.IFE$scvi_clusters_0.4))-1))
IFEMarkers.TFIDFList

FeaturePlot(HHA.IFE, features = c("SOX9"), reduction = "UMAP") &
  scale_color_gradientn(colors = GrayMagma) & NoAxes() & theme(aspect.ratio = 1)

#------Annotating Cells------#
IFEAnnotations <- c("Basal", #0
                    "Early_Granular", #1
                    "Spinous", #2
                    "Late_Granular") #4)
#-----Renaming Idents-----#
Idents(HHA.IFE) <- "scvi_clusters_0.4"
names(IFEAnnotations) <- levels(HHA.IFE)
HHA.IFE <- RenameIdents(HHA.IFE, IFEAnnotations)
HHA.IFE$IFEAnnotation <- HHA.IFE@active.ident
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "IFEAnnotation", label = F, raster = F, pt.size = 0.85) &
  theme(aspect.ratio = 1) & NoAxes()

#-------Plotting IFE Cluster DimPlot-----#
HHA.IFE$IFEAnnotationOrdered <- factor(HHA.IFE$IFEAnnotation, levels = c("Basal", "Spinous", "Early_Granular", "Late_Granular"))
#RColorBrewer::brewer.pal(9, "Blues")
IFEPal <- c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B")
DimPlot(HHA.IFE, reduction = "UMAP", group.by = "IFEAnnotationOrdered", pt.size = 0.6, raster = F, cols = IFEPal) &
  theme(aspect.ratio = 1) & NoAxes()


#------Saving Objects and Workspace------#
#HHA.Upper
save(HHA.Upper, file = paste0(SeuratDirectory, 'HHA_UpperHF_IFE_Subcluster_4_1_25.rds'))
#HHA.IFE
save(HHA.IFE, file = paste0(SeuratDirectory, 'HHA_IFE_Subcluster_4_1_25.rds'))
#Workspace
save.image(paste0(WorkspaceDirectory, "HHA_UpperHF_IFE_Subclustering_4_1_25.RData"))

#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
IFE.metadata <- HHA.IFE[[]] %>% mutate(UMAP_1 = HHA.IFE@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.IFE@reductions$UMAP@cell.embeddings[,2])
IFE.metadata$barcode <- rownames(IFE.metadata)
IFE.metadata <- IFE.metadata
write.csv(IFE.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_Metadata_4_1_25.csv"), quote = F, row.names = F)

#----Writing scVI Embeddings----#
write.csv(HHA.IFE@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_scVI_Embeddings_4_1_25.csv"), quote = F, row.names = F)

#----Writing Adjacency Matrices (Graphs)---#
#Nearest Neighbor Graph
RNAnn.adjmat <- as(object = HHA.IFE@graphs$RNA_nn, Class = "Matrix") #Converting to Matrix
RNAnn.adjmat <- as(object = RNAnn.adjmat, Class = "dgCMatrix") #Converting to Sparse
str(RNAnn.adjmat)
Matrix::writeMM(RNAnn.adjmat, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_NNAdjMat_4_1_25.mtx"))
#Shared Nearest Neighbor Graph
RNAsnn.adjmat <- as(object = HHA.IFE@graphs$RNA_snn, Class = "Matrix") #Converting to HHA.IFE
RNAsnn.adjmat <- as(object = RNAsnn.adjmat, Class = "dgCMatrix") #Converting to Matrix
str(RNAsnn.adjmat)
Matrix::writeMM(RNAsnn.adjmat, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_SNNAdjMat_4_1_25.mtx"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
IFE.logcounts <- HHA.IFE %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(IFE.logcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_Expression_LogCounts_4_1_25.mtx"))
#Pulling Raw Counts Data
IFE.rawcounts <- HHA.IFE %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(IFE.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_RawCounts_4_1_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(IFE.logcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_IFE_Conversion_GeneNames_4_1_25.csv"),
            quote = F, row.names = F, col.names = T)