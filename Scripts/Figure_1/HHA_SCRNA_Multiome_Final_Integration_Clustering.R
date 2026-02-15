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

#-------Loading Checkpoint------#
load(paste0(WorkspaceDirectory, "HHA_SCRNA_Multiome_Integration_Workspace_3_31_25.Rdata"))

#-------Loading in RNA/Multiome Data------#
#SCRNA
load("/projects/p31989/HHA_Spatial_SC_Atlas/SeuratData/HHA_SCVI_Annotated_10_26_24.rds")
HHA
#Multiome
load(paste0(SeuratDirectory, "HHA_Multiome_Integrated_Annotated_3_2_25.rds"))
HHA.Multiome

#-----Merging objects-----#
HHA.Full <- merge(HHA, HHA.Multiome) #128,125 cells
HHA.Full
#Removing individual objects
rm(HHA.Multiome, HHA)

#-------Scoring Dissociation Signature using AUCell-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation 
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Calculating stress response score with AUCell
Full.StressAUC <- HHA.Full %>% JoinLayers() %>% RunAUCell(list(DissocResponse = StressSig$Gene_Human))
#Adding to Seurat Object
HHA.Full$DissociationScore <- Full.StressAUC$DissocResponse
#Visualizing by Dataset with violin plot
VlnPlot(HHA.Full, features = "DissociationScore", group.by = "SampleID", alpha = 0, raster = F)
#Visualizing by Platform with violin plot
VlnPlot(HHA.Full, features = "DissociationScore", group.by = "Platform", alpha = 0, raster = F)

#-----Unintegrated Clustering/DR-----#
DefaultAssay(HHA.Full) <- "RNA"
#Selecting 2000 variable features for dimensionality reduction 
HHA.Full <- FindVariableFeatures(HHA.Full, selection.method = "vst", nfeatures = 2500)
HHA.Full.TopFeaturePlot <- VariableFeaturePlot(HHA.Full) %>% LabelPoints(points = head(VariableFeatures(HHA.Full), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.Full.TopFeaturePlot
HHA.Full <- ScaleData(HHA.Full)
HHA.Full <- RunPCA(HHA.Full)
DimHeatmap(HHA.Full, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA.Full <-FindNeighbors(HHA.Full, dims = 1:50, reduction = "pca")
HHA.Full <- RunUMAP(HHA.Full, dims = 1:50, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA.Full <- FindClusters(HHA.Full, resolution = 1, cluster.name = "unintegrated_clusters")

#------------Plotting Unintegrated Clustering---------#
#Plotting unintegrated UMAP
DimPlot(HHA.Full, reduction = "UMAP_Unintegrated", raster = F) + 
  labs(title = "Unintegrated Datasets: Louvain Clusters", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15), plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Sample
DimPlot(HHA.Full, reduction = "UMAP_Unintegrated", group.by = "orig.ident", raster = F, cols = MetBrewer::met.brewer("Signac", 16), shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting By Platform 
HHA.Full$Platform <- factor(HHA.Full$Platform, levels = c("V3.1", "V4", "Multiome"))
DimPlot(HHA.Full, reduction = "UMAP_Unintegrated", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting QC Metrics
FeaturePlot(HHA.Full, reduction = "UMAP_Unintegrated", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),raster = F) & NoAxes() & 
  scale_color_viridis_c(option = "inferno")
#Plotting Dissociation Score
FeaturePlot(HHA.Full, reduction = "UMAP_Unintegrated", features = c("DissociationScore"),raster = F) & NoAxes() & scale_color_viridis_c(option = "inferno")

#-----Removing Cell Cycle Genes and Stress Genes from HVGs-----#
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
VariableFeatures(HHA.Full) <- VariableFeatures(HHA.Full)[!VariableFeatures(HHA.Full) %in% CellCycleGenes] 
VariableFeatures(HHA.Full) %>% length() #1984 genes
#Removing Stress Signature Genes
VariableFeatures(HHA.Full) <- VariableFeatures(HHA.Full)[!VariableFeatures(HHA.Full) %in% StressSig$Gene_Human] 
VariableFeatures(HHA.Full) %>% length() #1977 genes

#----Preparing for Integration-----#
#-------Performing Scanpy Conversion for SCVI Integration-----#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#Subsetting Atlas for Variable Genes 
Diet.Full <- HHA.Full[VariableFeatures(HHA.Full)]
#----Writing Metadata----#
Full.metadata <- Diet.Full[[]] 
Full.metadata$barcode <- rownames(Full.metadata)
write.csv(Full.metadata, file = paste0(AnnDataDirectory, "SCRNA_Multiome_Metadata_PreIntegration_Final_3_31_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting count matrix
DefaultAssay(Diet.Full) <- "RNA"
DietFull.counts <- Diet.Full %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
#Writing to mtx file 
Matrix::writeMM(DietFull.counts, file = paste0(AnnDataDirectory, "SCRNA_Multiome_RawCounts_Final_3_31_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietFull.counts)), file = paste0(AnnDataDirectory, "SCRNA_Multiome_GeneNames_Final_3_31_25.csv"), 
            quote = F, row.names = F, col.names = T)

#-------Getting Stress Signature Associated Genes-----#
#Subsetting Seurat Object for genes in stress signature 
Full.StressSig <- HHA.Full[StressSig$Gene_Human,]
#Writing as csv
StressSig.GeneMat <- Full.StressSig %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
StressSig.GeneMat$Barcode <- rownames(StressSig.GeneMat)
#Writing to csv
write_csv(StressSig.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_StressSig_LogNorm_Final_3_31_25.csv"))

#-------Getting Stress Signature Associated Genes-----#
#Subsetting Seurat Object for genes in stress signature 
Full.CellCycle <- HHA.Full[CellCycleGenes,]
CellCycle.GeneMat <- Full.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "SCRNA_Multiome_CellCycle_LogNorm_Final_3_31_25.csv"))
#Removing superfluous objects
rm(Diet.Full, Full.metadata, DietFull.counts, Full.CellCycle, Full.StressSig, CellCycle.GeneMat, StressSig.GeneMat)

#-----Getting SCVI Integration-----#
Full.scvi <- read_csv(paste0(AnnDataDirectory, "SCRNA_Multiome_RNA_LatentRep_Final_3_31_25.csv")) %>% column_to_rownames("Barcode")
head(Full.scvi)
# put it back in our original Seurat object
HHA.Full[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Full.scvi), key = "scvi_")
rm(Full.scvi)

#-----Performing Final DR and Clustering in scVI Space-----#
HHA.Full <-FindNeighbors(HHA.Full, dims = 1:50, reduction = "scvi")
HHA.Full <- RunUMAP(HHA.Full, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Full<- FindClusters(HHA.Full, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA.Full<- FindClusters(HHA.Full, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA.Full<- FindClusters(HHA.Full, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Full<- FindClusters(HHA.Full, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(HHA.Full, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.85", raster = F, pt.size = 0.01) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Sample ID-----#
DimPlot(HHA.Full, reduction = "UMAP", label = F, shuffle = T, group.by = "SampleID", raster = F,
        cols = MetBrewer::met.brewer("Signac", 16), pt.size = 0.1) + theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)

#-----Plotting by Isolation Method-----#
HHA.Full$Isolation <- factor(HHA.Full$Isolation, levels = c("4mm", "FUE", "Bulb"))
DimPlot(HHA.Full, reduction = "UMAP", label = T, group.by = "scvi_clusters_1.0", split.by = "Isolation", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()

#-----Plotting by Platform-----#
#Split
DimPlot(HHA.Full, reduction = "UMAP", label = T, group.by = "scvi_clusters_1.0", split.by = "Platform", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()
#Grouped
DimPlot(HHA.Full, reduction = "UMAP", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick", "darkseagreen")) + 
  labs(title = "Integrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoAxes()

#-----Plotting by Sex-----#
#Split
DimPlot(HHA.Full, reduction = "UMAP", label = T, group.by = "scvi_clusters_1.0", split.by = "Sex", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()
#Grouped
DimPlot(HHA.Full, reduction = "UMAP", group.by = "Sex", raster = F, shuffle = T, cols = c("steelblue3", "firebrick")) + 
  labs(title = "Integrated Datasets: Sex", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2", color = "Sex") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)

#-----Plotting by scRNAseq Initial Annotation -----#
DimPlot(HHA.Full, reduction = "UMAP", group.by = "InitialAnnotation", raster = F, shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2", color = "Initial Annotation") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting by Multiome initial Annotation 
DimPlot(HHA.Full, reduction = "UMAP", group.by = "MultiomeAnnotationRNA", raster = F, shuffle = T) + 
  labs(title = "Multiome Initial Annotation", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Plotting QC Metrics
FeaturePlot(HHA.Full, reduction = "UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
            raster = F) & NoAxes() & scale_color_viridis_c(option = "inferno")
#Plotting QC Metrics
FeaturePlot(HHA.Full, reduction = "UMAP", features = c("DissociationScore"),
            raster = F) & NoAxes() & scale_color_viridis_c(option = "plasma") & theme(aspect.ratio = 1)
#Plotting Proliferative Markers
FeaturePlot(HHA.Full, reduction = "UMAP", features = c("MKI67"),
            raster = F) & NoAxes() & scale_color_gradientn(colors = GrayMagma)

#------Plotting Marker Genes-------#
#Major Groups
FeaturePlot(HHA.Full, features = c("KRT5", "VIM", "COL1A1", "PECAM1", "PTPRC", "RGS5", "DCT", "MPZ", 'KRT7'), reduction = "UMAP", raster = T) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Matrix
FeaturePlot(HHA.Full, features = c("KRT35", "KRT85", "GATA3", "FABP4", "ALDH3A1", "FOXQ1", "KRT23", "KRT75", "TCHH"), reduction = "UMAP", raster = T) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#ORS
FeaturePlot(HHA.Full, features = c("KRT15", "KRT17", "KRT6A", "DIO2"), reduction = "UMAP", raster =T) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Immune Cells 
FeaturePlot(HHA.Full, features = c("PTPRC", "CD3G", "FOXP3", "NKG7", "CD68", "CD207", "MZB1", "CD79A", "MS4A1"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)
#Upper HF and IFE
FeaturePlot(HHA.Full, features = c("KRT5", "KRT1", "KRTDAP", "IVL", "S100A8", "IGFBP3", "KRT1", "SOX9", "MGST1"), reduction = "UMAP", raster = F) & 
  scale_color_gradientn(colors = GrayMako)
#Sebaceous Gland
FeaturePlot(HHA.Full, features = c("AWAT2", "SCD", "PPARG", "MGST1"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Proliferating
FeaturePlot(HHA.Full, features = c("MKI67"), reduction = "UMAP") & 
  scale_color_gradientn(colors = GrayMagma)

#-----Finding Markers for Each Cluster-----#
#Getting markers for each cluster 
FullMarkers.TFIDF <- HHA.Full %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Full$scvi_clusters_1.0, N = 20, FDR = 0.01)
FullMarkers.TFIDF
#converting to list
FullMarkers.TFIDFList <- FullMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(FullMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Full$scvi_clusters_1.0))-1))
FullMarkers.TFIDFList[1:21]
FullMarkers.TFIDFList[22:41]
FullMarkers.TFIDFList[42:47]

#-----Annotating Clusters------#
FullAnnotations <- list(
  c("Suprabasal_ORS_I", "Lower_HF"), #0
  c("Bulge", "Lower_HF"), #1
  c("Early_Cortex", "Matrix_Derived"), #2
  c("Suprabasal_IFE", "Upper_HF_IFE"), #3
  c("Basal_IFN", "Upper_HF_IFE"), #4
  c("Late_Cortex", "Matrix_Derived"), #5
  c("Basal_IFE", "Upper_HF_IFE"), #6
  c("Basal_ORS", "Lower_HF"), #7
  c("T_Helper", "T_Cells"), #8
  c("Suprabasal_IFN", "Upper_HF_IFE"), #9
  c("Fibroblasts_I", "Fibroblasts"), #10
  c("Fibroblasts_II", "Fibroblasts"), #11
  c("Sebaceous", "Upper_HF_IFE"), #12
  c("Suprabasal_Isthmus", "Lower_HF"), #13
  c("Fibroblasts_III", "Fibroblasts"), #14
  c("SM_DS_I", "SM_Dermal_Sheath"), #15
  c("Myeloid_I", "APCs"), #16
  c("Inner_Root_Sheath", "Matrix_Derived"), #17
  c("Proximal_ORS", "Matrix_Derived"), #18
  c("Endothelial", "Endothelial_LV"), #19
  c("Cuticle", "Matrix_Derived"), #20
  c("SM_DS_II", "SM_Dermal_Sheath"), #21
  c("Upper_Companion", "Lower_HF"), #22
  c("Melanocytes_I", "Melanocytes"), #23
  c("Suprabasal_ORS_II", "Lower_HF"), #24
  c("T_Reg", "T_Cells"), #25
  c("NK_Cells", "T_Cells"), #26
  c("Plasma_B_Cells", "B_Cells"), #27
  c("Melanocytes_II", "Melanocytes"), #28
  c("Early_Cortex", "Matrix_Derived"), #29
  c("Pericytes", "Endothelial_LV"), #30
  c("Langerhans_Cells", "APCs"), #31
  c("Myeloid_II", "APCs"), #32
  c('Proximal_Companion', 'Matrix_Derived'), #33
  c("Eccrine_Myoepithelial", "Eccrine"), #34
  c("Schwann_I", "Glial"), #35
  c("Mast_Cells", "Mast_Cells"), #36
  c("Eccrine_Epithelial_I", "Eccrine"), #37
  c("Eccrine_Epithelial_II", "Eccrine"), #38
  c("Myeloid_III", "APCs"), #39
  c("Memory_B_Cells", "B_Cells"), #40
  c("Lymph_Vasculature", "Endothelial_LV"), #41
  c("Fibroblasts_IV", "Fibroblasts"), #42
  c("Arrector_Pili", "Arrector_Pili"), #43
  c("Dermal_Papilla", "Fibroblasts"), #44
  c("Schwann_II", "Glial"), #45
  c("Merkel_Cells", "Merkel_Cells") #44
) %>% map(~ data.frame(FineAnnotation = .x[1], BroadAnnotation= .x[2])) %>% bind_rows()

#--------Renaming Idents--------#
GlobalAnnotation <- FullAnnotations$FineAnnotation
BroadAnnotation <- FullAnnotations$BroadAnnotation
#Broad Annotation
Idents(HHA.Full) <- "scvi_clusters_1.0"
names(BroadAnnotation) <- levels(HHA.Full)
HHA.Full <- RenameIdents(HHA.Full, BroadAnnotation)
HHA.Full$BroadAnnotation <- HHA.Full@active.ident
DimPlot(HHA.Full, reduction = "UMAP", group.by = "BroadAnnotation", label = T, repel = T, raster = F)  & NoLegend() & NoAxes() & theme(aspect.ratio = 1)
#Fine Annotation
Idents(HHA.Full) <- "scvi_clusters_1.0"
names(GlobalAnnotation) <- levels(HHA.Full)
HHA.Full <- RenameIdents(HHA.Full, GlobalAnnotation)
HHA.Full$GlobalAnnotation <- HHA.Full@active.ident
DimPlot(HHA.Full, reduction = "UMAP", group.by = "GlobalAnnotation", label = T, repel = T, raster = F)  & 
  theme(aspect.ratio = 1) & NoLegend() & NoAxes()

#-----Factoring Levels-----#
HHA.Full$GlobalAnnotationOrdered <- factor(HHA.Full$GlobalAnnotation, levels = c("Early_Cortex", "Late_Cortex", "Cuticle", "Inner_Root_Sheath", "Proximal_ORS", "Proximal_Companion", "Upper_Companion",
                                                                                 "Bulge", "Basal_ORS", "Suprabasal_ORS_I", "Suprabasal_ORS_II", "Suprabasal_Isthmus",
                                                                                 "Basal_IFN", "Suprabasal_IFN", "Sebaceous", "Basal_IFE", "Suprabasal_IFE",
                                                                                 "Merkel_Cells", "Eccrine_Epithelial_I", "Eccrine_Epithelial_II", "Eccrine_Myoepithelial", "Arrector_Pili", "SM_DS_I", "SM_DS_II",  
                                                                                 "Endothelial", "Pericytes", "Lymph_Vasculature", "Fibroblasts_I", "Fibroblasts_II", "Fibroblasts_III", "Fibroblasts_IV", "Dermal_Papilla",
                                                                                 "Schwann_I", "Schwann_II", "Melanocytes_I", "Melanocytes_II",
                                                                                 "Myeloid_I", "Myeloid_II", "Myeloid_III", "Langerhans_Cells", "Mast_Cells",
                                                                                 "T_Helper", "T_Reg", "NK_Cells", "Memory_B_Cells", "Plasma_B_Cells"))
#Plotting Factored Data
DimPlot(HHA.Full, reduction = "UMAP", label = T, repel = T, raster = F, group.by = "GlobalAnnotationOrdered")  & NoLegend() &
  NoAxes() & theme(aspect.ratio = 1)

#-----Plotting Dot Plot-----#
DotPlot(HHA.Full,  group.by = "GlobalAnnotationOrdered", features = c("MSX2", "KRT35", "KRT31", "KRT32", "KRT82", "TCHH", "KRT28", "GATA3", "KRT75",'KRT5', "SOX9", "DIO2", "TCEAL2","LGR5", "CTNND2", "BARX2", "KRT6A", "S100A8", "S100A9", "SAA1", "IGFBP3", "WNT3", "KRT1", 'KRT20', "KRT18", "KRT7",
                                                                      "CHRM3", "DES", "ACTA1", "ACTA2", "ITGA8",  "RGS5", "PECAM1", "CDH5", "PROX1", "CCL21", "PDGFRA", "DCN", "RSPO3", "CDH19", "NRXN1", "MPZ", "DCT", "PMEL",
                                                                      "PTPRC", "FCER1G", "CD14", 'CD68', "CD207", "TPSB2", "CD3D", "CD3G", "FOXP3", "NKG7", "MS4A1", "CD79A", "MZB1",
                                                                      "IGHG1", 'VIM')) + RotatedAxis() + scale_color_gradientn(colors = GrayMagma) +
  theme(aspect.ratio = 1, panel.background = element_rect(color = "black"))

#------Updating Seurat Object with Final Annotations-------#
#Clusters renamed during subsequent analysis but assignments did not change. 
Full.Metadata <- read_csv("/projects/b1217/HHA/Multiome_Seurat/HHA_FullAtlas_Final_Metadata.csv") %>% column_to_rownames("...1")
Full.Metadata
#Checking to make sure columns in the same order
sum(rownames(Full.Metadata) == colnames(HHA.Full))
HHA.Full$FinalAnnotation <- Full.Metadata$FinalAnnotation
HHA.Full$GeneralAnnotation <- Full.Metadata$GeneralAnnotation
AnnotationPal <- c('#91921a', '#aeaeae', '#ad8fc3', '#103d5d', '#185a88', '#4a2d27',
    '#86acdd', '#aec7e8', '#d6e2f3', '#76cbda', '#6b4239', '#8c564b',
    '#aa6d60', '#bc8b81', '#ffa145', '#ffd5ab', '#1f77b4', '#ab1f20',
    '#df5152', '#da4daf', '#9edae5', '#df62b9', '#57a9e2', '#217821',
    '#37c837', '#f288b6', '#e377c2', '#e78ccb', '#eca1d5', '#bebebe',
    '#1294a1', '#31d7e8', '#ddd1e7', '#c6e9f0', '#2b93db', '#ff6663',
    '#bcbd22', '#dadb37', '#fce4ee', '#d0d0d0', '#e0e0e0')
#Plotting Clusters
DimPlot(HHA.Full, reduction = "UMAP", group.by = "FinalAnnotation", label = T, repel = T, raster = F, cols = AnnotationPal) + 
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()

#----Writing Multiome Metadata to csv-----#
Multiome.Metadata <- HHA.Full[[]] %>% filter (Platform == "Multiome")
Multiome.Metadata
write_csv(Multiome.Metadata, file = paste0(SeuratDirectory, "HHA_Multiome_Metadata_Final.csv"))

#------Saving Object and Workspace-----#
save(HHA.Full, file = paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))
save.image(paste0(WorkspaceDirectory, "HHA_SCRNA_Multiome_Integration_Workspace_3_31_25.Rdata"))

#-----Recording session information-----#
sessionInfo()
#writing to text file
writeLines(capture.output(sessionInfo()), "/home/cms2775/HHA_Spatial_scRNAseq_2024/sessionInfo/HHA_SCRNA_Multiome_Final_Integration_Clustering_3_31_25_sessionInfo.txt")
