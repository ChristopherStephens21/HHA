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
AnnDataDirectory <- "/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/"

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
stallion = c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC",
             "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8",  "#6E4B9E", "#0C727C", "#7E1416", "#D8A767")

#--------Loading Checkpoint-------#
load(paste0(WorkspaceDirectory, "HHA_Matrix_Multiome_Seurat_5_30_25.RData"))

#-------Loading Matrix Object-------#
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))

#--------Adding Column with All Clusters--------#
#Plotting UMAP
DimPlot(HHA.Matrix, reduction = "UMAP", group.by = "MatrixAnnotationFine", label = T, repel = T, raster = F, pt.size = 0.25, cols = stallion) 

#------Subsetting for Multiome Cells-----#
Matrix.Multiome <- HHA.Matrix %>% JoinLayers() %>% subset(Platform == "Multiome")
Matrix.Multiome #7,728 nuclei
DimPlot(Matrix.Multiome, reduction = "UMAP", group.by = "MatrixAnnotationFine", label = T, repel = T, raster = F, pt.size = 0.25, cols = stallion) 

#-----Reclustering after cleaning -----#
#Resplitting Layers
Matrix.Multiome[["RNA"]] <- split(Matrix.Multiome[["RNA"]], f = Matrix.Multiome$orig.ident)
Layers(Matrix.Multiome)
#Selecting 2000 variable features for dimensionality reduction
Matrix.Multiome <- FindVariableFeatures(Matrix.Multiome, selection.method = "vst", nfeatures = 2000)
Matrix.Multiome.TopFeaturePlot <- VariableFeaturePlot(Matrix.Multiome) %>% LabelPoints(points = head(VariableFeatures(Matrix.Multiome), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
Matrix.Multiome.TopFeaturePlot
#Scaling Data
Matrix.Multiome <- ScaleData(Matrix.Multiome)
#Running PCA
Matrix.Multiome <- RunPCA(Matrix.Multiome)
DimHeatmap(Matrix.Multiome, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
Matrix.Multiome <-FindNeighbors(Matrix.Multiome, dims = 1:30, reduction = "pca")
Matrix.Multiome <- RunUMAP(Matrix.Multiome, dims = 1:30, reduction = "pca", reduction.name = "UMAP_Unintegrated")
Matrix.Multiome <- FindClusters(Matrix.Multiome, resolution = 1, cluster.name = "unintegrated_clusters")

#------Getting Stress Signature Associated Genes-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Subsetting Seurat Object for genes in stress signature
Matrix.StressSig <- Matrix.Multiome[StressSig$Gene_Human,]
#------Pulling Cell Cycle Associated Genes-----#
#Tirosh et al Cell Cycle Genes
CanonicalCellCycle <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Subsetting Seurat Object for genes in stress signature
Matrix.CellCycle <- Matrix.Multiome[CanonicalCellCycle,]

#------Removing Dissociation Signature and Cell Cycle genes from Variable Features------#
#-----Removing Cell Cycle Genes and Stress Genes from HVGs-----#
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
VariableFeatures(Matrix.Multiome) <- VariableFeatures(Matrix.Multiome)[!VariableFeatures(Matrix.Multiome) %in% CellCycleGenes]
VariableFeatures(Matrix.Multiome) %>% length() #1986 genes
#Removing Stress Signature Genes
VariableFeatures(Matrix.Multiome) <- VariableFeatures(Matrix.Multiome)[!VariableFeatures(Matrix.Multiome) %in% StressSig$Gene_Human]
VariableFeatures(Matrix.Multiome) %>% length() #1970 genes genes

#-------Performing Scanpy Conversion for SCVI Integration: Cell Cycle and Stress Removed-----#
#Subsetting Atlas for Variable Genes
Diet.Matrix <- Matrix.Multiome[VariableFeatures(Matrix.Multiome)]
#----Writing Metadata----#
Matrix.metadata <- Diet.Matrix[[]]
Matrix.metadata$barcode <- rownames(Matrix.metadata)
write.csv(Matrix.metadata, file = paste0(AnnDataDirectory, "Multiome_Matrix_Metadata_5_30_25.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting log normalized counts
DefaultAssay(Diet.Matrix) <- "RNA"
#Writing log normalized counts to mtx
DietMatrix.counts <- Diet.Matrix %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(DietMatrix.counts, file = paste0(AnnDataDirectory, "Multiome_Matrix_RawCounts_5_30_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietMatrix.counts)), file = paste0(AnnDataDirectory, "Multiome_Matrix_GeneNames_5_30_25.csv"),
            quote = F, row.names = F, col.names = T)
#Removing superfluous objects
rm(Diet.Matrix, DietMatrix.counts, Matrix.metadata)
#------Writing Normalized Cell Cycle Gene Expression------#
#Writing as csv
CellCycle.GeneMat <- Matrix.CellCycle %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
CellCycle.GeneMat$Barcode <- rownames(CellCycle.GeneMat)
#Writing to csv
write_csv(CellCycle.GeneMat, file = paste0(AnnDataDirectory, "Multiome_Matrix_CellCycle_LogNorm_5_30_25.csv"))
#------Writing Normalized Stress Gene Expression------#
StressSig.GeneMat <- Matrix.StressSig %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data') %>% t() %>% as.data.frame()
StressSig.GeneMat$Barcode <- rownames(StressSig.GeneMat)
#Writing to csv
write_csv(StressSig.GeneMat, file = paste0(AnnDataDirectory, "Multiome_Matrix_StressSig_LogNorm_5_30_25.csv"))

#-----Getting SCVI Integration Results-----#
Matrix.scvi <- read_csv(paste0(AnnDataDirectory, "Multiome_RNA_Matrix_LatentRep_5_30_25.csv")) %>% column_to_rownames("Barcode")
head(Matrix.scvi)
# put it back in our original Seurat object
Matrix.Multiome[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(Matrix.scvi), key = "scvi_")
rm(Matrix.scvi)
#-----Performing Final DR and Clustering in scVI Space-----#
Matrix.Multiome <-FindNeighbors(Matrix.Multiome, dims = 1:50, reduction = "scvi")
Matrix.Multiome <- RunUMAP(Matrix.Multiome, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 1, cluster.name = "scvi_clusters_1.0")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.85, cluster.name = "scvi_clusters_0.85")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.35, cluster.name = "scvi_clusters_0.35")
Matrix.Multiome<- FindClusters(Matrix.Multiome, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#------Plotting UMAPs-----#
#-----Plotting Clusters-----#
DimPlot(Matrix.Multiome, reduction = "UMAP", group.by = "MatrixAnnotationFine", label = T, repel = T, raster = F, pt.size = 1, cols = stallion) 

#---------Saving Barcodes for ATAC analysis-------#
#Getting Barcodes
Matrix.ArchR <- Matrix.Multiome[[]] %>% dplyr::select(UnifiedBarcode, MatrixAnnotationBroad, MatrixAnnotationFine)
write_csv(Matrix.ArchR, paste0(SeuratDirectory, 'HHA_Multiome_Matrix_Barcodes_5_30_25.csv'), col_names = T)

#------Adding in ATAC UMAP to Seurat------#
#Loading in UMAP from ATAC
ATAC.UMAP <- read_csv(paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_UMAP_5_30.csv"))
sum(ATAC.UMAP$UnifiedBarcode == Matrix.Multiome$UnifiedBarcode)
ATAC.UMAP <- full_join(data.frame(UnifiedBarcode = Matrix.Multiome$UnifiedBarcode), ATAC.UMAP, by = "UnifiedBarcode")
sum(ATAC.UMAP$UnifiedBarcode == Matrix.Multiome$UnifiedBarcode)
ATAC.UMAP <- ATAC.UMAP %>% column_to_rownames("UnifiedBarcode")
rownames(ATAC.UMAP) <- rownames(Matrix.Multiome[[]])
#Adding DR
Matrix.Multiome[["UMAP_ATAC"]] <- CreateDimReducObject(embeddings = as.matrix(ATAC.UMAP), key = "UMAP_ATAC_")
#Visualizing
DimPlot(Matrix.Multiome, reduction = "UMAP_ATAC", label = T, group.by = "MatrixAnnotationBroad", raster = F, pt.size = 0.5)

#-----Adding Harmony Embedding to Seurat-----#
#Loading in Harmony from ATAC
ATAC.Harmony <- read_csv(paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_Harmony_5_30.csv"))
sum(ATAC.Harmony$UnifiedBarcode == Matrix.Multiome$UnifiedBarcode)
ATAC.Harmony <- full_join(data.frame(UnifiedBarcode = Matrix.Multiome$UnifiedBarcode), ATAC.Harmony, by = "UnifiedBarcode")
sum(ATAC.Harmony$UnifiedBarcode == Matrix.Multiome$UnifiedBarcode)
ATAC.Harmony <- ATAC.Harmony %>% column_to_rownames("UnifiedBarcode")
rownames(ATAC.Harmony) <- rownames(Matrix.Multiome[[]])
#Adding DR
Matrix.Multiome[["Harmony_ATAC"]] <- CreateDimReducObject(embeddings = as.matrix(ATAC.Harmony), key = "Harmony_ATAC_")

#--------Running WNN UMAP Computation------#
#Using Seurat WNN
Matrix.Multiome <- FindMultiModalNeighbors(Matrix.Multiome, reduction.list = list("scvi", "Harmony_ATAC"), dims.list = list(1:50, 1:30))
#Running UMAP on Multimodal Neighbors
Matrix.Multiome <- RunUMAP(Matrix.Multiome, nn.name = "weighted.nn", reduction.name = "UMAP_WNN", reduction.key = "WNNUMAP_")
#Visualizing
DimPlot(Matrix.Multiome, reduction = "UMAP_WNN", label = T, group.by = "MatrixAnnotationBroad", raster = F, pt.size = 0.5)
#Visualizing
DimPlot(Matrix.Multiome, reduction = "UMAP_WNN", label = T, group.by = "MatrixAnnotationFine", raster = F, cols = stallion)
#Writing to csv for archr
WNN.UMAP <- data.frame(UMAP_1 =Matrix.Multiome@reductions$UMAP_WNN@cell.embeddings[,1],
                       UMAP_2 = Matrix.Multiome@reductions$UMAP_WNN@cell.embeddings[,2],
                       UnifiedBarcode = Matrix.Multiome$UnifiedBarcode)
#writing to csv
write_csv(WNN.UMAP, file = paste0(SeuratDirectory, "HHA_Matrix_Multiome_WNN_UMAP_5_30.csv"))

#--------Saving Image-------#
save.image(paste0(WorkspaceDirectory, "HHA_Matrix_Multiome_Seurat_5_30_25.RData"))

#--------Exporting to Python for Downstream Analysis------#
#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
Matrix.metadata <- Matrix.Multiome[[]] %>% mutate(UMAP_1 = Matrix.Multiome@reductions$UMAP@cell.embeddings[,1], 
                                                  UMAP_2 = Matrix.Multiome@reductions$UMAP@cell.embeddings[,2],
                                                  WNN_UMAP_1 = Matrix.Multiome@reductions$UMAP_WNN@cell.embeddings[,1],
                                                  WNN_UMAP_2 = Matrix.Multiome@reductions$UMAP_WNN@cell.embeddings[,2])
Matrix.metadata$barcode <- rownames(Matrix.metadata)
Matrix.metadata <- Matrix.metadata
write.csv(Matrix.metadata, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_FullConversion_Metadata_5_30_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(Matrix.Multiome@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_FullConversion_scVI_Embeddings_5_30_25.csv"), quote = F, row.names = F)

#------Exporting SNN Graph Connectivities------#
#Getting long format DF containing all possible connections and their presence/absence
Matrix.SNN <- as.data.frame(summary(Matrix.Multiome@graphs$RNA_nn)) %>% #Matching Pythonic Indexing (starting at 0)
  mutate(i = i-1, j = j-1)
write_csv(Matrix.SNN, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_SNN_Connectivities_5_30_25.csv"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Matrix.logcounts <- Matrix.Multiome %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(Matrix.logcounts, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_FullConversion_Expression_LogCounts_5_30_25.mtx"))
#Pulling Raw Counts Data
Matrix.rawcounts <- Matrix.Multiome %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(Matrix.rawcounts, file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_FullConversion_RawCounts_5_30_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Matrix.logcounts)), file = paste0(AnnDataDirectory, "HHA_Multiome_Matrix_FullConversion_GeneNames_5_30_25.csv"),
            quote = F, row.names = F, col.names = T)

#-------Saving Multiome Matrix Object------#
save(Matrix.Multiome, file = paste0(SeuratDirectory, file = "HHA_Multiome_Matrix_Integrated_Annotated_5_30_25.rds"))

