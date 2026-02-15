#-----Loading in Packages-----#
library(readr)
library(dplyr)
library(ggplot2)
library(ArchR)
library(purrr)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)


#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Multiome_SeuratPlots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"
#Setting Working Directory to Store ArchR Outputs
ArchRDirectory <- "/projects/b1217/HHA/Multiome_ArchR/"
setwd(ArchRDirectory)

#-----Setting up Parallel Computing-----#
addArchRThreads(threads = 32) 

#-----Setting up HG38 Genome-----#
addArchRGenome("hg38")

#------Palettes-----#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
#From MetBrewer::met.brewer("Signac", 10) "#92c051" "#2b9b81" "#1f6e9c" "#633372" "#9f5691" "#e87b89" "#de597c" "#d8443c" "#fe9b00" "#f4c40f"
SampleIDPal <- c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f")
Metadata.Pal <- scale_fill_manual(values = c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f"))

#------Loading ArchR Project------#
HHA.archr <- loadArchRProject("Multiome_ATAC_Level_Analysis_5_24_25")
HHA.archr

#-----Loading in Dataframe with Bulb Multiome Cells------#
Matrix.Barcodes <- read_csv(paste0(SeuratDirectory, 'HHA_Multiome_Matrix_Barcodes_5_30_25.csv'), col_names = T)
Matrix.Barcodes

#--------Subsetting ArchR Project for IFE cells------#
Matrix.Cells <- HHA.archr$cellNames[which(HHA.archr$UnifiedBarcode %in% Matrix.Barcodes$UnifiedBarcode)]
Matrix.Cells %>% length() #7,728 nuclei

#-----Subsetting Project for Bulb Keratinocytes-----#
#Subsetting for selected cells 
Matrix.archr <- subsetArchRProject(ArchRProj = HHA.archr, cells = Matrix.Cells, outputDirectory = "Multiome_Matrix_5_30_25", dropCells = TRUE,
                                 threads = getArchRThreads(), force = T)
Matrix.archr

#------Adding Identities to Matrix.archr object-------#
#Reordering metadata to match ArchR
Matrix.Metadata <- full_join(data.frame(UnifiedBarcode = Matrix.archr$UnifiedBarcode), Matrix.Barcodes, by = "UnifiedBarcode")
sum(Matrix.Metadata$UnifiedBarcode == Matrix.archr$UnifiedBarcode)
#Adding Identities
Matrix.archr$MatrixAnnotationBroad <- Matrix.Metadata$MatrixAnnotationBroad
Matrix.archr$MatrixAnnotationFine <- Matrix.Metadata$MatrixAnnotationFine
Matrix.archr$DonorID <- as.character(Matrix.archr$DonorID)

#-------Saving ArchR Project-----#
Matrix.archr <- saveArchRProject(ArchRProj = Matrix.archr, outputDirectory = "Multiome_Matrix", load = T)
#Loading in archR object
Matrix.archr <- loadArchRProject("Multiome_Matrix_5_30_25")
Matrix.archr

#------Plotting FRIP by Sample------#
#Plotting by Sample
plotGroups(ArchRProj = Matrix.archr, groupBy = "Sample", alpha = 0.6,
           colorBy = "cellColData", name = "FRIP", plotAs = "violins", pal = SampleIDPal)
#Plotting by Donor
Matrix.archr$DonorID <- as.character(Matrix.archr$DonorID)
plotGroups(ArchRProj = Matrix.archr, groupBy = "DonorID", alpha = 0.6,
           colorBy = "cellColData", name = "FRIP", plotAs = "violins")
#Plotting by Broad Grouping
plotGroups(ArchRProj = Matrix.archr, groupBy = "MatrixAnnotationBroad", alpha = 0.6,
           colorBy = "cellColData", name = "FRIP", plotAs = "violins")

#--------Computing iterative LSI----------#
Matrix.archr <- addIterativeLSI(ArchRProj = Matrix.archr,
                                useMatrix = "PeakMatrix", 
                                name = "IterativeLSI", 
                                iterations = 4, 
                                clusterParams = list(resolution = c(0.1, 0.2, 0.4), n.start = 10), 
                                varFeatures = 25000, 
                                dimsToUse = 1:30, 
                                force = T)
#Adding clusters
Matrix.archr <- addClusters(input = Matrix.archr, reducedDims = "IterativeLSI",
                            method = "Seurat", name = "UnintegratedClusters", resolution = 0.5, force = T)
#Running UMAP
Matrix.archr <- addUMAP(ArchRProj = Matrix.archr, reducedDims = "IterativeLSI", name = "UMAP_Unintegrated", 
                        nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)

#------Plotting Clustering Results-----#
#Plotting Clusters vs Samples: Clusters are relatively well mixed by Sample
ggAlignPlots(plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "UnintegratedClusters", embedding = "UMAP_Unintegrated", size = 1),
             plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Unintegrated", size = 0.5), type = "h")
#Plotting Annotations: Clusters are relatively well mixed by Sample
ggAlignPlots(plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationFine", embedding = "UMAP_Unintegrated", size = 1),
             plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationBroad", embedding = "UMAP_Unintegrated", size = 1), type = "h")

#------Performing Integration with Harmony-----#
#integration
Matrix.archr <- addHarmony(ArchRProj = Matrix.archr, reducedDims = "IterativeLSI",
                           name = "Harmony", groupBy = c("Sample", "DonorID"), force = T)
#Adding clusters
Matrix.archr <- addClusters(input = Matrix.archr, reducedDims = "Harmony",
                            method = "Seurat", name = "HarmonyClusters_0.2", resolution = 0.2, force = T)
Matrix.archr <- addClusters(input = Matrix.archr, reducedDims = "Harmony",
                            method = "Seurat", name = "HarmonyClusters_0.3", resolution = 0.3, force = T)
Matrix.archr <- addClusters(input = Matrix.archr, reducedDims = "Harmony",
                            method = "Seurat", name = "HarmonyClusters_0.4", resolution = 0.4, force = T)
#Running UMAP
Matrix.archr <- addUMAP(ArchRProj = Matrix.archr, reducedDims = "Harmony", name = "UMAP_Harmony", 
                        nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)

#Plotting Embedding
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "HarmonyClusters_0.4", embedding = "UMAP_Harmony", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)
#Plotting Embedding
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)
#Plotting Embedding
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationBroad", embedding = "UMAP_Harmony", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)
#Plotting Embedding
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationFine", embedding = "UMAP_Harmony", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)

#-------Saving ArchR Project-----#
Matrix.archr <- saveArchRProject(ArchRProj = Matrix.archr, outputDirectory = "Multiome_Matrix_5_30_25", load = T)
#Loading in archR object
Matrix.archr <- loadArchRProject("Multiome_Matrix_5_30_25")
Matrix.archr

#------Pulling UMAP from Object and Sending to Seurat-----#
ATAC.UMAP <- getEmbedding(ArchRProj = Matrix.archr, embedding = "UMAP_Harmony")
ATAC.UMAP$Barcode <- str_replace(rownames(ATAC.UMAP), pattern = "#", replacement = "_")
colnames(ATAC.UMAP) <- c("UMAP_1", "UMAP_2", "UnifiedBarcode")
head(ATAC.UMAP)
write_csv(ATAC.UMAP, file = paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_UMAP_5_30.csv"))

#------Pulling Harmony Embedding from Object and Sending to Seurat-----#
ATAC.Harmony <- getReducedDims(ArchRProj = Matrix.archr, reducedDims = "Harmony") %>% as.data.frame()
ATAC.Harmony$UnifiedBarcode <- str_replace(rownames(ATAC.Harmony), pattern = "#", replacement = "_")
head(ATAC.Harmony)
write_csv(ATAC.Harmony, file = paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_Harmony_5_30.csv"))

#-------Pulling in WNN UMAP from Seurat------#
#Reading in Embedding
WNN.UMAP <- read_csv(file = paste0(SeuratDirectory, "HHA_Matrix_Multiome_WNN_UMAP_5_30.csv"))
WNN.UMAP
#Reordering to match ArchR
sum(WNN.UMAP$UnifiedBarcode == Matrix.archr$UnifiedBarcode)
WNN.UMAP <- full_join(data.frame(UnifiedBarcode = Matrix.archr$UnifiedBarcode), WNN.UMAP, by = "UnifiedBarcode")
sum(WNN.UMAP$UnifiedBarcode == Matrix.archr$UnifiedBarcode)
rownames(WNN.UMAP) <-Matrix.archr$cellNames
WNN.UMAP$UnifiedBarcode <- NULL
#Matching ArchR format
colnames(WNN.UMAP) <- c("WNN#UMAP1", "WNN#UMAP2")
#Adding to Archr Object
Matrix.archr@embeddings$WNNUMAP <- SimpleList(df = WNN.UMAP, params = list())
#Plotting Cells
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationBroad", embedding = "WNNUMAP", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)
#Plotting Cells
plotEmbedding(ArchRProj = Matrix.archr, colorBy = "cellColData", name = "MatrixAnnotationFine", embedding = "WNNUMAP", size = 1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_ArchR(legendTextSize = 10)

#-------Saving ArchR Project-----#
Matrix.archr <- saveArchRProject(ArchRProj = Matrix.archr, outputDirectory = "Multiome_Matrix_5_30_25", load = T)
#Loading in archR object
Matrix.archr <- loadArchRProject("Multiome_Matrix_5_30_25")
Matrix.archr

#-------Getting Peaks from ArchR----------#
#Pulling peak information from Seurat Object
Matrix.Peaks <- getPeakSet(Matrix.archr)
Matrix.Peaks
#Avoiding redundancy in colnames
names(Matrix.Peaks) <- paste0("Peak", 1:length(Matrix.Peaks))
#Pulling Peak Matrix from ArchR Object
Matrix.PeakCounts <- getMatrixFromProject(Matrix.archr, useMatrix='PeakMatrix')
Matrix.PeakCounts #297,443 peaks
#Extracting Matrix
Matrix.PeakMatrix <- as(assay(Matrix.PeakCounts, "PeakMatrix"), "dgCMatrix")
#Double Checking Peak order is the same
check1 <- as.data.frame(Matrix.Peaks)[,1:3]
rownames(check1) <- paste0("Peak", 1:length(Matrix.Peaks))
check2 <- as.data.frame(Matrix.PeakCounts@rowRanges)[,1:3]
rownames(check2) <- paste0("Peak", 1:length(Matrix.Peaks))
identical(check1, check2)

#-----Saving Data for Anndata Conversion-----#
#----Writing Metadata----#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#Adding metadata
Matrix.Metadata <- getCellColData(Matrix.archr) %>% as.data.frame() %>% cbind(data.frame(ATAC_UMAP1 = ATAC.UMAP$UMAP_1,
                                                                                         ATAC_UMAP2 = ATAC.UMAP$UMAP_2,
                                                                                         WNN_UMAP1 = WNN.UMAP$`WNN#UMAP1`,
                                                                                         WNN_UMAP2 = WNN.UMAP$`WNN#UMAP2`))
Matrix.Metadata$barcode <- rownames(Matrix.Metadata)
#Writing to df
sum(Matrix.Metadata$barcode == colnames(Matrix.PeakMatrix)) #Has to be reordered to match peak matrix
Matrix.Metadata <- full_join(data.frame(barcode = colnames(Matrix.PeakMatrix)), Matrix.Metadata, by = "barcode")
sum(Matrix.Metadata$barcode == colnames(Matrix.PeakMatrix)) #Now correctly ordered
write.csv(Matrix.Metadata, file = paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_5_30_25.csv"), quote = F, row.names = F)

#----Writing Peak Matrix and Region Names----#
#Writing peak matrix to mtx file 
Matrix::writeMM(Matrix.PeakMatrix, file = paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_PeakMatrix_5_30_25.mtx"))
#writing region names 
write_csv(as.data.frame(Matrix.Peaks), file = paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_PeakMatrix_RegionNames_5_30_25.csv"))

#-------Writing Harmony Embeddings-----#
ATAC.Harmony$barcode <- rownames(ATAC.Harmony)
sum(rownames(ATAC.Harmony) == colnames(Matrix.PeakMatrix)) #Has to be reordered to match peak matrix
ATACHarmony.Reordered <- full_join(data.frame(barcode = colnames(Matrix.PeakMatrix)), ATAC.Harmony, by = "barcode") 
sum(ATACHarmony.Reordered$barcode == colnames(Matrix.PeakMatrix)) #Has to be reordered to match peak matrix
write_csv(ATACHarmony.Reordered[,1:30], file = paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Harmony_5_30_25.csv"))

#------Writing Gene Scores-----#
#These won't really be used, only including to avoid potential issues with SEACells workflow if not included in the object
# Gene scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(Matrix.archr)
var_features <- Matrix.archr@reducedDims[["IterativeLSI"]]$LSIFeatures
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
Matrix.archr <- addGeneScoreMatrix(Matrix.archr, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)
#Exporting Results
gene.scores <- getMatrixFromProject(Matrix.archr)
scores <- as(assay(gene.scores, 'GeneScoreMatrix'), "dgCMatrix")
rownames(scores) <- rowData(gene.scores)$name
Matrix::writeMM(scores, file = paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_GeneScores_5_30_25.mtx"), quote=FALSE)

#-------Computing Marker Peaks------#
Matrix.markerPeaks <- getMarkerFeatures(ArchRProj = Matrix.archr,
                                      useMatrix = "PeakMatrix", 
                                      groupBy = "MatrixAnnotationBroad",
                                      bias = c("TSSEnrichment", "log10(nFrags)"),
                                      testMethod = "wilcoxon")
#pulling markers into list
Matrix.markerList <- getMarkers(Matrix.markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
#Plotting Heatmap
Matrix.heatmapPeaks <- markerHeatmap(
  seMarker = Matrix.markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE)
draw(Matrix.heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#------Determining if marker peak differences can be found between the upper and lower compartments-----#
COL17.markerTest <- getMarkerFeatures(
  ArchRProj = Matrix.archr, 
  useMatrix = "PeakMatrix",
  groupBy = "MatrixAnnotationFine",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Upper COL17",
  bgdGroups = "Lower COL17")
#Plotting MA Plot of Markers
COL17.MA <- plotMarkers(seMarker = COL17.markerTest, name = "Upper COL17", cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5", plotAs = "MA")
COL17.MA
#Plotting MA Plot of Markers
COL17.Volcano <- plotMarkers(seMarker = COL17.markerTest, name = "Upper COL17", cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5", plotAs = "volcano")
COL17.Volcano

#--------Retrieving DA Peaks--------#
COL17.Markers <- getMarkers(COL17.markerTest, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")[[1]]
rownames(COL17.Markers) <- paste0(COL17.Markers$seqnames, ":", COL17.Markers$start, "-", COL17.Markers$end)
COL17.Markers$peak <- rownames(COL17.Markers)
COL17.Markers

#-------Writing DA Peaks to csv-------#
write_csv(as.data.frame(COL17.Markers), file = paste0(SeuratDirectory, "COL17_Upper_vs_Lower_DA_Peaks.csv"))
