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

#------Loading Checkpoint------#
load(paste0(SeuratDirectory, file = "HHA_Matrix_DORC_Scoring_1_7_26.RData"))

#---------Loading Objects------#
#ATAC Only
load(paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_Signac_6_16_25.rds"))
Matrix.ATAC
#Multiome RNA/ATAC
load(paste0(SeuratDirectory, "HHA_Multiome_ATAC_Matrix_Integrated_Annotated_6_16_25.rds"))
Matrix.Multiome
#Full RNA
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))

#------Ordering Clusters-------#
MatrixPal <- c("#89C75F", "#3BBCA8",  "#208A42", "#0C727C", "#9ECAE1", "#4292C6", "#08306B", "#E6C2DC", "#C06CAB",  "#89288F", "#D8A767", "#F47D2B", "#F37B7D",  "#7E1416", "#D24B27")
Matrix.Multiome$MatrixAnnotationOrdered <- factor(Matrix.Multiome$MatrixAnnotationFine, levels = c("Lower COL17", "Upper COL17", "LPC", "Medulla", 
                                                                                                   "Early_Cortex", "Middle_Cortex", "Late_Cortex", 
                                                                                                   "Early_Cuticle", "Middle_Cuticle", "Late_Cuticle",
                                                                                                   "Early_IRS_I", "Early_IRS_II", "IRS_Henle", "IRS_Huxley", "IRS_Cuticle"))
#Plotting
DimPlot(Matrix.Multiome, reduction = "UMAP", group.by = "MatrixAnnotationOrdered", raster = F, cols = MatrixPal) &
  NoAxes() & theme(aspect.ratio = 1) & labs(title = "")
#Looking at Multiome Clusters
DimPlot(Matrix.Multiome, reduction = "UMAP_WNN", label = T, group.by = "MatrixAnnotationOrdered", raster = F, cols = MatrixPal, pt.size = 0.5)

#------Reading in P2G Link Counts per Gene-----#
#Peak2Gene Linkages per Gene from 100 kb run 
P2G100.counts <- read_csv("/projects/b1217/HHA/Matrix_ATAC_SEACells/SEACells_Peak2GeneCounts_6_5_25.csv")
#Peak2Gene Linkages per Gene from 500 kb run 
P2G500.counts <- read_csv("/projects/b1217/HHA/Matrix_ATAC_SEACells/SEACells_HRG_Peak2GeneCounts_6_5_25.csv")

#------Plotting Gene Ranks from P2G Linkage Analysis and Identifying HRGs--------#
P2G100.rankedcounts <- P2G100.counts %>% arrange(n_peaks) %>% mutate(rank = 1:nrow(P2G100.counts), HRG = case_when(n_peaks >= 8 ~ "Highly Regulated",
                                                                                                                   .default = "Not Highly Regulated"))
P2G500.rankedcounts <- P2G500.counts %>% arrange(n_peaks) %>% mutate(rank = 1:nrow(P2G500.counts), HRG = case_when(n_peaks >= 8 ~ "Highly Regulated",
                                                                                                                   .default = "Not Highly Regulated"))
#---------Plotting Significant Peaks per Gene: 100 kb Search Space-------#
#Plotting Top HRGs
Top10HRG100 <- P2G100.rankedcounts %>% tail(10)
#Plotting Peaks per Gene
P2G100.RankPlot <- ggplot(P2G100.rankedcounts, aes(x = rank, y = n_peaks, color = HRG)) + geom_point() +
  geom_hline(yintercept = 7.5, linetype = "dashed") + scale_color_manual(values = c("cornflowerblue", "firebrick")) +
  labs(x = "Gene Rank", y = "Number of Linked Peaks", title = "Peak Linkages per Gene: 100 kB Search Space") + theme_bw() + 
  ggrepel::geom_text_repel(data = Top10HRG100, aes(x = rank, y = n_peaks, label = gene), color = "black", force = 5) +
  theme(aspect.ratio = 1)
P2G100.RankPlot
#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HRG_100kb_RankPlot.tiff"), res = 600, height = 6, width = 8, units = "in") 
P2G100.RankPlot
dev.off()

#---------Plotting Significant Peaks per Gene: 500 kb Search Space-------#
#500 kb search space for initially identified HRGs
#Plotting Top HRGs
Top10HRG500 <- P2G500.rankedcounts %>% tail(10)
#Plotting Peaks per Gene
P2G500.RankPlot <- ggplot(P2G500.rankedcounts, aes(x = rank, y = n_peaks, color = HRG)) + geom_point() +
  geom_hline(yintercept = 7.5, linetype = "dashed") + scale_color_manual(values = c("cornflowerblue", "firebrick")) +
  labs(x = "Gene Rank", y = "Number of Linked Peaks", title = "Peak Linkages per Gene: 500 kB Search Space") + theme_bw() + 
  ggrepel::geom_text_repel(data = Top10HRG500, aes(x = rank, y = n_peaks, label = gene), color = "black", force = 5) +
  theme(aspect.ratio = 1)
P2G500.RankPlot
#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HRG_500kb_RankPlot.tiff"), res = 600, height = 6, width = 8, units = "in") 
P2G500.RankPlot
dev.off()

#------Running GO Term Analysis on HRGs-------#
HRG.Genes <- P2G500.rankedcounts %>% filter(HRG == "Highly Regulated") %>% pull(gene)
HRG.Genes #1493 genes pass threshold

#-------Getting Enriched GO Terms for HRGs-----#
#Using all genes tested for peaks as the background gene set to minimize bias
HRG.GO <- enrichGO(HRG.Genes, OrgDb = "org.Hs.eg.db", keyType = 'SYMBOL', readable = T, universe = P2G100.rankedcounts$gene,
                   ont = "MF", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH")
#Plotting Top GO Terms
dotplot(HRG.GO, x = "p.adjust", color = "p.adjust", showCategory = 10)  & 
  scale_fill_viridis_c(option = "magma") & labs(x = "Adjusted P Value", y = "GOMF Term") & ggtitle("Highly Regulated Genes")
HRG.GOSim <- enrichplot::pairwise_termsim(HRG.GO)
#Plotting as Tree
enrichplot::treeplot(HRG.GOSim)

#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HRG_GO_Treeplot.tiff"), res = 600, height = 6, width = 8, units = "in") 
enrichplot::treeplot(HRG.GOSim)
dev.off()

#----Removing Stress Response and Cell Cycle Genes from DORCs-----#
#Loading in stress response gene list (van den Brink et al (2017))
StressSig <- read_csv("/projects/b1217/HHA/Dissociation_Signature/van_den_Brink_2017_Dissociation_DEGS.csv")
#Converting from mouse to human notation
StressSig$Gene_Human <- str_to_upper(StressSig$Gene)
#Tirosh et al Cell Cycle Genes
CanonicalCellCycle <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Filtering out DORCs from stress response genes 
HRG.StressCC <- HRG.Genes[HRG.Genes %in% c(StressSig$Gene_Human, CanonicalCellCycle)]
HRG.StressCC #2 cell cycle and 6 stress respoinse genes
#Removing from gene list
HRG.Genes <- HRG.Genes[HRG.Genes %in% HRG.StressCC == F]
HRG.Genes #1485 genes 

#--------Reading in Full Peak-Gene Information for HRGs-----#
P2G500.SigLinks <- read_csv("/projects/b1217/HHA/Matrix_ATAC_SEACells/Matrix_HRG_PG_Links_6_5_25.csv")
#Filtering for genes retained after cell cycle/stress removal 
P2G500.SigLinks <- P2G500.SigLinks %>% filter(gene %in% HRG.Genes)
#Writing filtered SigLinks to csv
write_csv(P2G500.SigLinks, "/projects/b1217/HHA/Matrix_ATAC_SEACells/Matrix_HRG_PG_Links_Stress_CC_Filtered_6_5_25.csv")

#-----Computing Significant Links to HRGs by Peak-----#
DORC.GeneCounts <- P2G500.SigLinks %>% group_by(peak) %>% dplyr::summarize(n_genes = n()) %>% arrange(desc(n_genes))
DORC.GeneCounts

#-------Getting Peaks in GRanges Format------#
#Pulling DORC Peaks
DORC.Peaks <- unique(P2G500.SigLinks$peak) #25,955 peaks are linked to HRGs
DORCPeaks.GR <- StringToGRanges(DORC.Peaks, sep = c(":", "-"))
#Adding Peak Names in Signac Format
DORCPeaks.GR$name <- str_replace(DORC.Peaks, pattern = ":", replacement = "-")
names(DORCPeaks.GR) <- DORCPeaks.GR$name
DORCPeaks.GR
#Saving as rds
save(DORCPeaks.GR, file = "/projects/b1217/HHA/Matrix_ATAC_SEACells/Matrix_DORC_Peak_GR.rds")

#----Reading in Peak Matrix and Region Names----#
#Reading in ATAC Cell Metadata
Matrix.ATACMetadata <- read_csv(paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_5_30_25.csv"))
Matrix.ATACMetadata$UnifiedBarcode
#Getting Barcodes in the correct order
RNA.Barcodes <- data.frame(UnifiedBarcode = Matrix.Multiome$UnifiedBarcode,
                           Barcode = colnames(Matrix.Multiome))
RNA.Barcodes$UnifiedBarcode %in% Matrix.ATACMetadata$UnifiedBarcode
#Getting barcodes matching RNA barcodes in the original ATAC order matching the peak matrix 
ATAC.Barcodes <- full_join(data.frame(UnifiedBarcode = Matrix.ATACMetadata$UnifiedBarcode), RNA.Barcodes, by = "UnifiedBarcode")
sum(ATAC.Barcodes$UnifiedBarcode == Matrix.ATACMetadata$UnifiedBarcode)
#Reading in Gene Names
Matrix.RegionNames <- read_csv(paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_PeakMatrix_RegionNames_5_30_25.csv"))
Matrix.RegionNames <- Matrix.RegionNames %>% mutate(Peak_Name = paste0(seqnames, ":", start, "-", end))
#Reading in Peak Matrix
Matrix.PeakMatrix <- Matrix::readMM(paste0(AnnDataDirectory, "HHA_Matrix_Multiome_ATAC_Metadata_PeakMatrix_5_30_25.mtx"))
rownames(Matrix.PeakMatrix) <- Matrix.RegionNames$Peak_Name
colnames(Matrix.PeakMatrix) <- ATAC.Barcodes$Barcode
str(Matrix.PeakMatrix)
#Subsetting Peak Matrix for DORC Peaks
DORC.PeakMatrix <- Matrix.PeakMatrix[DORC.Peaks,]

#------Normalizing Peaks by Total Fragments in Peaks-----#
#Normalizing to 1,000,000 fragments in peaks per cell
NormalizeByFIP <- function(CellExp) {return(CellExp/FIP * 1000000)}
Matrix.FragsInPeaks <- colSums(Matrix.PeakMatrix)
#Removing Peak Matrix
rm(Matrix.PeakMatrix)
#-------Normalizing Peaks-------#
FIP <- Matrix.FragsInPeaks
DORC.NormMat <- t(apply(DORC.PeakMatrix, MARGIN = 1, FUN=NormalizeByFIP))

#-------Retrieving Peaks for each DORC-------#
#Function that returns a list of the peaks associated with each DORC 
GetPeaksbyDORC <- function(HRGs, SigLinks) {
  DORCPeaks <- list()
  for (Gene in HRGs) {
    Peaks <- SigLinks %>% filter(gene == Gene) %>% pull(peak)
    DORCPeaks <- c(DORCPeaks, list(Peaks))}
  names(DORCPeaks) <- HRGs
  return(DORCPeaks)}
#Grouping Peaks by DORCs
PeaksByDORC <- GetPeaksbyDORC(HRG.Genes, P2G500.SigLinks)
PeaksByDORC

#-------Computing DORC Scores-------#
#Function for Computing DORC Scores
GetDORCScores <- function(PeakMat, DORCPeaks) {
  DORCScores.List <- list()
  i <- 1
  cat("Calculating DORC scores for DORC", i, "\n")
  for (PeakSet in DORCPeaks) {
    i <- i+ 1
    if (i %% 200 == 0) {cat("Calculating DORC scores for DORC", i, "\n")}
    #Filtering peak matrix for peaks in DORC
    DORCPeakMat <- PeakMat[PeakSet,]
    #Summing normalized counts across all peaks in DORC in each cell 
    DORCScores <- colSums(DORCPeakMat)
    DORCScores.List <- c(DORCScores.List, list(DORCScores))}
  names(DORCScores.List) <- names(DORCPeaks)
  cat("Creating DORC Matrix \n")
  #This generates the matrix, but messes up some of the gene names by replacing dashes with periods
  DORCScoreMat <- as.data.frame(DORCScores.List) %>% as.matrix() %>% t()
  rownames( DORCScoreMat) <- names(DORCScores.List) #Fixing gene names
  return(DORCScoreMat)}
#Computing DORC scores
DORCScore.Mat <- GetDORCScores(DORC.NormMat, PeaksByDORC)
#Adding Seurat Assay
Matrix.Multiome[["DORC_Scores"]] <- CreateAssayObject(data = DORCScore.Mat, key = "DORC")

#-----Calculating DADs for the COL17A1 Layers------#
#Logarithmizing DORC Matrix
Matrix.Multiome[["DORC_LogNorm"]] <- CreateAssayObject(data = log1p(DORCScore.Mat), key = "DORC_LogNorm")
#Subsetting for COL17A1+ Layer
DefaultAssay(Matrix.Multiome) <- "RNA"
HHA.COL17 <- Matrix.Multiome %>% JoinLayers() %>% subset(MatrixAnnotationFine %in% c("Lower COL17", "Upper COL17"))
table(HHA.COL17$MatrixAnnotationFine)
#Calculating Marker DORCs with wilcoxon
Idents(HHA.COL17) <- "MatrixAnnotationFine"
DefaultAssay(HHA.COL17) <- "DORC_LogNorm"
#Calculating DEGs between lower and upper COL17A1
DORC.LowervsUpper <- FindMarkers(HHA.COL17, ident.1 = "Upper COL17", ident.2 = "Lower COL17", test.use = "wilcox", assay = 'DORC_LogNorm', only.pos = F, verbose = T, logfc.threshold = 0, min.pct = 0)  %>% 
  mutate(diffexp = factor(case_when(p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Upregulated", 
                                    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Downregulated",
                                    .default = "Not Significant"), levels = c("Upregulated", "Downregulated", "Not Significant")))
DORC.LowervsUpper$gene <- rownames(DORC.LowervsUpper)

#-----Plotting Volcano Plot-----#
Top20Up <- DORC.LowervsUpper %>% filter(diffexp == "Upregulated") %>% arrange(p_val_adj) %>% head(20) %>% pull(gene)
Top20Down <- DORC.LowervsUpper %>% filter(diffexp == "Downregulated") %>% arrange(p_val_adj) %>% head(20) %>% pull(gene)
TopMarkers <- DORC.LowervsUpper %>% filter(gene %in% c(Top20Up, Top20Down))
SelectedMarkers <- 
  #Plotting
  DORC.VolcanoPlot <- ggplot(DORC.LowervsUpper, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp)) + geom_point() + geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')+
  geom_vline(xintercept = c(-0.25, 0.25), col = "black", linetype = 'dashed') +
  ggrepel::geom_text_repel(data = TopMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), color = "black", force = 1, max.overlaps = Inf) +
  theme_classic() + scale_color_manual(values = c("cornflowerblue", "firebrick", "gray"), labels = c("Enriched in Upper COL17", "Enriched in Lower COL17", "Not Significant")) + 
  labs(x = "Average Log2Fold Change", y = "- Log10 P Value", color = "Differential Regulation", title = "Differential DORC Accessibility: Lower vs Upper COL17A1 Layer")
DORC.VolcanoPlot

#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Upper_vs_Lower_DORCs_Volcano.png"), res = 600, height = 6, width = 8, units = "in") 
DORC.VolcanoPlot
dev.off()

#-------Reading in DA Peaks from ArchR-------#
COL17.DAPeaks <- read_csv(file = paste0(SeuratDirectory, "COL17_Upper_vs_Lower_DA_Peaks.csv"))
COL17.DAPeaks
#Getting Differential Peaks per Gene
COL17.UpperDAP <- COL17.DAPeaks %>% filter(Log2FC > 0)
COL17.UpperDAP
COL17.LowerDAP <- COL17.DAPeaks %>% filter(Log2FC < 0)
COL17.LowerDAP

#-------Getting Number of (same direction) DAPs per DA DORC-----#
#Upper DAPs
Upper.DORCs <- DORC.LowervsUpper %>% filter(avg_log2FC > 0) %>% pull(gene)
UpperDORC.SigLinks <- P2G500.SigLinks %>% filter(gene %in% Upper.DORCs)
#Getting Counts for All DORCs
UpperDORC.DALinks <- UpperDORC.SigLinks %>% filter(peak %in% COL17.UpperDAP$peak) %>% 
  group_by(gene) %>% dplyr::summarize(N_DAPs = n()) %>% arrange(desc(N_DAPs))
UpperDORC.DALinks
UpperDORC.NoLinks <- data.frame(gene = Upper.DORCs) %>% filter(gene %in% UpperDORC.DALinks$gene == F) %>% 
  mutate(N_DAPs = 0)
UpperDORC.NoLinks
#Lower DAPs
Lower.DORCs <- DORC.LowervsUpper %>% filter(avg_log2FC < 0) %>% pull(gene)
LowerDORC.SigLinks <- P2G500.SigLinks %>% filter(gene %in% Lower.DORCs)
LowerDORC.DALinks <- LowerDORC.SigLinks %>% filter(peak %in% COL17.LowerDAP$peak) %>% 
  group_by(gene) %>% dplyr::summarize(N_DAPs = n()) %>% arrange(desc(N_DAPs))
LowerDORC.DALinks
LowerDORC.NoLinks <- data.frame(gene = Lower.DORCs) %>% filter(gene %in% LowerDORC.DALinks$gene == F) %>% 
  mutate(N_DAPs = 0)
LowerDORC.NoLinks
#Concatenating
COL17.NDAPs <- as.data.frame(rbind(UpperDORC.DALinks, LowerDORC.DALinks, UpperDORC.NoLinks, LowerDORC.NoLinks))
rownames(COL17.NDAPs) <- COL17.NDAPs$gene
COL17.NDAPs
#Reordering to match DORC.LowervsUpper
COL17.NDAPs <- COL17.NDAPs[DORC.LowervsUpper$gene,]
sum(COL17.NDAPs$gene == DORC.LowervsUpper$gene)
#Adding Number of Linked Differentially Accessible Peaks per Gene to DETest
DORC.LowervsUpper$N_DAPs <- COL17.NDAPs$N_DAPs
DORC.LowervsUpper %>% arrange(desc(N_DAPs)) %>% head()

#-------Writing to csv-------#
write_csv(DORC.LowervsUpper, file = '/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_5/Plots/Upper_vs_Lower_DA_DORC_Genes.csv')

#-----Plotting Volcano Plot with N DAPs-----#
Top20Up <- DORC.LowervsUpper %>% filter(diffexp == "Upregulated") %>% arrange(p_val_adj) %>% head(30) %>% pull(gene)
Top20Down <- DORC.LowervsUpper %>% filter(diffexp == "Downregulated") %>% arrange(p_val_adj) %>% head(30) %>% pull(gene)
TopMarkers <- DORC.LowervsUpper %>% filter(gene %in% c(Top20Up, Top20Down))
SelectedMarkers <- DORC.LowervsUpper %>% filter(gene %in% c("DRAM1", "NTN1", "SPINK5", "SAMD5", "GATA3", "NFIB", "BACH2", "SHH", "FRY",
                                                            "CPA6", "KRT35", "KRT36", "SNRK", "TLE4", "DNASE1L2"))
#Adding NDAPs as size aesthetic
DORC.NDAPVolcanoPlot <- ggplot(DORC.LowervsUpper, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp, size = N_DAPs)) + geom_point(alpha = 0.75) + geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')+
  geom_vline(xintercept = c(-0.25, 0.25), col = "black", linetype = 'dashed') +
  ggrepel::geom_text_repel(data = TopMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), color = "black", size = 4, force = 1, max.overlaps = Inf) +
  theme_classic() +  scale_color_manual(values = c("cornflowerblue", "firebrick", "gray"), labels = c("Enriched in Upper COL17", "Enriched in Lower COL17", "Not Significant"))  + 
  labs(x = "Average Log2Fold Change", y = "- Log10 P Value", color = "Differential Regulation", size = "N DA Peaks", title = "Differential DORC Accessibility: Lower vs Upper COL17A1 Layer")
DORC.NDAPVolcanoPlot
#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Upper_vs_Lower_DORCs_NDAPs_Volcano.png"), res = 600, height = 6, width = 8, units = "in") 
DORC.NDAPVolcanoPlot
dev.off()

#-------Saving Plot without Labels-------#
#Adding NDAPs as size aesthetic
DORC.NDAPVolcanoPlotUnlabeled <- ggplot(DORC.LowervsUpper, aes(x = avg_log2FC, y = -log10(p_val_adj), color = diffexp, size = N_DAPs)) + geom_point(alpha = 0.75) + geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed')+
  geom_vline(xintercept = c(-0.25, 0.25), col = "black", linetype = 'dashed') + theme_classic() +  scale_color_manual(values = c("cornflowerblue", "firebrick", "gray"), labels = c("Enriched in Upper COL17", "Enriched in Lower COL17", "Not Significant"))  + 
  labs(x = "Average Log2Fold Change", y = "- Log10 P Value", color = "Differential Regulation", size = "N DA Peaks", title = "Differential DORC Accessibility: Lower vs Upper COL17A1 Layer")
DORC.NDAPVolcanoPlotUnlabeled

#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0('/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_5/Plots/', "COL17_Upper_vs_Lower_DORCs_NDAPs_Volcano_Unlabeled.png"), res = 600, height = 6, width = 8, units = "in") 
DORC.NDAPVolcanoPlotUnlabeled
dev.off()

#-------Finding Clusters of Peak Expression for each DA DORC Gene----------#
#------Calculating AUCell Module Scores for DAD genes-------#
COL17.SigDAGs <- DORC.LowervsUpper %>% filter(diffexp %in% c("Downregulated", "Upregulated"))
COL17.SigDAGs
#Aggregating Gene Expression Across Clusters
Matrix.Agg <- HHA.Matrix %>% JoinLayers() %>% 
  AggregateExpression(assays = "RNA", group.by = "MatrixAnnotationOrdered", normalization.method = "LogNormalize")
str(Matrix.Agg)
#Subsetting for DORC Genes
Matrix.AggDORC <- Matrix.Agg$RNA[COL17.SigDAGs$gene,]
str(Matrix.AggDORC)
colnames(Matrix.AggDORC) <- str_replace(colnames(Matrix.AggDORC), pattern = "-", replacement = "_")
#Function that returns top expressing clusters for each gene
GetTopExp <- function(AggMat) {
  TopClusters <- c()
  for (Gene in rownames(AggMat)) {
    print(Gene)
    SortExp <- rev(sort(AggMat[Gene,]))
    #Getting Top Cluster
    TopCluster <- names(SortExp)[1]
    #Adding to vector
    TopClusters <- c(TopClusters, TopCluster)}
  #returning as data frame
  TopClusters.DF <- data.frame(Gene = rownames(AggMat), TopCluster = TopClusters)
  rownames(TopClusters.DF) <- TopClusters.DF$Gene
  return(TopClusters.DF)}
#Retrieving Dataframe of Top Expressing Clusters
SigDAGs.TopExp <- GetTopExp(Matrix.AggDORC)
table(SigDAGs.TopExp$TopCluster)
#Adding to DORCs Lower vs Upper 
DORC.LowervsUpper <- DORC.LowervsUpper %>% mutate(TopCluster = factor(case_when(gene %in% SigDAGs.TopExp$Gene ~ SigDAGs.TopExp[gene,]$TopCluster,
                                                                                .default = "Not DE"), levels = c("Not DE", levels(Matrix.Multiome$MatrixAnnotationOrdered))))
table(DORC.LowervsUpper$TopCluster)

#-------Plotting Heatmap------#
#Aggregating Gene Expression Across Clusters and scaling
Matrix.DORCAggScaled <- HHA.Matrix %>% JoinLayers() %>% subset(MatrixAnnotationFine %in% c("Medulla", "LPC") ==F) %>% 
  AggregateExpression(assays = "RNA", group.by = "MatrixAnnotationOrdered", normalization.method = "LogNormalize", return.seurat = T) %>%
  ScaleData(features = COL17.SigDAGs$gene) %>% GetAssayData(assay = "RNA", slot = "scale.data")
colnames(Matrix.DORCAggScaled) <- str_replace_all(colnames(Matrix.DORCAggScaled), pattern = "-", replacement = "_")
#Reordering 
str(Matrix.DORCAggScaled)

#-----------Plotting Heatmap--------#
HeatmapOrder <- c("Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Henle", "IRS_Huxley", "IRS_Cuticle",
                  "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex","Early_Cuticle", "Middle_Cuticle", "Late_Cuticle")
HeatmapPal <- c("#89C75F", "#D8A767", "#F47D2B", "#F37B7D",  "#7E1416", "#D24B27", "#3BBCA8", "#9ECAE1", "#4292C6", "#08306B", "#E6C2DC", "#C06CAB",  "#89288F")
names(HeatmapPal) <- HeatmapOrder
GEPSplit <- c(rep("Lower", 6), rep("Upper", 7))
#Plotting
DORC.AggHeatmap <- Heatmap(Matrix.DORCAggScaled[COL17.SigDAGs$gene, HeatmapOrder], 
                           cluster_rows = T,
                           column_names_rot = 45,
                           cluster_columns = F,
                           show_row_names = T,
                           row_names_gp = gpar(fontsize = 3),
                           row_split = COL17.SigDAGs$diffexp,
                           column_split = GEPSplit,
                           #left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#00A08A", "#5BBCD6")))),
                           top_annotation =  HeatmapAnnotation(Cluster = HeatmapOrder, col = list(Cluster = HeatmapPal)),
                           heatmap_legend_param = list(title = "Z Score"),
                           col = colorRamp2(seq(from= -2,to=2,length=9), rev(RColorBrewer::brewer.pal(9, "RdBu"))),
                           row_dend_reorder = T)
draw(DORC.AggHeatmap, padding = unit(c(1, 1, 1, 1), "cm"))

#-------Saving as tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_COL17_DAG_Cluster_Heatmap_Large.tiff"), res = 300, height = 18, width = 24, units = "in") 
draw(DORC.AggHeatmap, padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()

#-------Saving as tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_COL17_DAG_Cluster_Heatmap.tiff"), res = 300, height = 6, width = 8, units = "in") 
draw(DORC.AggHeatmap, padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()

#------Calculating AUCell Module Scores for DAD genes-------#
Lower.DAGs <- DORC.LowervsUpper %>% filter(diffexp == "Downregulated") %>% pull(gene)
Upper.DAGs <- DORC.LowervsUpper %>% filter(diffexp == "Upregulated") %>% pull(gene)
COL17.DAGs <- list("Lower_DA" = Lower.DAGs,
                   "Upper_DA" = Upper.DAGs)
COL17.DAGs

#------Scoring DAGs in scRNA-seq data------#
#Running AUCell
COL17DAG.MultiomeScores <- RunAUCell(JoinLayers(Matrix.Multiome), COL17.DAGs)
#Adding to Seurat
Matrix.Multiome <- AddMetaData(Matrix.Multiome, COL17DAG.MultiomeScores)
Matrix.Multiome[[]]
#Plotting Violins
VlnPlot(Matrix.Multiome, features = c('Lower_DA', 'Upper_DA'), group.by = "MatrixAnnotationOrdered", cols = MatrixPal)
FeaturePlot(Matrix.Multiome, features = c('Lower_DA', 'Upper_DA'), reduction = "UMAP_WNN", pt.size = 1) & scale_color_viridis_c(option = "inferno")

#--------Retrieving Spatial Data from Matrix-------#
#Loading Full Spatial Object
load("/projects/b1217/HHA/Spatial/Seurat_Data/HHA_Integrated_Spatial_Initial_Annotations_12_11_24.rds")
#Reading in relevant barcodes
Spatial.Metadata <- read_csv('/projects/b1217/HHA/Bulb_Spatial/20250516_HHA_Spatial_TACCO_Matrix_Metadata_6_11_25.csv') %>% column_to_rownames("...1")
Spatial.Metadata #Metadata from spatial object after TACCO Mapping/subsetting for Matrix
#Extracting Cells with TACCO Mapping to Matrix Clusters
HHA.Merged <- JoinLayers(HHA.Merged)
#Subsetting for relevant sample
Spatial.Matrix <- subset(HHA.Merged, orig.ident == "EL_S2")
Spatial.Matrix$OriginalBarcode <- str_split_i(colnames(Spatial.Matrix), pattern = "_2", i = 1)
#Subsetting for cells annotated as Matrix 
Spatial.Matrix <- subset(Spatial.Matrix, OriginalBarcode %in% Spatial.Metadata$barcode)
colnames(Spatial.Matrix) <- Spatial.Matrix$OriginalBarcode
Spatial.Matrix #10,584 spots
rm(HHA.Merged)

#------Scoring DAGs in spatial data------#
#Running AUCell
COL17DAG.SpatialScores <- RunAUCell(Spatial.Matrix, COL17.DAGs, assay = "Spatial.008um")
#Adding to Seurat
Spatial.Matrix <- AddMetaData(Spatial.Matrix, COL17DAG.SpatialScores)
Spatial.Matrix[[]]
#Adding SpatialScores to Metadata
SpatialMetadata.COL17DAGScores <- cbind(Spatial.Metadata, COL17DAG.SpatialScores[rownames(Spatial.Metadata),])
SpatialMetadata.COL17DAGScores
#Plotting Violins
ggplot(SpatialMetadata.COL17DAGScores, aes(x = MatrixAnnotationFine, y = Lower_DA)) + geom_boxplot() & RotatedAxis()
#Plotting Violins
ggplot(SpatialMetadata.COL17DAGScores, aes(x = MatrixAnnotationFine, y = Upper_DA)) + geom_boxplot() & RotatedAxis()

#-------Saving Results--------#
write_csv(SpatialMetadata.COL17DAGScores, file = "/projects/b1217/HHA/Bulb_Spatial/HHA_Spatial_Lower_Upper_DA_AUCell_6_11_25.csv")

#--------Saving Image------#
save.image(paste0(SeuratDirectory, file = "HHA_Matrix_DORC_Scoring_1_7_26.RData"))














