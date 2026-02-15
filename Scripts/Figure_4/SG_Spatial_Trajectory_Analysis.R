#In this script, we use TradeSeq to calculate gene trends across pseudotime using a NB GAM model, 
#identify pseudospatially variable genes, and group them into functional modules using Louvain clustering. 

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
library(tradeSeq)
library(ComplexHeatmap)
library(circlize)
library(grid)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Spatial/Seurat_Data/"
PlotsDirectory <- "/projects/b1217/HHA/Spatial/Plots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
AnnDataDirectory <- "/projects/b1217/HHA/SG_Trajectory/AnnData/"

#-------Helper Functions-----#
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

#------Beautified Heatmap Visualization------#
#Custom function for Ordering Genes by timepoint of max expression. 
OrderGenes <- function(Smoothed) {
  #Getting max expression values per gene
  MaxExp <- rowMax(Smoothed)
  names(MaxExp) <- rownames(Smoothed)
  #Getting row indices for max expression values
  MaxIndices <- c()
  for (i in 1:length(MaxExp)) {
    #Pulling smoothed expression into vector
    SmoothedExp <- Smoothed[i, ]
    MaxIndex <- as.numeric(str_split_i(names(SmoothedExp[SmoothedExp == MaxExp[i]]), pattern = "_", i = 2))
    MaxIndices <- c(MaxIndices, MaxIndex)}
  Sorting.DF <- data.frame(Gene = rownames(Smoothed), MaxIndex = MaxIndices)
  GeneOrder <- Sorting.DF %>% arrange(MaxIndex) %>% pull(Gene)
  return(GeneOrder)}

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

#-----Loading Checkpoint-----#
load(paste0(WorkspaceDirectory, "SG_Spatial_Trajectory_Analysis_5_27_25.RData"))

#-------Loading Seurat Object--------#
load(paste0(SeuratDirectory, "HHA_Spatial_SG_4_30_25.rds"))

#------Splitting Seurat Objects by Cluster-----#
#Remapping Cluster names
Spatial.SG$SGAnnotation <- plyr::mapvalues(Spatial.SG$SGAnnotationOrdered, from = levels(Spatial.SG$SGAnnotationOrdered),
                                           to = c("Duct_Basal", "Peripheral Zone", "SEB_I", "SEB_II", "SEB_III", "SEB_IV",
                                                  "Duct_Suprabasal"))
#Plotting UMAP
DimPlot(Spatial.SG, reduction = 'UMAP', group.by = "SGAnnotation")

#-------Pulling Pseudotime Values from csv-------#
PseudotimeDir <- "/projects/b1217/HHA/SG_Trajectory/Pseudotime_Duct_SG/Palantir/"
#SG. 
SG.Pseudotime <- read_csv("/projects/b1217/HHA/SG_Trajectory/Palantir/HHA_SG_palantir_pseudotime_results_4_30_25.csv")

#------Removing poor quality spots removed in original pseudotime analysis and adding in pseudotime values -----#
Spatial.SG$barcode <- colnames(Spatial.SG)
Spatial.SG <- subset(Spatial.SG, barcode %in% SG.Pseudotime$barcode)
DimPlot(Spatial.SG, reduction = 'UMAP', group.by = "SGAnnotation")
#Reordering to match object
SG.Pseudotime <- full_join(data.frame(barcode = colnames(Spatial.SG)), SG.Pseudotime, by = "barcode")
#Adding in Pseudotime
Spatial.SG$SGPseudotime <- SG.Pseudotime$PalantirPseudotime
#Plotting UMAP
FeaturePlot(Spatial.SG, reduction = 'UMAP', features = "PalantirPseudotime") & scale_color_viridis_c(option = "plasma")

#-------Getting Spatially Variable Features-----#
SG.DuctTraj <- subset(Spatial.SG, SGAnnotation != "Duct_Basal")
SG.DuctTraj
#Visualizing Results: Clusters
DimPlot(SG.DuctTraj, reduction = "UMAP", label = F, group.by = "SGAnnotation", cols = c("#2a7185", "#a64027", "#fbdf72","#60824f", "#022336", "#725ca5")) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "") & theme(aspect.ratio = 1) & NoAxes()
#Visualizing Results: Pseudotime
FeaturePlot(SG.DuctTraj, reduction = "UMAP", features = "PalantirPseudotime") + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "") & theme(aspect.ratio = 1) & NoAxes() & scale_color_viridis_c()

#------Setting Up TradeSeq Run------##Setting all cell weights to be equal
SG.CellWeights <- rep(1, length(colnames(SG.DuctTraj))) %>% as.matrix()
rownames(SG.CellWeights) <- colnames(SG.DuctTraj)
#Extracting Pseudotime Values
SG.Pseudotime <- as.matrix(SG.DuctTraj$PalantirPseudotime)
rownames(SG.Pseudotime) <- colnames(SG.DuctTraj)

#-------Finding Expression Thresholds for GAM calculation------#
#Getting raw count matrix from Seurat Object
SG.counts <- SG.DuctTraj %>% JoinLayers() %>% GetAssayData(assay = "Spatial.008um", slot = "counts")
str(SG.counts)
SG.counts[1:10, 1:10]
#Calculating percentage cells with expression > 0 for each gene
SG.PercentExpressed <- data.frame(gene = rownames(SG.counts), PercentExpressing = rowMeans(SG.counts >0))
SG.PercentExpressed %>% arrange(desc(PercentExpressing)) %>% head(20)
#Plotting Histogram
ggplot(data.frame(SG.PercentExpressed, PercentExpressing = SG.PercentExpressed), aes(x = PercentExpressing)) +
  geom_histogram(bins = 100) + geom_vline(xintercept = 0.05, linetype = "dashed") + theme_classic()
#-------Filtering Expression Matrix------#
#This seems like a reasonable threshold due to noise in spatial data. 
GenesKept <- rowMeans(SG.counts >0) > 0.05
SG.filteredcounts <- SG.counts[GenesKept, ]
str(SG.filteredcounts) #2,969 genes kept after filtering

#-------Determining optimal number of knots for model----------#
#------Setting up Parallelization-----#
BiocParallel::register(BiocParallel::MulticoreParam(workers = 16))
BiocParallel::bpparam() #using 16 workers in total
set.seed(1701)
#Testing number of knots across 500 randomly selected genes. 
SG.icMat <- evaluateK(counts = SG.filteredcounts, pseudotime = SG.Pseudotime, cellWeights = SG.CellWeights, k = 3:15, 
                      family = "nb", nGenes = 500, verbose = T, parallel = T, BPPARAM = BiocParallel::bpparam()) #results plateau at 13 knots

#-------Fitting GAM for each gene-------#
#Fitting GAM to each gene using 13 knots
set.seed(1701)
#Can't parallelize due to errors in package 
Models <- fitGAM(SG.filteredcounts, pseudotime = SG.Pseudotime, cellWeights = SG.CellWeights, nknots = 13,
                 verbose = T, parallel = F)
Models
#Checking for Convergence
table(rowData(Models)$tradeSeq$converged) #All genes converged

#-----Saving Checkpoint-----#
save.image(paste0(WorkspaceDirectory, "SG_Spatial_Trajectory_Analysis_5_27_25.RData"))

#------Calculating Association Test--------#
SG.assoRes <- associationTest(Models, inverse="eigen")
head(SG.assoRes)
#Calculating Bonferroni Adjusted P Values
SG.assoRes$pvalue_adj <- p.adjust(SG.assoRes$pvalue, method = "bonferroni")
SG.assoRes %>% arrange(desc(waldStat)) %>% head(10)
table(SG.assoRes$pvalue_adj < 0.01)

#-------Plotting Individual Genes-----#
plotSmoothers(Models, SG.filteredcounts, gene = "KRT79")
plotSmoothers(Models, SG.filteredcounts, gene = "VIM")
plotSmoothers(Models, SG.filteredcounts, gene = "PPARG")
plotSmoothers(Models, SG.filteredcounts, gene = "KRT5")

#-------Getting Top 1500 Genes by Wald Statistic-------#
TrajGenes <- SG.assoRes %>% arrange(desc(waldStat)) %>% head(1500) %>% rownames()
TrajGenes

#-----Calculating Smoothed Models for each of the top 1500 genes-----#
set.seed(1701)
TrajGenes.Smoothed <- predictSmooth(Models, gene = TrajGenes, nPoints = 500, tidy = F)
#Scaling to Unit Variance
TrajGenes.SmoothedScaled <- t(scale(t(TrajGenes.Smoothed)))
TrajGenes.Order <- OrderGenes(TrajGenes.Smoothed)

#-------Louvain Clustering of Trajectory Associated Genes-------#
#Creating Seurat Object for Running Louvain Clustering 
Pseudotime.Seurat <- CreateSeuratObject(t(TrajGenes.SmoothedScaled))
#Adding the smoothed gene expression values as a dimensional reduction
Pseudotime.Seurat[["SmoothedExp"]] <- CreateDimReducObject(TrajGenes.SmoothedScaled, key = "lineage1_")
Pseudotime.Seurat@assays$RNA@layers$scale.data <- CreateAssayObject(t(TrajGenes.SmoothedScaled))
#Finding Neighbors
Pseudotime.Seurat <- FindNeighbors(Pseudotime.Seurat, reduction = "SmoothedExp", dims = 1:500)
#Running UMAP
Pseudotime.Seurat <- RunUMAP(Pseudotime.Seurat, reduction = "SmoothedExp", reduction.name = "UMAP_Smooth", dims = 1:500)
#Running Louvain Clustering
Pseudotime.Seurat <- FindClusters(Pseudotime.Seurat, cluster.name = "GeneModule", resolution = 0.3)
DimPlot(Pseudotime.Seurat, reduction = "UMAP_Smooth", group.by = "GeneModule", label = T)

LouvainClusters <- factor(plyr::mapvalues(Pseudotime.Seurat$GeneModule, from = levels(factor(Pseudotime.Seurat$GeneModule, levels = c(1, 4, 0, 2, 3))),
                                          to = c("Peripheral", "Early Maturation", "Mid Maturation", "Late Maturation", "Duct")), levels = c("Peripheral", "Early Maturation", "Mid Maturation", "Late Maturation", "Duct"))

#-----Plotting Heatmap-----#
set.seed(1701)
#Generating Heatmap
SG.PseudotimeHeatmap <- Heatmap(TrajGenes.SmoothedScaled, name = "z-score",
                                show_row_names = TRUE,
                                show_column_names = FALSE,
                                row_names_gp = gpar(fontsize = 0),
                                row_title_rot = 0,
                                cluster_rows = F,
                                cluster_row_slices = F,
                                row_split = LouvainClusters,
                                row_order = TrajGenes.Order,
                                cluster_columns = FALSE,
                                col = colorRamp2(seq(from=-2.5,to=2.5,length=9), solarExtra),
                                left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")))),
                                show_row_dend = F)
print(SG.PseudotimeHeatmap)

#-------Plotting as tiff--------#
ModuleDirectory <- "/projects/b1217/HHA/SG_Trajectory/Gene_Modules/"
#Plotting as tiff
pdf(file = paste0(ModuleDirectory,  "SG_Gene_Modules_Heatmap.pdf"),
               width = 8, height = 7)
print(SG.PseudotimeHeatmap)
dev.off()


#------Gene Module Analysis-------#
#Getting gene modules in data frame
GeneModules.DF <- data.frame(gene = rownames(TrajGenes.SmoothedScaled), GeneModule = LouvainClusters)
#Adding Entrez IDs for each gene
GeneModules.EntrezIDs <- bitr(GeneModules.DF$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#Mapping Back into Dataframe
GeneModules.DF$EntrezID <- plyr::mapvalues(GeneModules.DF$gene, from = GeneModules.EntrezIDs$SYMBOL, to = GeneModules.EntrezIDs$ENTREZID)
#Getting list of genes in each module
#Gene Symbols
GeneModules.List <- GeneModules.DF %>% group_by(GeneModule) %>% group_split() %>% map(~.x %>% pull(gene))
names(GeneModules.List) <- levels(LouvainClusters)
#ENTREZ IDs
GeneModulesEntrez.List <- GeneModules.DF %>% filter(EntrezID %in% GeneModules.EntrezIDs$ENTREZID) %>% group_by(GeneModule) %>% group_split() %>% map(~.x %>% pull(EntrezID))
names(GeneModulesEntrez.List) <- levels(LouvainClusters)
GeneModules.List

#---------Saving Gene Modules to csv file-----#
write_csv(GeneModules.DF, file = paste0(ModuleDirectory, "SG_Trajectory_GeneModules_5_27_25.csv"))

#-------Getting Enriched GO Terms for each module-----#
GeneModules.GO <- map(GeneModules.List, ~enrichGO(.x, OrgDb = "org.Hs.eg.db", keyType = 'SYMBOL', readable = T,
                                                  ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH"))
names(GeneModules.GO) <- names(GeneModules.List)
View(GeneModules.GO[[5]]@result %>% filter(p.adjust < 0.05))

#-------Writing GO Term enrichment results to csv------#
for (i in 1:length(GeneModules.GO)) {
  write_csv(GeneModules.GO[[i]]@result, file = paste0(ModuleDirectory, "SG_Trajectory_Module_", 
                                                      str_replace_all(names(GeneModules.GO[i]), pattern = " ", replacement = "_"),
                                                      "_GO_Results_5_27_25.csv"))}

#-------Plotting GO Term Dot Plots for Each Module-----#
#Module 1
dotplot(GeneModules.GO[[1]], x = "p.adjust", color = "p.adjust", showCategory = 10)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "p.adjust", y = "Adjusted P Value") & ggtitle("Module 1")
#Module 2
dotplot(GeneModules.GO[[2]], x = "p.adjust", color = "p.adjust", showCategory = 10)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Module 2")
#Module 3
dotplot(GeneModules.GO[[3]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Module 3")
#Module 4
dotplot(GeneModules.GO[[4]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Module 4")
#Suprabasal II Markers
dotplot(GeneModules.GO[[5]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Module 5")

#-------Scoring Modules-----#
#Joining Layers
Spatial.SG <- JoinLayers(Spatial.SG)
#Calculating module scores with AUCell
GeneModules.SpatialScores <- RunAUCell(Spatial.SG, GeneModules.List)
colnames(GeneModules.SpatialScores) <- paste0("AUCell_", colnames(GeneModules.SpatialScores))
#Adding to Object
Spatial.SG <- AddMetaData(Spatial.SG, GeneModules.SpatialScores)
#Plotting on UMAP
FeaturePlot(Spatial.SG, features = colnames(GeneModules.SpatialScores), reduction = "UMAP") & 
  scale_color_viridis_c(option = "magma")
#Plotting ViolinPlot
VlnPlot(Spatial.SG, features = colnames(GeneModules.SpatialScores), group.by = "SGAnnotationOrdered", stack =T, flip = T)

#-------Writing Module Scores to csv------#
write_csv(GeneModules.SpatialScores, file = paste0(ModuleDirectory, "SG_Module_AUCell_Scores_5_27_25.csv"))

#-----Saving Checkpoint-----#
save.image(paste0(WorkspaceDirectory, "SG_Spatial_Trajectory_Analysis_5_27_25.RData"))

read_csv(paste0(ModuleDirectory, "SG_Module_AUCell_Scores_5_27_25.csv"))