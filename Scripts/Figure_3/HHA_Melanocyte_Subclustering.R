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
library(ggpubr)
library(rstatix)
library(lme4)
library(car)
library(emmeans)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Melanocyte_Plots/"
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
AnnotationPal <- c('#91921a', '#aeaeae', '#ad8fc3', '#103d5d', '#185a88', '#4a2d27',
                   '#86acdd', '#aec7e8', '#d6e2f3', '#76cbda', '#6b4239', '#8c564b',
                   '#aa6d60', '#bc8b81', '#ffa145', '#ffd5ab', '#1f77b4', '#ab1f20',
                   '#df5152', '#da4daf', '#9edae5', '#df62b9', '#57a9e2', '#217821',
                   '#37c837', '#f288b6', '#e377c2', '#e78ccb', '#eca1d5', '#bebebe',
                   '#1294a1', '#31d7e8', '#ddd1e7', '#c6e9f0', '#2b93db', '#ff6663',
                   '#bcbd22', '#dadb37', '#fce4ee', '#d0d0d0', '#e0e0e0')
Feathered <- c("#606372", "#79A8A4", "#B2AD8F", "#AD8082", "#DEC18C", "#92A185") 
#MoMAColors::Abbot + gray
Abbot <- c('#950404FF', '#E04B28FF', '#C38961FF', '#9F5630FF', '#388F30FF', '#0F542FFF', '#007D82FF', '#004042FF', "#606372")
Bay <- c('#00496FFF', '#0F85A0FF', '#EDD746FF', '#ED8B00FF', '#DD4124FF')
CafeDeNuit <- c('#467326FF', '#8AA676FF', '#D9B23DFF', '#BF793BFF', '#A63333FF')
EdvardMunch1 <- c('#272A2AFF', '#E69253FF', '#EDB931FF', '#E4502EFF', '#4378A0FF')

#-----------Loading Workspace------------#
load(paste0(WorkspaceDirectory, "HHA_Melanocyte_Subclustering_11_21_25.RData"))

#------Loading Full scRNA/Multiome Atlas Object------#
load(paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))
HHA.Full

#Plotting Clusters
DimPlot(HHA.Full, reduction = "UMAP", group.by = "FinalAnnotation", label = T, repel = T, raster = F, cols = AnnotationPal) + 
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()

#----------Pulling Melanocytes for Subclustering----------#
HHA.Melano <- subset(HHA.Full, FinalAnnotation %in% c("Melanocytes I", "Melanocytes II"))
HHA.Melano

#-----ReClustering in Original SCVI Space-----#
#Running a new SCVI model appears to be less effective at subclustering due to low cell counts
HHA.Melano <-FindNeighbors(HHA.Melano, dims = 1:50, reduction = "scvi")
HHA.Melano <- RunUMAP(HHA.Melano, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.3, cluster.name = "scvi_clusters_0.3")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.25, cluster.name = "scvi_clusters_0.25")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.2, cluster.name = "scvi_clusters_0.2")
HHA.Melano<- FindClusters(HHA.Melano, resolution = 0.1, cluster.name = "scvi_clusters_0.1")

#------Plotting UMAPs-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "FinalAnnotation", raster = F, pt.size = 3) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting Clusters-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.2", raster = F, pt.size = 3) + 
  labs(x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Sample ID-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = F, shuffle = T, group.by = "SampleID", raster = F,
        cols = MetBrewer::met.brewer("Signac", 16), pt.size = 3) + theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#-----Plotting by Donor ID-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = F, shuffle = T, group.by = "DonorID", raster = F,
        cols = MetBrewer::met.brewer("Signac", 9), pt.size = 3) + theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)

#-----Plotting by Isolation Method-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = F, shuffle = T, group.by = "Isolation", raster = F, pt.size = 3) + theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
HHA.Melano$Isolation <- factor(HHA.Melano$Isolation, levels = c("4mm", "FUE", "Bulb"))
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.2", split.by = "Isolation", raster = F, pt.size = 3) & NoLegend() & NoAxes()
#-----Plotting by Platform-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = F, shuffle = T, group.by = "Platform", raster = F, pt.size = 3) + theme(aspect.ratio = 1) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
#Split
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.2", split.by = "Platform", raster = F, pt.size = 3) & NoLegend() & NoAxes()
#-----Plotting by Sex-----#
#Split
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.2", split.by = "Sex", raster = F, pt.size = 3) & NoLegend() & NoAxes()
#-----Plotting by Phase-----#
DimPlot(HHA.Melano, reduction = "UMAP", label = T, group.by = "Phase", raster = F, pt.size = 3) & NoLegend() & NoAxes()

#-----Finding Markers for Each Cluster-----#
#Getting markers for each cluster 
MelanoMarkers.TFIDF <- HHA.Melano %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA.Melano$scvi_clusters_0.2, N = 20, FDR = 0.01)
MelanoMarkers.TFIDF
#converting to list
MelanoMarkers.TFIDFList <- MelanoMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(MelanoMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA.Melano$scvi_clusters_0.2))-1))
MelanoMarkers.TFIDFList

#-----------Plotting Top Markers-----------#
#Cluster 0 Markers 
FeaturePlot(HHA.Melano, features = c("MGP", "NTRK2", 'FMNL2', 'CREB5'), reduction = "UMAP", pt.size = 1) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Cluster 1 Markers 
FeaturePlot(HHA.Melano, features = c("TPPP3", "STXBP6", 'TFAP2B', 'PDE4B'), reduction = "UMAP", raster = F, pt.size = 1) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Cluster 2 Markers 
FeaturePlot(HHA.Melano, features = c("FSTL5", "LINC02099", 'NRG3', 'MET'), reduction = "UMAP", raster = F, pt.size = 1) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Differentiation Markers
FeaturePlot(HHA.Melano, features = c("DCT", "PMEL", "MITF", "MLANA","SLC45A2", "NES"), reduction = "UMAP", raster = F, pt.size = 1) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Differentiation Markers
FeaturePlot(HHA.Melano, features = c("SOX10", "DCT", "PAX3", "MITF", "MLANA", "PMEL", "TYR", "SLC25A3", "SLC45A2"), reduction = "UMAP", raster = F, pt.size = 1) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()

#--------Assigning Temporary Identities--------------#
#MSCs I and II are DCT, NFIB+ 
MelanoAnnotations <- c("Melanocytes_I", "Melanocytes_II", "Melanocytes_III")
#Fine Annotation
Idents(HHA.Melano) <- "scvi_clusters_0.2"
names(MelanoAnnotations) <- levels(HHA.Melano)
HHA.Melano <- RenameIdents(HHA.Melano, MelanoAnnotations)
HHA.Melano$MelanoAnnotation <- HHA.Melano@active.ident

#----------Plotting UMAP-----------#
#Plotting Clusters
MelanoPal <- c('#EDB931FF', '#E4502EFF', '#4378A0FF')
Melano.UMAP <- DimPlot(HHA.Melano, reduction = "UMAP", label = F, shuffle = T, group.by = "MelanoAnnotation", raster = F,
        cols = MelanoPal, pt.size = 2)  + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
Melano.UMAP

#-------Saving As PNG----------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melanocyte_Subclustering_UMAP.png"), res = 600, height = 6, width = 8, units = "in") 
Melano.UMAP
dev.off()

#-------------Plotting Melanocyte Differentiation Genes by Cluster----------#
#Plotting Markers
DiffViolin <- VlnPlot(HHA.Melano, features = c("SOX10", "KIT", "MITF", "DCT", "MLANA", "PMEL", "TYRP1", "TYR", "SLC45A2"),
        alpha = 0.1, group.by = "MelanoAnnotation", stack = T, flip = T, cols = Abbot) & 
  theme(aspect.ratio = 0.25) & NoLegend() & labs(x = "")
DiffViolin

#-------Saving As PNG----------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melanocyte_Differentiation_Gene_Violin.png"), res = 600, height = 6, width = 8, units = "in") 
DiffViolin
dev.off()

#------------Plotting Ratios by Covariate------------#
#Sex
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, fill = Sex)) + geom_bar(position = "fill") + theme_bw()
#Platform
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, fill = Platform)) + geom_bar(position = "fill") + theme_bw()
#Final Annotation versus Melanocyte annotation
ggplot(HHA.Melano[[]], aes(x = FinalAnnotation, fill = MelanoAnnotation)) + geom_bar(position = "fill") + theme_bw() + scale_fill_manual(values = MelanoPal)
#Final Annotation versus Melanocyte annotation
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, fill = FinalAnnotation)) + geom_bar(position = "fill") + theme_bw()
#Annotated Cell Cycle Phase
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, fill = Phase)) + geom_bar(position = "fill") + theme_bw()

#----------Plotting Ratios by Isolation Method-----------#
IsolationPlot <- ggplot(HHA.Melano[[]], aes(x = Isolation, fill = MelanoAnnotation)) + geom_bar(position = "fill") + 
  theme_bw()  + scale_fill_manual(values = MelanoPal) + labs(x = "Isolation Method", y = "Cell Type Proportion", fill = "Population")
IsolationPlot

#-------Saving As PNG----------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melanocyte_Isolation_Proportion_Barplot.png"), res = 600, height = 6, width = 8, units = "in") 
IsolationPlot
dev.off()

#--------------Plotting Gene Expression Heatmap-------------#
HHA.Melano <-  ScaleData(HHA.Melano, features = rownames(HHA.Melano))
#Pulling Top 10 Marker Genes by Cluster
Melano.HeatMarkers <- MelanoMarkers.TFIDFList %>% map(~(head(.x, 10))) %>% bind_rows() %>% pull(gene) 
Melano.HeatMarkers
Melano.Heatmap <- DoHeatmap(HHA.Melano, features = Melano.HeatMarkers, group.by = "MelanoAnnotation", group.colors = MelanoPal) + 
  NoLegend() + scale_fill_viridis_c(option = "cividis") + theme(legend.position = "right")
Melano.Heatmap

#------Saving as tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melano_Heatmap.png"), res = 600, height = 8, width = 12, units = "in") 
Melano.Heatmap
dev.off()

#----------------Plotting Gene Expression Violins---------------#
#Plotting Marker Violin Plot: Vertical
Melano.MarkerVlnPlot <- VlnPlot(HHA.Melano, group.by = "MelanoAnnotation", features = c("TPPP3", "TFAP2B", "STXBP6", "MGP", "NPAS3", "FMNL2", "FSTL5", "NRG3", "MET"),
        fill.by = "ident", cols = MelanoPal, stack = T, flip = T) & theme(aspect.ratio =  0.25)
Melano.MarkerVlnPlot
#Plotting Marker Violin Plot: Horizontal
Melano.MarkerVlnPlotFlipped <-VlnPlot(HHA.Melano, group.by = "MelanoAnnotation", features = c("TPPP3", "TFAP2B", "STXBP6", "MGP", "NPAS3", "FMNL2", "FSTL5", "NRG3", "MET"),
        fill.by = "ident", cols = MelanoPal, stack = T, flip = F) & theme(aspect.ratio =  4)
Melano.MarkerVlnPlotFlipped

#------Saving as tiff-----#
#Vertical 
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melano_Marker_Violins.png"), res = 600, height = 12, width = 8, units = "in") 
Melano.MarkerVlnPlot
dev.off()
#Horizontal
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melano_Marker_Violins_Flipped.png"), res = 600, height = 12, width = 8, units = "in") 
Melano.MarkerVlnPlot
dev.off()

#-------------Plotting Melanocytes on Original UMAP------------#
Melano.Barcodes <- HHA.Melano[[]] %>% rownames_to_column("Barcode") %>% select(Barcode, MelanoAnnotation)
Melano.Barcodes
#Adding Barcodes to original object
HHA.Full$Barcode <- colnames(HHA.Full)
#Subsetting out Melanocytes again 
Melano.OrigUMAP <- HHA.Full %>% JoinLayers() %>% subset(FinalAnnotation %in% c("Melanocytes I", "Melanocytes II"))
Melano.OrigUMAP
#Adding Barcode Map
Melano.BM <- Melano.Barcodes %>% column_to_rownames("Barcode")
Melano.BM <- Melano.BM[colnames(HHA.Melano),]
#Reordering to match the original ordering
Melano.OrigUMAP$MelanoAnnotation  <- Melano.BM

#-------Plotting Melanocyte Annotations on Original UMAP----------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melanocyte_Subclustering_Original_UMAP.png"), res = 600, height = 6, width = 8, units = "in") 
DimPlot(Melano.OrigUMAP, group.by = "MelanoAnnotation", reduction = "UMAP", cols = MelanoPal)  + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
dev.off()

#--------Pigmentation/Proliferation Analysis----------#
#Grabbing Seurat Cell Cycle Genes
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Reading in Gene Set (from Science Bajpai et al, 2023, A genome-wide genetic screen uncovers determinants of human pigmentation)
PigGenes <- readxl::read_xlsx(paste0(MetadataDirectory, "science.ade6289_table_s1.xlsx"))
#Subsetting for genes significantly associated with pigmentation using qval from paper 
PigGenes.use <- PigGenes %>% filter(q_value < 0.1) %>% pull(Symbol)
PigGenes.use #169 genes associated with pigmentation 

#---------Running AUCell------------#
#Calculating stress response score with AUCell
Melano.PigmentAUC <- HHA.Melano %>% JoinLayers() %>% RunAUCell(list(Pigmentation_AUCell = PigGenes.use, Proliferation_AUCell = CellCycleGenes))
Melano.PigmentAUC
#Adding to Seurat Object
HHA.Melano$PigmentationScore <- Melano.PigmentAUC$Pigmentation_AUCell
HHA.Melano$ProliferationScore <- Melano.PigmentAUC$Proliferation_AUCell
#Visualizing by Cluster with violin plot
VlnPlot(HHA.Melano, features = "PigmentationScore", group.by = "MelanoAnnotation", cols = MelanoPal, alpha = 0, raster = F)
#Visualizing by Melanocytes I vs II
VlnPlot(HHA.Melano, features = "PigmentationScore", group.by = "FinalAnnotation", alpha = 0, raster = F)
VlnPlot(HHA.Melano, features = "ProliferationScore", group.by = "MelanoAnnotation", alpha = 0, raster = F)

#--------------Plotting Results---------#
#Plotting on Feature Plot
FeaturePlot(HHA.Melano, features = c("PigmentationScore"), reduction = "UMAP", pt.size = 3) & 
  scale_color_gradientn(colors = GrayMagma) & NoAxes()
#Plotting as BoxPlot
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, y = PigmentationScore, color = MelanoAnnotation)) + 
  theme_bw() + geom_boxplot() + scale_color_manual(values = MelanoPal) + geom_point(position = "jitter", alpha = 0.25)
#Plotting as BoxPlot
ggplot(HHA.Melano[[]], aes(x = FinalAnnotation, y = PigmentationScore, color = FinalAnnotation)) + 
  theme_bw() + geom_boxplot() + scale_color_manual(values = MelanoPal) + geom_point(position = "jitter", alpha = 0.25)
#Plotting as BoxPlot
ggplot(HHA.Melano[[]], aes(x = MelanoAnnotation, y = ProliferationScore, color = MelanoAnnotation)) + 
  theme_bw() + geom_boxplot() + scale_color_manual(values = MelanoPal) + geom_point(position = "jitter", alpha = 0.25)

#---------Checking for Platform Specific Effects------------#
PlatformPal <- c("gray40", "#E4502EFF", "#4378A0FF")
#Multiome data shows consistently lower scores versus the other two platforms. 
ggboxplot(HHA.Melano[[]], x = "MelanoAnnotation", y = "PigmentationScore", color = "Platform", legend = "right", 
          palette = PlatformPal, add = "jitter", add.params = list(alpha = 0.6, size = 0.35)) +
  labs(x = "", y = "Seurat G2M Score") + theme(aspect.ratio = 2, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis()

#------------Computing Statistics: Pigmentation--------------#
V4.PigmentationScores <- HHA.Melano[[]] %>% filter(Platform == "V4") %>% select(MelanoAnnotation, PigmentationScore)
V4.PigmentationScores
#Computing Kruskall-Wallis Test
Pigmentation.KW <- kruskal_test(V4.PigmentationScores, PigmentationScore ~ MelanoAnnotation)
Pigmentation.KW
#Writing to csv
write_csv(Pigmentation.KW, file = paste0(PlotsDirectory, "HHA_Melano_Pigmentation_Score_Kruskall_Wallis.csv"))
#Computing Pairwise Wilcoxon
Pigmentation.Wilcox <- V4.PigmentationScores %>%
  pairwise_wilcox_test(PigmentationScore ~ MelanoAnnotation) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")
Pigmentation.Wilcox 
#Writing to csv
write_csv(Pigmentation.Wilcox, file = paste0(PlotsDirectory, "HHA_Melano_Pigmentation_Score_Wilcoxon.csv"))

#------------Plotting BoxPlots: Pigmentation--------------#
Pigmentation.BoxPlot <- ggboxplot(V4.PigmentationScores, x = "MelanoAnnotation", y = "PigmentationScore", fill = "MelanoAnnotation", legend = "none", palette = MelanoPal, add = "jitter",
                                  add.params = list(alpha = 0.1, size = 0.35)) + labs(x = "", y = "Pigmentation Score", fill = "Annotation") + 
  theme(aspect.ratio = 1.2, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis() +
  guides(color = guide_legend(override.aes = list(shape = 16,size  = 4,alpha = 1, linetype = 0))) 
Pigmentation.BoxPlot

#----------------Plotting with Statistics--------------#
#Adding Y Positions for Plotting 
Pigmentation.WilcoxPlot <- Pigmentation.Wilcox %>% mutate(y.position = max(V4.PigmentationScores$PigmentationScore) + row_number() * 0.03) 
#Plotting with Stats
Pigmentation.StatsPlot <- ggboxplot(V4.PigmentationScores, x = "MelanoAnnotation", y = "PigmentationScore", fill = "MelanoAnnotation", legend = "none", palette = MelanoPal, add = "jitter",
                                    add.params = list(alpha = 0.2, size = 0.35)) + labs(x = "", y = "Pigmentation Score", fill = "Annotation") + # Add pairwise comparisons p-value
  stat_pvalue_manual(Pigmentation.WilcoxPlot,label = "p.adj.signif", y.position = "y.position", tip.length = 0.01, hide.ns = F) + theme_classic() +
  theme(aspect.ratio = 1.2, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis()
Pigmentation.StatsPlot 

#------Saving as PNGs-----#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melano_Pigmentation_Boxplot_V4.png"), res = 600, height = 6, width = 8, units = "in") 
Pigmentation.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melano_Pigmentation_Boxplot_V4_Stats.png"), res = 600, height = 6, width = 8, units = "in") 
Pigmentation.StatsPlot
dev.off()

#------------DEG Analysis between Clusters---------#
#1 Vs All 
MelanoMarkers.wilcox <- HHA.Melano %>% JoinLayers() %>% FindAllMarkers(test.use = "wilcox", only.pos = T)
#Splitting to List
MelanoMarkers.wilcoxList <- MelanoMarkers.wilcox %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>% group_by(cluster) %>% group_split()
names(MelanoMarkers.wilcoxList) <- c("Melanocytes_I", "Melanocytes II", "Melanocytes III")
#Runig GO Term Aalysis
MelanoMarkers.wilcoxList
MelanoMarkers.GeneList <- MelanoMarkers.wilcoxList %>% map(~.x %>% pull(gene))

#----------Running GO Term Analysis for Melanocyte Markers-------#
MelanoMarkers.GO <- map(MelanoMarkers.GeneList, ~enrichGO(.x, OrgDb = "org.Hs.eg.db", keyType = 'SYMBOL', readable = T,
                                                          ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH"))
names(MelanoMarkers.GO) <- names(MelanoMarkers.GeneList)

#-------Plotting GO Term Dot Plots for Each Marker Group-----#
#Module 1
dotplot(MelanoMarkers.GO[[1]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "p.adjust", y = "Adjusted P Value") & ggtitle("Melanocytes 1")
#Module 2
dotplot(MelanoMarkers.GO[[2]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Melanocytes 2")
#Module 3
dotplot(MelanoMarkers.GO[[3]], x = "p.adjust", color = "p.adjust", showCategory = 12)  & 
  scale_fill_viridis_c(option = "mako") & labs(x = "Gene Ratio", y = "Adjusted P Value") & ggtitle("Melanocytes 3")

#--------------Exporting to Python--------------#
#------------Exporting Barcodes------------#
Melano.Barcodes <- HHA.Melano[[]] %>% rownames_to_column("Barcode") %>% select(Barcode, MelanoAnnotation)
Melano.Barcodes
#Writing to csv
write_csv(Melano.Barcodes, file = paste0(MetadataDirectory, "HHA_Melanocyte_Annotations_11_18_25.csv"))

#----------Exporting Barcodes for Melanocytes/Schwann Cells---------#
Schwann.Barcodes <- HHA.Full[[]] %>% filter(FinalAnnotation == "Schwann") %>% rownames_to_column("Barcode") %>% select(Barcode, MelanoAnnotation = FinalAnnotation)
Schwann.Barcodes
#Concatenating with Melanocyte Barcodes
MelanoSchwann.Barcodes <- rbind(Melano.Barcodes, Schwann.Barcodes)
MelanoSchwann.Barcodes
#Writing to csv
write_csv(MelanoSchwann.Barcodes, file = paste0(MetadataDirectory, "HHA_Melanocyte_Schwann_Annotations_11_18_25.csv"))

#--------Exporting to Python for Downstream Analysis------#
AnnDataDirectory <- "/projects/b1217/HHA/Multiome_Scanpy_Conversion/Full_Atlas/"
#------Full Conversion to Scanpy------#
#----Convert seurat object to anndata object----#
#----Writing Metadata----#
Melano.metadata <- HHA.Melano[[]] %>% mutate(UMAP_1 = HHA.Melano@reductions$UMAP@cell.embeddings[,1], UMAP_2 = HHA.Melano@reductions$UMAP@cell.embeddings[,2])
Melano.metadata$barcode <- rownames(Melano.metadata)
Melano.metadata <- Melano.metadata
write.csv(Melano.metadata, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Melano_FullConversion_Metadata_11_18_25.csv"), quote = F, row.names = F)
#----Writing scVI Embeddings----#
write.csv(HHA.Melano@reductions$scvi@cell.embeddings, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Melano_FullConversion_scVI_Embeddings_11_18_25.csv"), quote = F, row.names = F)

#------Exporting SNN Graph Connectivities------#
#Getting long format DF containing all possible connections and their presence/absence
Melano.SNN <- as.data.frame(summary(HHA.Melano@graphs$RNA_nn)) %>% #Matching Pythonic Indexing (starting at 0)
  mutate(i = i-1, j = j-1)
write_csv(Melano.SNN, file = paste0(AnnDataDirectory, "HHA_Multiome_Melano_SNN_Connectivities_5_22_25.csv"))

#----Writing Count Matrix and Gene Names----#
#Pulling Normalized Counts Data
Melano.logcounts <- HHA.Melano %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'data')
Matrix::writeMM(Melano.logcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Melano_FullConversion_Expression_LogCounts_11_18_25.mtx"))
#Pulling Raw Counts Data
Melano.rawcounts <- HHA.Melano %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
Matrix::writeMM(Melano.rawcounts, file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Melano_FullConversion_RawCounts_11_18_25.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(Melano.logcounts)), file = paste0(AnnDataDirectory, "HHA_SCRNA_Multiome_Melano_FullConversion_GeneNames_11_18_25.csv"),
            quote = F, row.names = F, col.names = T)

#---------------Analyzing Melanocyte Localization--------------#
#B2C Data
Melano.SpatialMetadata <- read_csv("/projects/b1217/HHA/Melanocyte_Plots/b2c_Melano_Metadata.csv") %>% column_to_rownames("Barcode")
Melano.SpatialMetadata #Eccrine SG mappings are nonsense, due to high proportion of Schwann cells in Eccrine glands, + secretory environment. They lack melanocyte markers. 
#Plotting Distribution of Isolations
Melano.LocalizationPlot <- Melano.SpatialMetadata %>% filter(Localization != "Eccrine_SG") %>% 
  ggplot(aes(x = Localization, fill = MelanoAnnotation)) + geom_bar(position = "fill") + 
  theme_bw()  + scale_fill_manual(values = MelanoPal) + labs(x = "Localization Site", y = "Cell Type Proportion", fill = "Population")
Melano.LocalizationPlot 

#-------Saving As PNG----------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Melanocyte_Mapping_Proportion_Barplot.png"), res = 600, height = 6, width = 8, units = "in") 
Melano.LocalizationPlot 
dev.off()

#----------Computing Population Proportions per Localization Category-------------#
Melano.LocalizationRates <- Melano.SpatialMetadata %>%
  group_by(Localization, MelanoAnnotation) %>%
  summarize(NSegmentations = n(), .groups = "drop_last") %>%
  mutate(Proportion = NSegmentations / sum(NSegmentations)) %>% ungroup()
Melano.LocalizationRates
#Writing to csv
write_csv(Melano.LocalizationRates, file = paste0(PlotsDirectory, "Melanocyte_Proportions_IFE_Bulb.csv"))

#---------------Saving Image------------#
save.image(paste0(WorkspaceDirectory, "HHA_Melanocyte_Subclustering_11_21_25.RData"))





