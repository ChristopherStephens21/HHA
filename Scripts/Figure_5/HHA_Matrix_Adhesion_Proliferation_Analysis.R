#------Pre-Processing and Analysis of HHA Multiomics Data-----#
#This script contains our workflow for initial processing and integration of our single cell multiomics data. 
#We will generate multimodal Seurat objects, perform QC, run scDblFinder based on the RNA data and then integrate with the scRNA-seq data to get clusters. 
#-----Loading Packages-----#
library(Seurat)
library(scCustomize)
library(scDblFinder)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(AUCell)
library(clusterProfiler)
library(CellChat)
#Stats
library(ggpubr)
library(rstatix)

#Avoiding memory error
options(future.globals.maxSize = 3e+09)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Revisions/Matrix_Adhesion_Proliferation_Plots/"
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
MatrixPal <- c("#89C75F", "#3BBCA8",  "#208A42", "#0C727C", "#9ECAE1", "#4292C6", "#08306B", "#E6C2DC", "#C06CAB",  "#89288F", "#D8A767", "#F47D2B", "#F37B7D",  "#7E1416", "#D24B27")

#---------Loading Checkpoint-------#
load(paste0(SeuratDirectory, "HHA_Matrix_Adhesion_Proliferation_Analysis_1_26_26.RData"))

#--------Loading Full HHA Matrix Object------#
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))
HHA.Matrix

#-----Loading Global Object-----#
load(paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))
HHA.Full

#-------Selecting DP for Analysis------#
HHA.DP <- HHA.Full %>% JoinLayers() %>% subset(GlobalAnnotation == "Dermal_Papilla")
HHA.DP
#-----Merging with Matrix-----#
HHA.Matrix <- JoinLayers(HHA.Matrix)
HHA.MatrixDP <- merge(HHA.Matrix, HHA.DP)
HHA.MatrixDP <- JoinLayers(HHA.MatrixDP)
#Combining idents
HHA.MatrixDP$CombinedAnnotation <- case_when(HHA.MatrixDP$MatrixAnnotationFine %in% unique(HHA.Matrix$MatrixAnnotationFine) ~ HHA.MatrixDP$MatrixAnnotationFine,
                                             .default = "Dermal_Papilla")
Idents(HHA.Matrix) <- "CombinedAnnotation"
table(HHA.MatrixDP@active.ident)
#------Creating CellChat Object----#
MatrixDP.CC <- createCellChat(object = HHA.MatrixDP, group.by = "ident", assay = "RNA")
#Setting CellChat database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
#Subset for ECM receptor interactions
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor", key = "annotation")
#Setting Database for CellChat object
MatrixDP.CC@DB <- CellChatDB.use

#-------Running CellChat---------#
#Subsetting data for signaling genes
MatrixDP.CC <- subsetData(MatrixDP.CC)
#Running CellChat
MatrixDP.CC <- identifyOverExpressedGenes(MatrixDP.CC)
MatrixDP.CC <- identifyOverExpressedInteractions(MatrixDP.CC)
#Computing probabilities for receptor ligand pairs
MatrixDP.CC <- computeCommunProb(MatrixDP.CC, type = "triMean")
#Computing probabilities for pathways
MatrixDP.CC <- computeCommunProbPathway(MatrixDP.CC)
#Aggregating network
MatrixDP.CC <- aggregateNet(MatrixDP.CC)

#------Visualizing Interaction Numbers/Strengths-----#
groupSize <- as.numeric(table(MatrixDP.CC@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(MatrixDP.CC@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(MatrixDP.CC@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#Looking at set of enriched pathways
SigPathways <- MatrixDP.CC@netP$pathways
SigPathways
#Plotting relative gene contributions to pathway
pathways.show = "COLLAGEN"
netAnalysis_contribution( MatrixDP.CC, signaling = pathways.show)
# Hierarchy plot
par(mfrow=c(1,1))
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(MatrixDP.CC, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
netVisual_aggregate(MatrixDP.CC, signaling = pathways.show, layout = "circle")
# Chord diagram
netVisual_aggregate(MatrixDP.CC, signaling = pathways.show, layout = "chord")
# Heatmap
netVisual_heatmap(MatrixDP.CC, signaling = pathways.show, color.heatmap = "Reds")

#-----Computing Network Centrality for Identified Pathways--------#
MatrixDP.CC <- netAnalysis_computeCentrality(MatrixDP.CC, slot.name = "netP")
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(MatrixDP.CC, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(MatrixDP.CC)
# Signaling role analysis on the cell-cell communication networks of interest
netAnalysis_signalingRole_scatter(MatrixDP.CC, signaling = c("COLLAGEN", "LAMININ"))
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#Outgoing
netAnalysis_signalingRole_heatmap(MatrixDP.CC, pattern = "outgoing", signaling = c("COLLAGEN", "LAMININ", "FN1", "THBS", "TENASCIN", "HSPG"))
#Incoming
netAnalysis_signalingRole_heatmap(MatrixDP.CC, pattern = "incoming", signaling = c("LAMININ", "COLLAGEN", "FN1", "THBS"))
#Incoming
netAnalysis_signalingRole_heatmap(MatrixDP.CC, pattern = "incoming", signaling = c("COLLAGEN", 'LAMININ'))
#Incoming
netAnalysis_signalingRole_heatmap(MatrixDP.CC, pattern = "incoming")
#show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(MatrixDP.CC, remove.isolate = FALSE)

#--------Calculating Adhesion Scores---------#
#Extracting Significant Interaction Pairs
MatrixDP.LR <- extractEnrichedLR(MatrixDP.CC, signaling = c("COLLAGEN", "LAMININ", "FN1", "THBS", "TENASCIN", "HSPG"), geneLR.return = T,
                                 enriched.only = TRUE, thresh = 0.05)
unique(MatrixDP.LR$geneLR)
#Getting list of receptors 
AdhesionReceptors <- c("CD44", "SDC1", "SDC4", "ITGA1", "ITGA2", "ITGA3", "ITGA9", "ITGAV", "ITGB1", "ITGB8", "SV2C", "DAG1", "ITGA6", "ITGB4", "ITGA5", "ITGA8", "ITGB6", "CD47", "TNC", "HSPG2")
#Calculating adhesion and proliferation scores with AUCell
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
Matrix.AdhesionScores <- HHA.Matrix %>% JoinLayers() %>% RunAUCell(list(AdhesionScore = AdhesionReceptors,
                                                                        ProliferationScore = CellCycleGenes))
HHA.Matrix <- AddMetaData(HHA.Matrix, Matrix.AdhesionScores)
#Plotting Adhesion Score
FeaturePlot(HHA.Matrix, features = 'AdhesionScore', reduction = "UMAP", pt.size = 1.5) & scale_color_viridis_c(option = "turbo") & NoAxes()
#Plotting Proliferation Score
FeaturePlot(HHA.Matrix, features = 'ProliferationScore', reduction = "UMAP", pt.size = 1.5) & scale_color_viridis_c(option = "turbo") & NoAxes()

#-------Saving Proliferation and AdhesionScores to csv--------#
AdhesionProliferationScores.DF <- Matrix.AdhesionScores %>% mutate(Barcode = rownames(Matrix.AdhesionScores))
write_csv(AdhesionProliferationScores.DF, file = paste0(PlotsDirectory, "Matrix_Adhesion_Proliferation_Scores.csv"))

#-------Plotting Score Violin Plots------#
#Ordering by Early vs Late Differentiation 
HHA.Matrix$DifferentiationOrder <- factor(HHA.Matrix$MatrixAnnotationFine, levels = c("LPC", "Lower COL17", "Upper COL17", "Early_Cortex", "Early_Cuticle", "Early_IRS_I",
                                                                                      "Early_IRS_II", "Medulla", "Middle_Cortex",  "Late_Cortex", "Middle_Cuticle", "Late_Cuticle",  "IRS_Cuticle", "IRS_Henle", "IRS_Huxley"))
#Palette Reordering 
DifferentiationPal <- c("#208A42", "#89C75F", "#3BBCA8", "#9ECAE1", "#E6C2DC", "#D8A767", "#F47D2B", "#0C727C", "#4292C6", "#08306B","#C06CAB" ,"#89288F" , "#F37B7D" ,"#7E1416" ,"#D24B27")
#Plotting Adhesion
Adhesion.ViolinPlot <- VlnPlot(HHA.Matrix, features = "AdhesionScore", group.by = "DifferentiationOrder", alpha = 0, raster = F, cols = DifferentiationPal) & labs(x = "") &
  theme(aspect.ratio = 1) & NoLegend()
Adhesion.ViolinPlot
#Plotting Proliferation
Proliferation.ViolinPlot <- VlnPlot(HHA.Matrix, features = "ProliferationScore", group.by = "DifferentiationOrder", alpha = 0, raster = F, cols = DifferentiationPal) & labs(x = "") &
  theme(aspect.ratio = 1) & NoLegend()
Proliferation.ViolinPlot

#------------Saving Plots-----------#
#Adhesion
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_Adhesion_Scores_by_Cluster_Violins.png"), res = 600, height = 6, width = 8, units = "in") 
Adhesion.ViolinPlot
dev.off()
#Proliferation
ragg::agg_tiff(filename = paste0(PlotsDirectory, "Matrix_Proliferation_Scores_by_Cluster_Violins.png"), res = 600, height = 6, width = 8, units = "in") 
Proliferation.ViolinPlot
dev.off()

#-------Subsetting for COL17 Populations-------#
#Pulling populations from the V4 kit for analysis
V4.COL17 <- HHA.Matrix %>% JoinLayers() %>% subset(MatrixAnnotationFine %in% c("Lower COL17", "Upper COL17") & Platform == "V4")
#Adding Simplified Factored Annotation Column 
V4.COL17$COL17Annotation <- factor(V4.COL17$MatrixAnnotationFine, levels = c("Lower COL17", "Upper COL17"))
#Plotting Annotations per Cluster
VlnPlot(V4.COL17, features = "AdhesionScore", group.by = "MatrixAnnotationOrdered", alpha = 0, raster = F, cols = c("#89C75F", "#3BBCA8", "#0C727C")) & labs(x = "")

#----------------Computing Adhesion Terciles-------------------------------#
V4.AdhesionBinning <- V4.COL17[[]] %>% group_by(MatrixAnnotationOrdered) %>% mutate(AdQ50 = quantile(AdhesionScore, prob = 0.5),
                                                                                    ClusterAdhesionBin = factor(case_when(AdhesionScore <= AdQ50 ~"Low_Adhesion",
                                                                                                                          AdhesionScore > AdQ50 ~"High_Adhesion"),
                                                                                                                levels = c("High_Adhesion", "Low_Adhesion"))) %>%
  select(Barcode, MatrixAnnotationOrdered, AdQ50, ClusterAdhesionBin, ProliferationScore, S.Score, G2M.Score) %>% ungroup()
#Adding Combined Annotation Adhesion Column
V4.AdhesionBinning$AnnotAdhesion <- factor(paste(str_replace(V4.AdhesionBinning$MatrixAnnotationOrdered, pattern = " ", replacement = "_"), V4.AdhesionBinning$ClusterAdhesionBin, sep = "_"),
                                           levels = c("Lower_COL17_High_Adhesion", "Lower_COL17_Low_Adhesion",
                                                      "Upper_COL17_High_Adhesion", "Upper_COL17_Low_Adhesion"))
#Adding scores back to Seurat
V4.COL17$ClusterAdhesionGroup <- V4.AdhesionBinning$ClusterAdhesionBin 


#---------------Plotting Adhesion and Proliferation------------#
AdProlif.BoxPlot <- ggboxplot(V4.AdhesionBinning, x = "MatrixAnnotationOrdered", y = "ProliferationScore", fill = "ClusterAdhesionBin", palette = c("#4A5D73", "#C84A2A"),
                              legend = "none", add = "jitter", add.params = list(alpha = 1, size  = 2, color = "ClusterAdhesionBin", width = 0.2), alpha = 0.6) +
  labs(x = "", y = "Proliferation Score") + theme_classic() + theme(aspect.ratio = 0.75, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  RotatedAxis() + NoLegend()
AdProlif.BoxPlot

#---------------Computing Stats------------#
AdProlif.Comps <- list(c("Lower_COL17_High_Adhesion","Lower_COL17_Low_Adhesion"),
                       c("Upper_COL17_High_Adhesion","Upper_COL17_Low_Adhesion"))
#-----------Plotting with Stats-----------------#
#Calculating Statistics
#Kruskall-Wallis
AdProlif.KW <- kruskal_test(V4.AdhesionBinning, ProliferationScore ~ AnnotAdhesion) %>% add_significance("p")
AdProlif.KW
#Pairwise Wilcoxon
AdProlif.Wilcox <- V4.AdhesionBinning %>% wilcox_test(ProliferationScore ~ AnnotAdhesion, comparisons = AdProlif.Comps) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")
AdProlif.Wilcox
#Writing to csv
write_csv(AdProlif.KW, file = paste0(PlotsDirectory, "HHA_COL17_AdhesionTerciles_Proliferation_Kruskall_Wallis.csv"))
write_csv(AdProlif.Wilcox, file = paste0(PlotsDirectory, "HHA_COL17_AdhesionTerciles_Proliferation_Wilcoxon.csv"))

#-------------Plotting with Stats---------------#
AdProlif.WilcoxPlot <- AdProlif.Wilcox %>% mutate(y.position = max(V4.AdhesionBinning$ProliferationScore) + row_number() * 0.03) 
AdProlif.StatsPlot <- ggboxplot(V4.AdhesionBinning, x = "AnnotAdhesion", y = "ProliferationScore", fill = 'ClusterAdhesionBin', palette = c("#4A5D73", "#C84A2A"), alpha = 0.5, legend = "none", add = "jitter", add.params = list(alpha = 1, size  = 2, color = "ClusterAdhesionBin"), alpha = 0.6) +
  stat_pvalue_manual(AdProlif.WilcoxPlot,label = "p.adj.signif", y.position = "y.position", tip.length = 0.01, hide.ns = F) +
  labs(x = "", y = "Proliferation Score") + theme_classic() + theme(aspect.ratio = 0.75, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  RotatedAxis() + NoLegend()
AdProlif.StatsPlot

#----------------Saving as PNG------------------$
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Proliferation_Adhesion_BoxPlot.png"), res = 600, height = 8, width = 8, units = "in") 
AdProlif.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Proliferation_Adhesion_BoxPlot_Stats.png"), res = 600, height = 8, width = 8, units = "in") 
AdProlif.StatsPlot
dev.off()

#---------------Plotting Proliferation Stats Across Inner Populations----------------#
#Subsetting for V4 Kit, Upper/Lower COl17, Medulla
V4.InnerProlif <- HHA.Matrix[[]] %>% filter(MatrixAnnotationOrdered %in% c("Upper COL17", "Lower COL17", "Medulla") & Platform == "V4") %>% mutate(COL17Annotation = factor(MatrixAnnotationOrdered, levels = c("Lower COL17", "Upper COL17", "Medulla")))
InnerProlif.BoxPlot <- ggboxplot(V4.InnerProlif, x = "COL17Annotation", y = "ProliferationScore", fill = "COL17Annotation", palette = c("#89C75F", "#3BBCA8", "#0C727C"),
                                 legend = "none", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = COL17Annotation), width = 0.05, height = 0, size = 2, alpha = 0.75) +
  labs(x = "", y = "Proliferation Score") + theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  RotatedAxis() + NoLegend()
InnerProlif.BoxPlot

#-----------Running Stats-----------------#
#Calculating Statistics
#Kruskall-Wallis
InnerProlif.KW <- kruskal_test(V4.InnerProlif, ProliferationScore ~ COL17Annotation)
InnerProlif.KW
#Pairwise Wilcoxon
InnerProlif.Wilcox <- V4.InnerProlif %>% pairwise_wilcox_test(ProliferationScore ~ COL17Annotation) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")
InnerProlif.Wilcox
#Writing to csv
write_csv(InnerProlif.KW, file = paste0(PlotsDirectory, "HHA_COL17_Medulla_Proliferation_Kruskall_Wallis.csv"))
write_csv(AdProlif.Wilcox, file = paste0(PlotsDirectory, "HHA_COL17_Medulla_Proliferation_Wilcoxon.csv"))

#-----------Plotting with Stats-----------------#
InnerProlif.WilcoxPlot <- InnerProlif.Wilcox %>% mutate(y.position = max(V4.InnerProlif$ProliferationScore) + row_number() * 0.03) 
InnerProlif.StatsPlot <- ggboxplot(V4.InnerProlif, x = "COL17Annotation", y = "ProliferationScore", fill = "COL17Annotation", palette = c("#89C75F", "#3BBCA8", "#0C727C"),
                                   legend = "none", alpha = 0.5, outlier.shape = NA) + 
  stat_pvalue_manual(InnerProlif.WilcoxPlot,label = "p.adj.signif", y.position = "y.position", tip.length = 0.01, hide.ns = F) +
  geom_jitter(aes(color = COL17Annotation), width = 0.05, height = 0, size = 2, alpha = 0.75) +
  labs(x = "", y = "Proliferation Score") + theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
  RotatedAxis() + NoLegend()
InnerProlif.StatsPlot

#----------------Saving as PNG------------------$
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Medulla_Proliferation_BoxPlot.png"), res = 600, height = 8, width = 8, units = "in") 
InnerProlif.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_Medulla_Proliferation_BoxPlot_Stats.png"), res = 600, height = 8, width = 8, units = "in") 
InnerProlif.StatsPlot
dev.off()

#----------------Computing Adhesion Groups: COL17 Layer-------------------------------#
HHA.COL17 <- HHA.Matrix %>% JoinLayers() %>% subset(MatrixAnnotationOrdered %in% c("Upper COL17", "Lower COL17"))
#Computing Quantiles by Platform and Cluster
Full.AdhesionBinning <- HHA.COL17[[]] %>% group_by(Platform, MatrixAnnotationOrdered) %>% mutate(AdQ50 = quantile(AdhesionScore, prob = 0.5),
                                                                                                 ClusterAdhesionBin = factor(case_when(AdhesionScore <= AdQ50 ~"Low_Adhesion",
                                                                                                                                       AdhesionScore > AdQ50 ~"High_Adhesion"),
                                                                                                                             levels = c("High_Adhesion", "Low_Adhesion"))) %>%
  select(Barcode, MatrixAnnotationOrdered, AdQ50, ClusterAdhesionBin, ProliferationScore, S.Score, G2M.Score) %>% ungroup()
Full.AdhesionBinning
#Adding scores back to Seurat
HHA.COL17$ClusterAdhesionGroup <- Full.AdhesionBinning$ClusterAdhesionBin
#Validating Groupings
ggplot(HHA.COL17[[]], aes(x = MatrixAnnotationOrdered, y = AdhesionScore, fill = ClusterAdhesionGroup)) +
  geom_violin() + facet_wrap(~Platform) + theme_classic()


#-----------------Validation Plots-----------------#
#Checking group Assignments
VlnPlot(HHA.COL17, features = "AdhesionScore", group.by = "MatrixAnnotationOrdered", split.by = "ClusterAdhesionGroup", alpha = 0, raster = F, cols = c("#4A5D73", "#C84A2A")) & labs(x = "")
#Checking group Assignments
VlnPlot(HHA.COL17, features = "ProliferationScore", group.by = "MatrixAnnotationOrdered", split.by = "ClusterAdhesionGroup", alpha = 0, raster = F, cols = c("#4A5D73", "#C84A2A")) & labs(x = "")
#Checking for similar count distributions
VlnPlot(HHA.COL17, features = "nFeature_RNA", group.by = "MatrixAnnotationOrdered", split.by = "ClusterAdhesionGroup", alpha = 0, raster = F, cols = c("#4A5D73", "#C84A2A")) & labs(x = "")

#---------------Generating DEGs between COL17A1 Populations (Including Medulla)-----------------#
#Subsetting for COL17A1 + Medulla
HHA.GL <- HHA.Matrix %>% JoinLayers() %>% subset(MatrixAnnotationFine %in% c("Upper COL17", 'Lower COL17', 'Medulla'))
HHA.GL
#----------Finding Markers-------------#
HHA.GL$COL17Annotation <- factor(HHA.GL$MatrixAnnotationFine, levels = c('Lower COL17', "Upper COL17", 'Medulla'))
Idents(HHA.GL) <- "COL17Annotation"
#Running Marker Test
GL.Markers <- FindAllMarkers(HHA.GL, test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.1, only.pos = T) 
GL.Markers <- GL.Markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5)
GL.Markers
#Writing to csv
write_csv(GL.Markers, file = paste0(PlotsDirectory, "HHA_COL17_Medulla_SigMarkers.csv"))



#---------Saving Workspace-------#
save.image(paste0(SeuratDirectory, "HHA_Matrix_Adhesion_Proliferation_Analysis_1_26_26.RData"))


