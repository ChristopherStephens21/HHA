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
PlotsDirectory <- "/projects/b1217/HHA/Revisions/Bulge_Revisions_Plots/"
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

#-------------Loading Checkpoint-----------#
load(paste0(SeuratDirectory, "HHA_ORS_Bulge_Quiescence_Proliferation_Analysis_1_27_26.RData"))

#---------Loading Full Object-------#
load(paste0(SeuratDirectory, "HHA_SCRNA_Multiome_Integrated_Annotated_3_31_25.rds"))
HHA.Full

#----------------Plotting Clusters----------------#
#GlobalAnnotation
DimPlot(HHA.Full, reduction = "UMAP", group.by = "FinalAnnotation", label = T, repel = T, raster = F) +
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()
#Broad Annotation
DimPlot(HHA.Full, reduction = "UMAP", group.by = "GeneralAnnotation", label = T, repel = T, raster = F) +
  labs(title = "The Human Scalp Atlas: Global Populations", subtitle = "128,185 Cells", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1) & NoLegend()

#-------Subsetting Out Bulge Cells-----#
HHA.ORS <- subset(HHA.Full, FinalAnnotation %in%  c("Bulge", "ORS_Basal", "ORS_Suprabasal", "Isthmus_Suprabasal"))
HHA.ORS #31,954 cells, 7,806 multiome 
DimPlot(HHA.ORS, reduction = "UMAP", group.by = "FinalAnnotation", label = T, repel = T, raster = F)
table(HHA.ORS$Platform)
table(HHA.ORS$FinalAnnotation, HHA.ORS$Platform)

#-----Recomputing Neighbors and UMAP in original scVI Space-----#
HHA.ORS <- FindNeighbors(HHA.ORS, dims = 1:50, reduction = "scvi")
HHA.ORS <- RunUMAP(HHA.ORS, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
#Plotting
DimPlot(HHA.ORS, reduction = "UMAP", group.by = "FinalAnnotation", cols = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'), pt.size = 1, label = T, repel = T, raster = F)

#------------Calculating Proliferation and Quiescence Scores with AUCell-------------#
#From Rando & Chung (2013) Molecular regulation of stem cell quiescence
Quiescence <- read_csv("/projects/b1217/HHA/Dissociation_Signature/Chung_Rando_2013_Quiescence_Sig.csv")
#Pulling Gene Set and converting to human notation 
Quiescence.Genes <- str_to_upper(Quiescence$Gene)
Quiescence.Genes[Quiescence.Genes %in% rownames(HHA.ORS) ==F] #GLTSCR2 is absence from scrnaseq data, will not be used
Quiescence.Genes 
#Pulling Cell Cycle Genes
CellCycleGenes <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
#Running AUCell
ORS.QuiescenceScores <- HHA.ORS %>% JoinLayers() %>% RunAUCell(list(ProliferationScore = CellCycleGenes,
                                                                    QuiescenceScore = Quiescence.Genes))
HHA.ORS <- AddMetaData(HHA.ORS, ORS.QuiescenceScores)
#Plotting Quiescence Score
FeaturePlot(HHA.ORS, features = 'QuiescenceScore', reduction = "UMAP", pt.size = 1) & scale_color_viridis_c(option = "mako") & NoAxes()
#Plotting Proliferation Score
FeaturePlot(HHA.ORS, features = 'ProliferationScore', reduction = "UMAP", pt.size = 1) & scale_color_viridis_c(option = "plasma") & NoAxes()

#--------Subsetting for V4 Kit & Quantifying--------------#
#Pulling V4 Kit Scores 
V4.QuiescenceScores <- HHA.ORS[[]] %>% filter(Platform == "V4") %>% select(FinalAnnotation, ProliferationScore, QuiescenceScore) %>% mutate(FinalAnnotation = factor(FinalAnnotation, levels = c("Bulge", "Isthmus_Suprabasal", "ORS_Basal", "ORS_Suprabasal")))
V4.QuiescenceScores

#------------Computing Statistics--------------#
#Computing Kruskall-Wallis Test
Quiescence.KW <- kruskal_test(V4.QuiescenceScores, QuiescenceScore ~ FinalAnnotation)
Quiescence.KW
#Writing to csv
write_csv(Quiescence.KW, file = paste0(PlotsDirectory, "HHA_ORS_Populations_Quiescence_Score_Kruskall_Wallis.csv"))
#Computing Pairwise Wilcoxon
Quiescence.Wilcox <- V4.QuiescenceScores %>%
  pairwise_wilcox_test(QuiescenceScore ~ FinalAnnotation) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")
Quiescence.Wilcox 
#Writing to csv
write_csv(Quiescence.Wilcox, file = paste0(PlotsDirectory, "HHA_ORS_Populations_Quiescence_Score_Wilcoxon.csv"))

#------------Plotting BoxPlots: Quiescence--------------#
Quiescence.BoxPlot <- ggboxplot(V4.QuiescenceScores, x = "FinalAnnotation", y = "QuiescenceScore", fill = "FinalAnnotation", legend = "none", palette =c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'), add = "jitter",
                                add.params = list(alpha = 0.025, size = 0.35)) + labs(x = "", y = "Quiescence Score", fill = "Annotation") + coord_flip() + 
  theme(aspect.ratio = 0.66, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis() +
  guides(color = guide_legend(override.aes = list(shape = 16,size  = 4,alpha = 1, linetype = 0))) 
Quiescence.BoxPlot

#----------------Plotting with Statistics--------------#
#Adding Y Positions for Plotting 
Quiescence.WilcoxPlot <- Quiescence.Wilcox %>% mutate(y.position = max(V4.QuiescenceScores$QuiescenceScore) + row_number() * 0.03) 
#Plotting with Stats
Quiescence.StatsPlot <- ggboxplot(V4.QuiescenceScores, x = "FinalAnnotation", y = "QuiescenceScore", fill = "FinalAnnotation", legend = "none", palette =c(c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda')), add = "jitter",
                                add.params = list(alpha = 0.025, size = 0.35)) + labs(x = "", y = "Quiescence Score", color = "Cluster") + # Add pairwise comparisons p-value
  stat_pvalue_manual(Quiescence.WilcoxPlot,label = "p.adj.signif", y.position = "y.position", tip.length = 0.01, hide.ns = F) + theme_classic() + coord_flip() +
  theme(aspect.ratio = 0.66, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis() + NoLegend()
Quiescence.StatsPlot 

#-------------Saving as PNGs------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Quiescence_Score_BoxPlot.png"), res = 600, height = 6, width = 8, units = "in") 
Quiescence.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Quiescence_Score_Stats_BoxPlot.png"), res = 600, height = 6, width = 8, units = "in") 
Quiescence.StatsPlot
dev.off()

#-----------------Quantifying Number of Proliferating Cells--------------#
#Restricting to V4 Data
V4.ORS <- HHA.ORS %>% JoinLayers %>% subset(Platform == "V4")
V4.ORS
#Pulling raw count matrix
V4.counts <- GetAssayData(V4.ORS, assay = "RNA", slot = "counts")
V4.counts
#Counting all cells with Ki67 > 1 as positive 
V4.ORS$KI67Positive <- V4.counts['MKI67',] > 0
#Ki67 positive cells reflect true cycling population 
VlnPlot(V4.ORS, group.by = "FinalAnnotation", split.by = "KI67Positive", features = c("ProliferationScore", "QuiescenceScore", "S.Score", "G2M.Score"), alpha = 0.01, ncol = 2)
#Subsetting for 4 mm Isolations
V4.Metadata <- V4.ORS[[]] %>% filter(Isolation == "4mm")
table(V4.Metadata$SampleID, V4.Metadata$DonorID) # 5 samples from 5 distinct donors 
V4.Metadata
#Computing Proliferation Percentages by Donor 
Prolif.Rates <- V4.Metadata %>% group_by(SampleID, FinalAnnotation) %>% summarize(Percent_MKI67 = mean(KI67Positive) * 100,
                                                                                  N_Ki67 = sum(KI67Positive),
                                                                                  N_Total = n(),
                                                                                  ) %>% ungroup()
Prolif.Rates$FinalAnnotation <- factor(Prolif.Rates$FinalAnnotation, levels = c("Bulge", "Isthmus_Suprabasal", "ORS_Basal", "ORS_Suprabasal"))
Prolif.Rates

#--------Plotting RawStats as BoxPlot------------#
ggplot(Prolif.Rates, aes(x = Percent_MKI67, fill = FinalAnnotation, y = FinalAnnotation)) + geom_boxplot() + scale_fill_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda')) + theme_bw() +
  labs(x = "", y = "", title = "Percent MKI67 Positive") + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "None")

#------------Computing Proliferation Stats: Binomial Mixed Model-----------#
#Annotation as fixed effect, sample as random effect
Ki67.Model <- glmer(cbind(N_Ki67, N_Total - N_Ki67) ~ FinalAnnotation + (1 | SampleID),
  family = binomial, data = Prolif.Rates)
#Computing marginal means
Ki67.EMM <- emmeans(Ki67.Model, ~ FinalAnnotation, type = "response")
Ki67.EMM
#Computing Pairwise comparison p values
Ki67.Wald  <- pairs(Ki67.EMM, adjust = "BH")
Ki67.Wald
#Saving EMMs to csv
write_csv(as.data.frame(Ki67.EMM), file = paste0(PlotsDirectory, "HHA_ORS_Populations_Ki67_LogReg_Estimated_Marginal_Means.csv"))
#Saving contrasts to csv
write_csv(as.data.frame(Ki67.Wald), file = paste0(PlotsDirectory, "HHA_ORS_Populations_Ki67_Pairwise_Contrasts_Wald.csv"))

#------------Formatting Results for Plotting-------------#
#Getting estimated proportions of cycling cells 
Ki67.EMMPlot <- as.data.frame(Ki67.EMM) %>%
  mutate(FinalAnnotation = factor(FinalAnnotation, levels = FinalAnnotation),
    pct = 100 * prob, pct_lcl = 100 * asymp.LCL, pct_ucl = 100 * asymp.UCL)
#Formatting Statistics for Plotting 
Ki67.WaldPlot <- as.data.frame(Ki67.Wald) %>% tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") %>% 
  mutate(p.adj = p.value, p.adj.signif = case_when(p.adj < 0.001 ~ "***",
                                                   p.adj < 0.05  ~ "*",
                                                   TRUE ~ "ns"),
         group1 = factor(group1, levels = levels(Prolif.Rates$FinalAnnotation)),
         group2 = factor(group2, levels = levels(Prolif.Rates$FinalAnnotation))) %>%
  add_y_position(data = Prolif.Rates, formula = Percent_MKI67 ~ FinalAnnotation, step.increase = 0.05)

#-------------Plotting Model------------#
#Plotting without Stats
Ki67.ModelPlot <- ggplot() + 
  #Observed percentages
  geom_point(data = Prolif.Rates, aes(x = FinalAnnotation, y = Percent_MKI67, color = FinalAnnotation, size = N_Total), alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) + #Observed percentages
  #Model with Confidence interval
  geom_errorbar(data = Ki67.EMMPlot, aes(x = FinalAnnotation, ymin = pct_lcl, ymax = pct_ucl), width = 0.15) +
  #EMMs
  geom_point(data = Ki67.EMMPlot, aes(x = FinalAnnotation, y = pct), size = 2) + labs(x = "", y = "% MKI67 Positivity", size = "Total Cells", title = "% MKI67 Positivity", color = "Population") + coord_flip() +
  theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.5) + scale_color_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'))
Ki67.ModelPlot
#Plotting with Stats
Ki67.ModelStatsPlot <- ggplot() + 
  #Observed percentages
  geom_point(data = Prolif.Rates, aes(x = FinalAnnotation, y = Percent_MKI67, color = FinalAnnotation, size = N_Total), alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) + #Observed percentages
  #Model with Confidence interval
  geom_errorbar(data = Ki67.EMMPlot, aes(x = FinalAnnotation, ymin = pct_lcl, ymax = pct_ucl), width = 0.15) +
  #EMMs
  geom_point(data = Ki67.EMMPlot, aes(x = FinalAnnotation, y = pct), size = 2) +
  #Significance
  ggpubr::stat_pvalue_manual(Ki67.WaldPlot, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) + coord_flip() + labs(x = "", y = "% MKI67 Positivity", size = "Total Cells", title = "% MKI67 Positivity", color = "Population") +
  theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.5) + scale_color_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'))
Ki67.ModelStatsPlot

#----------Saving as PNG------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_KI67_LogReg_Raincloud.png"), res = 600, height = 8, width = 12, units = "in") 
Ki67.ModelPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_KI67_LogReg_Raincloud_Stats.png"), res = 600, height = 6, width = 12, units = "in") 
Ki67.ModelStatsPlot
dev.off()

#-----------Plotting Raw Score BoxPlot---------#
#Without Stats
Prolif.BoxPlot <- ggpubr::ggboxplot(Prolif.Rates, x = "FinalAnnotation", y = "Percent_MKI67", fill = "FinalAnnotation",
                                    palette = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'), legend = "none", add = "jitter", add.params = list(alpha = 0.6, size = 3)) +
  coord_flip() + theme_classic() + labs(x = "", y = "% MKI67 Positivity", title = "% MKI67 Positivity") + theme(plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio = 0.66) + NoLegend()
Prolif.BoxPlot
#With Stats
Prolif.StatsPlot <- ggpubr::ggboxplot(Prolif.Rates, x = "FinalAnnotation", y = "Percent_MKI67", fill = "FinalAnnotation",
  palette = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda'), legend = "none", add = "jitter", add.params = list(alpha = 0.6, size = 2)) +
  ggpubr::stat_pvalue_manual(Ki67.WaldPlot, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) +
  coord_flip() + theme_classic() + labs(x = "", y = "% MKI67 Positivity", title = "% MKI67 Positivity") + theme(plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio = 0.66) + NoLegend()
Prolif.StatsPlot

#----------Saving as PNG------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_KI67_BoxPlot.png"), res = 600, height = 6, width = 8, units = "in") 
Prolif.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_KI67_BoxPlot_Stats.png"), res = 600, height = 6, width = 8, units = "in") 
Prolif.StatsPlot
dev.off()

#--------Loading Full HHA Matrix Object------#
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))
HHA.Matrix

#--------------Calculating Matrical Proliferation Rates------------#
#Extracting Matrical Cells from V4 Kit
V4.Matrix <- HHA.Matrix %>% JoinLayers() %>% subset(Platform == "V4")
V4.Matrix
#Plotting 
DimPlot(V4.Matrix, group.by = "MatrixAnnotationFine", reduction = "UMAP")
table(V4.Matrix$DonorID, V4.Matrix$SampleID)

#-----------------Quantifying Number of Proliferating Cells--------------#
#Adding GeneralAnnotation Column
V4.Matrix$GeneralAnnotation <- "HF_Matrix"
#Pulling raw count matrix
Matrix.counts <- GetAssayData(V4.Matrix, assay = "RNA", slot = "counts")
Matrix.counts
#Counting all cells with Ki67 > 1 as positive 
V4.Matrix$KI67Positive <- Matrix.counts['MKI67',] > 0
#Computing Proliferation Percentages by Donor/Sample
Matrix.ProlifRates <- V4.Matrix[[]] %>% group_by(SampleID, GeneralAnnotation) %>% summarize(Percent_MKI67 = mean(KI67Positive) * 100,
                                                                                  N_Ki67 = sum(KI67Positive),
                                                                                  N_Total = n()) %>% ungroup() %>% select(FinalAnnotation = GeneralAnnotation, everything())
Matrix.ProlifRates
#Combining with ORS Measurements 
Combined.ProlifRates <- Prolif.Rates %>% bind_rows(Matrix.ProlifRates) %>% mutate(FinalAnnotation = factor(FinalAnnotation, levels = c("Bulge", "Isthmus_Suprabasal", "ORS_Basal", "ORS_Suprabasal", "HF_Matrix")))
Combined.ProlifRates
#--------Plotting RawStats as BoxPlot------------#
ggplot(Combined.ProlifRates, aes(x = Percent_MKI67, fill = FinalAnnotation, y = FinalAnnotation)) + geom_boxplot() + scale_fill_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', '#185a88')) + theme_bw() +
  labs(x = "", y = "", title = "Percent MKI67 Positive") + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "None")

#------------Computing Proliferation Stats: Binomial Mixed Model-----------#
#Annotation as fixed effect, sample as random effect
Combined.Model <- glmer(cbind(N_Ki67, N_Total - N_Ki67) ~ FinalAnnotation + (1 | SampleID),
                    family = binomial, data = Combined.ProlifRates)
#Computing marginal means
Combined.EMM <- emmeans(Combined.Model, ~ FinalAnnotation, type = "response")
Combined.EMM
#Computing Pairwise comparison p values
Combined.Wald  <- pairs(Combined.EMM, adjust = "BH")
Combined.Wald
#Saving EMMs to csv
write_csv(as.data.frame(Ki67.EMM), file = paste0(PlotsDirectory, "HHA_ORS_Populations__Matrix_Ki67_LogReg_Estimated_Marginal_Means.csv"))
#Saving contrasts to csv
write_csv(as.data.frame(Ki67.Wald), file = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_Ki67_Pairwise_Contrasts_Wald.csv"))

#------------Formatting Results for Plotting-------------#
#Getting estimated proportions of cycling cells 
Combined.EMMPlot <- as.data.frame(Combined.EMM) %>%
  mutate(FinalAnnotation = factor(FinalAnnotation, levels = FinalAnnotation),
         pct = 100 * prob, pct_lcl = 100 * asymp.LCL, pct_ucl = 100 * asymp.UCL)
#Formatting Statistics for Plotting 
Combined.WaldPlot <- as.data.frame(Combined.Wald) %>% tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") %>% 
  mutate(p.adj = p.value, p.adj.signif = case_when(p.adj < 0.001 ~ "***",
                                                   p.adj < 0.05  ~ "*",
                                                   TRUE ~ "ns"),
         group1 = factor(group1, levels = levels(Combined.ProlifRates$FinalAnnotation)),
         group2 = factor(group2, levels = levels(Combined.ProlifRates$FinalAnnotation))) %>%
  add_y_position(data = Combined.ProlifRates, formula = Percent_MKI67 ~ FinalAnnotation, step.increase = 0.05)

#-------------Plotting Model------------#
#Plotting without Stats
Combined.ModelPlot <- ggplot() + 
  #Observed percentages
  geom_point(data = Combined.ProlifRates, aes(x = FinalAnnotation, y = Percent_MKI67, color = FinalAnnotation, size = N_Total), alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) + #Observed percentages
  #Model with Confidence interval
  geom_errorbar(data = Combined.EMMPlot, aes(x = FinalAnnotation, ymin = pct_lcl, ymax = pct_ucl), width = 0.15) +
  #EMMs
  geom_point(data = Combined.EMMPlot, aes(x = FinalAnnotation, y = pct), size = 2) + labs(x = "", y = "% MKI67 Positivity", size = "Total Cells", title = "% MKI67 Positivity", color = "Population") + coord_flip() +
  theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.5) + scale_color_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', '#185a88'))
Combined.ModelPlot
#Plotting with Stats
Combined.ModelStatsPlot <- ggplot() + 
  #Observed percentages
  geom_point(data = Combined.ProlifRates, aes(x = FinalAnnotation, y = Percent_MKI67, color = FinalAnnotation, size = N_Total), alpha = 0.85, position = position_jitter(width = 0.12, height = 0)) + #Observed percentages
  #Model with Confidence interval
  geom_errorbar(data = Combined.EMMPlot, aes(x = FinalAnnotation, ymin = pct_lcl, ymax = pct_ucl), width = 0.15) +
  #EMMs
  geom_point(data = Combined.EMMPlot, aes(x = FinalAnnotation, y = pct), size = 2) +
  #Significance
  ggpubr::stat_pvalue_manual(Combined.WaldPlot, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) + coord_flip() + labs(x = "", y = "% MKI67 Positivity", size = "Total Cells", title = "% MKI67 Positivity", color = "Population") +
  theme_classic() + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.5) + scale_color_manual(values = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', '#185a88'))
Combined.ModelStatsPlot

#----------Saving as PNG------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_KI67_LogReg_Raincloud.png"), res = 600, height = 6, width = 12, units = "in") 
Combined.ModelPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_KI67_Matrix_LogReg_Raincloud_Stats.png"), res = 600, height = 6, width = 12, units = "in") 
Combined.ModelStatsPlot
dev.off()

#-----------Plotting Raw Score BoxPlot---------#
#Without Stats
CombinedProlif.BoxPlot <- ggpubr::ggboxplot(Combined.ProlifRates, x = "FinalAnnotation", y = "Percent_MKI67", fill = "FinalAnnotation",
                                    palette = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', '#185a88'), legend = "none", add = "jitter", add.params = list(alpha = 0.6, size = 3)) +
  coord_flip() + theme_classic() + labs(x = "", y = "% MKI67 Positivity", title = "% MKI67 Positivity") + theme(plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio = 0.66) & NoLegend()
CombinedProlif.BoxPlot
#With Stats
CombinedProlif.StatsPlot <- ggpubr::ggboxplot(Combined.ProlifRates,, x = "FinalAnnotation", y = "Percent_MKI67", fill = "FinalAnnotation",
                                      palette = c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', '#185a88'), legend = "none", add = "jitter", add.params = list(alpha = 0.6, size = 2)) +
  ggpubr::stat_pvalue_manual(Combined.WaldPlot, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) +
  coord_flip() + theme_classic() + labs(x = "", y = "% MKI67 Positivity", title = "% MKI67 Positivity") + theme(plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio = 0.66) & NoLegend()
CombinedProlif.StatsPlot

#----------Saving as PNG------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_KI67_BoxPlot.png"), res = 600, height = 6, width = 8, units = "in") 
CombinedProlif.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_KI67_BoxPlot_Stats.png"), res = 600, height = 6, width = 8, units = "in") 
CombinedProlif.StatsPlot
dev.off()



#------------------------Analyzing Quiescence Gene Expression with Matrical Cells--------------------#
#Running AUCell
Matrix.QuiescenceScores <-V4.Matrix %>% JoinLayers() %>%  RunAUCell(list(ProliferationScore = CellCycleGenes,
                                                                    QuiescenceScore = Quiescence.Genes))
HHA.Matrix <- AddMetaData(V4.Matrix, Matrix.QuiescenceScores)
#Pulling Quiescence Scores into DataFrame
Matrix.Quiescence <- HHA.Matrix[[]] %>% select(FinalAnnotation = GeneralAnnotation, QuiescenceScore, ProliferationScore)
Combined.Quiescence <- bind_rows(V4.QuiescenceScores, Matrix.Quiescence)  %>% mutate(FinalAnnotation = factor(FinalAnnotation, levels = c("Bulge", "Isthmus_Suprabasal", "ORS_Basal", "ORS_Suprabasal", "HF_Matrix")))
Combined.Quiescence


#------------Computing Statistics--------------#
#Computing Kruskall-Wallis Test
MatrixQuiescence.KW <- kruskal_test(Combined.Quiescence, QuiescenceScore ~ FinalAnnotation)
MatrixQuiescence.KW
#Writing to csv
write_csv(MatrixQuiescence.KW, file = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_Quiescence_Score_Kruskall_Wallis.csv"))
#Computing Pairwise Wilcoxon
MatrixQuiescence.Wilcox <- Combined.Quiescence %>%
  pairwise_wilcox_test(QuiescenceScore ~ FinalAnnotation) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")
MatrixQuiescence.Wilcox 
#Writing to csv
write_csv(MatrixQuiescence.Wilcox, file = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_Quiescence_Score_Wilcoxon.csv"))

#------------Plotting BoxPlots: Quiescence--------------#
MatrixQuiescence.BoxPlot <- ggboxplot(Combined.Quiescence, x = "FinalAnnotation", y = "QuiescenceScore", fill = "FinalAnnotation", legend = "none", palette =c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', "#185a88"), add = "jitter",
                                add.params = list(alpha = 0.025, size = 0.35)) + labs(x = "", y = "Quiescence Score", fill = "Annotation") + coord_flip() + theme_classic() +
  theme(aspect.ratio = 0.66, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis() +
  guides(color = guide_legend(override.aes = list(shape = 16,size  = 4,alpha = 1, linetype = 0))) + NoLegend()
MatrixQuiescence.BoxPlot

#----------------Plotting with Statistics--------------#
#Adding Y Positions for Plotting 
MatrixQuiescence.WilcoxPlot <- MatrixQuiescence.Wilcox %>% mutate(y.position = max(Combined.Quiescence$QuiescenceScore) + row_number() * 0.03) 
#Plotting with Stats
MatrixQuiescence.StatsPlot <- ggboxplot(Combined.Quiescence, x = "FinalAnnotation", y = "QuiescenceScore", fill = "FinalAnnotation", legend = "none", palette =c(c('#ad8fc3', '#ddd1e7', '#1294a1', '#76cbda', "#185a88")), add = "jitter",
                                  add.params = list(alpha = 0.025, size = 0.35)) + labs(x = "", y = "Quiescence Score", color = "Cluster") + # Add pairwise comparisons p-value
  stat_pvalue_manual(MatrixQuiescence.WilcoxPlot,label = "p.adj.signif", y.position = "y.position", tip.length = 0.01, hide.ns = F) + theme_classic() + coord_flip() +
  theme(aspect.ratio = 0.66, plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) + RotatedAxis() + NoLegend()
MatrixQuiescence.StatsPlot 


#-----------Saving Plots---------------#
#Without Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_Quiescence_Score_BoxPlot.png"), res = 600, height = 6, width = 8, units = "in") 
MatrixQuiescence.BoxPlot
dev.off()
#With Stats
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HHA_ORS_Populations_Matrix_Quiescence_Score_BoxPlot_Stats.png"), res = 600, height = 6, width = 8, units = "in") 
MatrixQuiescence.StatsPlot
dev.off()

#------------Saving Workspace Image--------------#
save.image(paste0(SeuratDirectory, "HHA_ORS_Bulge_Quiescence_Proliferation_Analysis_1_27_26.RData"))



















