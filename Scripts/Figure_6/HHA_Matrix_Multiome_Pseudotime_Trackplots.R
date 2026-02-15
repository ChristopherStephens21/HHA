#------Pre-Processing and Analysis of HHA Multiomics Data-----#
#This script will involve GRN analysis with Pando
#-----Loading in Packages-----#
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86) 
library(JASPAR2020)
library(TFBSTools)
library(universalmotif)
library(chromVAR)
library(Pando)
#For Heatmaps
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(scCustomize)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Bulb_Seurat_Plots/Figure_6/Trackplots/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#--------------Plotting--------------#
#------Palettes-----#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
#From MetBrewer::met.brewer("Signac", 10) "#92c051" "#2b9b81" "#1f6e9c" "#633372" "#9f5691" "#e87b89" "#de597c" "#d8443c" "#fe9b00" "#f4c40f"
SampleIDPal <- c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f")
#Plotting Pseudotime Quartiles
Cortex.Pal <- c("#89C75F", "#3BBCA8", "#9ECAE1", "#4292C6", "#08306B")
IRS.Pal <- c( "#3BBCA8", "#89C75F", "#D8A767", "#F47D2B", "#F37B7D")

#--------Loading Multiome Object-----#
load(paste0(SeuratDirectory, file = "HHA_Matrix_Multiome_Motifs_7_21_25.rds"))
Matrix.Multiome

#--------Loading Full HHA------#
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))
HHA.Matrix

#----------Examining GEP Trends Across Pseudotime: Cortex-------#
#Reading in Pseudotime Results: Cortex
#Pseudotime
Cortex.Pseudotime <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_Cortex_Cuticle_Pseudotime.csv")  %>% 
  select(barcode = "barcode...1", "palantir_pseudotime_CoCu", "palantir_entropy_CoCu")
Cortex.Pseudotime
#Branch Assignments 
Cortex.Branches <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_Cortex_Cuticle_Branches.csv") %>% 
  select(barcode = "barcode...1", "Cortex_Terminal", "Cuticle_Terminal")
Cortex.BranchMetadata <- Cortex.Branches %>% mutate(Barcode = barcode) %>% left_join(HHA.Matrix[[]], by = "Barcode")
Cortex.BranchMetadata
#Selecting Multiome nuclei within the IRS Cuticle Branch
Cortex.BranchPseudotime <- full_join(Cortex.BranchMetadata, Cortex.Pseudotime, by = "barcode") %>% filter(Platform == "Multiome") %>% 
  select(Barcode, Cortex_Terminal, Cuticle_Terminal, palantir_pseudotime_CoCu, palantir_entropy_CoCu)
Cortex.Branch <- Cortex.BranchPseudotime %>% filter(Cortex_Terminal == T) %>% select(Barcode)
#Adding column for whether cells are in Cortex branch
Matrix.Multiome$CortexBranch  <- Matrix.Multiome$Barcode %in% Cortex.Branch$Barcode 
#Adding Pseudotime
Cortex.PseudotimeReordered <- full_join(Matrix.Multiome[[]], Cortex.BranchPseudotime, by = "Barcode") %>% select(Barcode, palantir_pseudotime_CoCu)
Matrix.Multiome$palantir_pseudotime_CoCu <- Cortex.PseudotimeReordered$palantir_pseudotime_CoCu
#Plotting Pseudotime
FeaturePlot(Matrix.Multiome, features = "palantir_pseudotime_CoCu", reduction = "UMAP_WNN") & scale_color_viridis_c(option = "turbo")

#------Calculating Pseudotime Quantiles-----#
Cortex.Subset <- Matrix.Multiome[[]] %>% filter(CortexBranch == T & MatrixAnnotationFine != "Upper COL17") 
Cortex.Terciles <- quantile(Cortex.Subset$palantir_pseudotime_CoCu, c(0.33, 0.66))
Matrix.Multiome$CortexPseudotimeGrouping <- case_when(Matrix.Multiome$CortexBranch == F & Matrix.Multiome$MatrixAnnotationFine == "Lower COL17" ~ "Lower COL17",
                                                      Matrix.Multiome$CortexBranch == F & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" ~ "Other",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine == "Upper COL17" ~ "Upper COL17",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" & Matrix.Multiome$palantir_pseudotime_CoCu <= Cortex.Terciles[1] ~ "Early Cortex",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" & Matrix.Multiome$palantir_pseudotime_CoCu > Cortex.Terciles[1] & Matrix.Multiome$palantir_pseudotime_CoCu <= Cortex.Terciles[2]  ~ "Middle Cortex",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" & Matrix.Multiome$palantir_pseudotime_CoCu > Cortex.Terciles[2] ~ "Late Cortex",) %>% 
  factor(levels = c("Lower COL17", "Upper COL17", "Early Cortex", "Middle Cortex", "Late Cortex", "Other"))
Cortex.Median <- median(Cortex.Subset$palantir_pseudotime_CoCu)
Matrix.Multiome$CortexPseudotimeGrouping <- case_when(Matrix.Multiome$CortexBranch == F & Matrix.Multiome$MatrixAnnotationFine == "Lower COL17" ~ "Lower COL17",
                                                      Matrix.Multiome$CortexBranch == F & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" ~ "Other",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine == "Upper COL17" ~ "Upper COL17",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" & Matrix.Multiome$palantir_pseudotime_CoCu <= Cortex.Median ~ "Early Cortex",
                                                      Matrix.Multiome$CortexBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" & Matrix.Multiome$palantir_pseudotime_CoCu > Cortex.Median ~ "Late Cortex",) %>% 
  factor(levels = c("Lower COL17", "Upper COL17", "Early Cortex", "Late Cortex", "Other"))
DimPlot(Matrix.Multiome, group.by = "CortexPseudotimeGrouping", reduction = "UMAP_WNN") 
#Plotting Bins
Matrix.Multiome[[]] %>% filter(CortexPseudotimeGrouping %in% c("Other", "Lower COL17") == F) %>% 
  ggplot(aes(x = palantir_pseudotime_CoCu, y = CortexPseudotimeGrouping, color = CortexPseudotimeGrouping)) + 
  ggbeeswarm::geom_quasirandom() + theme_bw() & scale_color_manual(values = Cortex.Pal)

#----------Examining GEP Trends Across Pseudotime: IRSCuticle-------#
#Reading in Pseudotime Results: IRSCuticle
#Pseudotime
IRSCuticle.Pseudotime <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_IRS_Pseudotime.csv")  
IRSCuticle.Pseudotime
#Branch Assignments 
IRSCuticle.Branches <- read_csv("/projects/b1217/HHA/Bulb_Recluster_5_22_AnnData/Palantir_IRS_Branches.csv") 
IRSCuticle.BranchMetadata <- IRSCuticle.Branches %>% mutate(Barcode = barcode) %>% left_join(HHA.Matrix[[]], by = "Barcode")
IRSCuticle.BranchMetadata
#Selecting Multiome nuclei within the IRS Cuticle Branch
IRSCuticle.BranchPseudotime <- full_join(IRSCuticle.BranchMetadata, IRSCuticle.Pseudotime, by = "barcode") %>% filter(Platform == "Multiome") %>% 
  select(Barcode, IRS_Cuticle_Terminal, IRS_Henle_Terminal, IRS_Huxley_Terminal, palantir_pseudotime_IRS, palantir_entropy_IRS)
IRSCuticle.Branch <- IRSCuticle.BranchPseudotime %>% filter(IRS_Cuticle_Terminal == T) %>% select(Barcode)
#Adding column for whether cells are in IRSCuticle branch
Matrix.Multiome$IRSCuticleBranch  <- Matrix.Multiome$Barcode %in% IRSCuticle.Branch$Barcode 
#Adding Pseudotime
IRSCuticle.PseudotimeReordered <- full_join(Matrix.Multiome[[]], IRSCuticle.BranchPseudotime, by = "Barcode") %>% select(Barcode, palantir_pseudotime_IRS)
Matrix.Multiome$palantir_pseudotime_IRS <- IRSCuticle.PseudotimeReordered$palantir_pseudotime_IRS
#Plotting Pseudotime
FeaturePlot(Matrix.Multiome, features = "palantir_pseudotime_IRS", reduction = "UMAP_WNN") & scale_color_viridis_c(option = "turbo")

#------Calculating Pseudotime Quantiles-----#
IRSCuticle.Subset <- Matrix.Multiome[[]] %>% filter(IRSCuticleBranch == T & MatrixAnnotationFine != "Lower COL17") 
IRSCuticle.Terciles <- quantile(IRSCuticle.Subset$palantir_pseudotime_IRS, c(0.33, 0.66))
Matrix.Multiome$IRSCuticlePseudotimeGrouping <- case_when(Matrix.Multiome$IRSCuticleBranch == F & Matrix.Multiome$MatrixAnnotationFine == "Upper COL17" ~ "Upper COL17",
                                                          Matrix.Multiome$IRSCuticleBranch == F & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" ~ "Other",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine == "Lower COL17" ~ "Lower COL17",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" & Matrix.Multiome$palantir_pseudotime_IRS <= IRSCuticle.Terciles[1] ~ "Early IRS Cuticle",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" & Matrix.Multiome$palantir_pseudotime_IRS > IRSCuticle.Terciles[1] & Matrix.Multiome$palantir_pseudotime_IRS <= IRSCuticle.Terciles[2]  ~ "Middle IRS Cuticle",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" & Matrix.Multiome$palantir_pseudotime_IRS > IRSCuticle.Terciles[2] ~ "Late IRS Cuticle") %>% 
  factor(levels = c("Upper COL17", "Lower COL17", "Early IRS Cuticle", "Middle IRS Cuticle", "Late IRS Cuticle", "Other"))

IRSCuticle.Median <- median(IRSCuticle.Subset$palantir_pseudotime_IRS)
Matrix.Multiome$IRSCuticlePseudotimeGrouping <- case_when(Matrix.Multiome$IRSCuticleBranch == F & Matrix.Multiome$MatrixAnnotationFine == "Upper COL17" ~ "Upper COL17",
                                                          Matrix.Multiome$IRSCuticleBranch == F & Matrix.Multiome$MatrixAnnotationFine != "Upper COL17" ~ "Other",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine == "Lower COL17" ~ "Lower COL17",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" & Matrix.Multiome$palantir_pseudotime_IRS <= IRSCuticle.Median ~ "Early IRS Cuticle",
                                                          Matrix.Multiome$IRSCuticleBranch == T & Matrix.Multiome$MatrixAnnotationFine != "Lower COL17" & Matrix.Multiome$palantir_pseudotime_IRS > IRSCuticle.Median ~ "Late IRS Cuticle") %>% 
  factor(levels = c("Upper COL17", "Lower COL17", "Early IRS Cuticle", "Late IRS Cuticle", "Other"))
DimPlot(Matrix.Multiome, group.by = "IRSCuticlePseudotimeGrouping", reduction = "UMAP_WNN") 

#Plotting Bins
Matrix.Multiome[[]] %>% filter(IRSCuticlePseudotimeGrouping %in% c("Other", "Upper COL17") == F) %>% 
  ggplot(aes(x = palantir_pseudotime_IRS, y = IRSCuticlePseudotimeGrouping, color = IRSCuticlePseudotimeGrouping)) + 
  ggbeeswarm::geom_quasirandom() + theme_bw() + scale_color_manual(values = IRS.Pal)

#--------Plotting TrackPlots by Pseudotime Bin-----#
#-----Loading in Significant Peak2Gene Links-------#
Matrix.P2G <- read_csv("/projects/b1217/HHA/Matrix_ATAC_SEACells/Matrix_HRG_PG_Links_Stress_CC_Filtered_6_5_25.csv")
#Converting Peak names to match signac format
Matrix.P2G$peak <- str_replace(Matrix.P2G$peak, pattern = ":", replacement = "-")
Matrix.P2G

#-----KRT36------#
DefaultAssay(Matrix.Multiome) <- "peaks"
KRT36.Peaks <- Matrix.P2G %>% filter(gene == "KRT36") %>% pull(peak)
KRT36.PeakGR <- StringToGRanges(regions = KRT36.Peaks, sep = c("-", "-"))
KRT36.Coords <- LookupGeneCoords(Matrix.Multiome, "KRT36")
KRT36.Bounds <- paste(seqnames(KRT36.Coords), start(KRT36.Coords) - 10000, end(KRT36.Coords) + 10000, sep = "-") %>% 
  StringToGRanges(sep = c("-", "-"))
KRT36.Highlight <- subsetByOverlaps(KRT36.PeakGR, KRT36.Bounds)
DefaultAssay(Matrix.Multiome) <- "peaks"
Idents(Matrix.Multiome) <- "CortexPseudotimeGrouping"
#Plotting
KRT36.LinkagePlot <- CoveragePlot(
  object = Matrix.Multiome,
  region = c("KRT36"),
  region.highlight = KRT36.Highlight,
  features = "KRT36",
  idents = c("Lower COL17", "Upper COL17", "Early Cortex", "Late Cortex"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1)  & scale_fill_manual(values = Cortex.Pal)
KRT36.LinkagePlot
#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "KRT36_Cortex_PseudotimeGrouping_Early_Late_LinkagePlot.tiff"), res = 600, height = 8, width = 12, units = "in") 
KRT36.LinkagePlot
dev.off()

#----------Adding IRS Pseudotime to Seurat Object-------#
#Refactoring Annotations for IRS Plotting 
Matrix.Multiome$AnnotationRefactoredIRS <- factor(case_when(Matrix.Multiome$MatrixAnnotationFine %in% c("Early_IRS_I", "Early_IRS_II") ~ "Early_IRS",
                                                            .default = Matrix.Multiome$MatrixAnnotationFine), 
                                                  levels = c("Upper COL17", "Lower COL17", "Early_IRS", "IRS_Henle", "IRS_Huxley", "IRS_Cuticle",
                                                             "LPC", "Medulla", "Early_Cortex", "Middle_Cortex", "Late_Cortex", "Early_Cuticle", "Middle_Cuticle", "Late_Cuticle"))

#-----KRT73------#
KRT73.Peaks <- Matrix.P2G %>% filter(gene == "KRT73") %>% pull(peak)
KRT73.PeakGR <- StringToGRanges(regions = KRT73.Peaks, sep = c("-", "-"))
KRT73.Coords <- LookupGeneCoords(Matrix.Multiome, "KRT73")
KRT73.Bounds <- paste(seqnames(KRT73.Coords), start(KRT73.Coords) - 10000, end(KRT73.Coords) + 10000, sep = "-") %>% 
  StringToGRanges(sep = c("-", "-"))
KRT73.Highlight <- subsetByOverlaps(KRT73.PeakGR, KRT73.Bounds)
DefaultAssay(Matrix.Multiome) <- "peaks"
Idents(Matrix.Multiome) <- "IRSCuticlePseudotimeGrouping"
#Plotting
KRT73.LinkagePlot <- CoveragePlot(
  object = Matrix.Multiome,
  region = c("KRT73"),
  region.highlight = KRT73.Highlight,
  idents = c("Upper COL17", "Lower COL17", "Early IRS Cuticle", "Late IRS Cuticle"),
  features = "KRT73",
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1)  & scale_fill_manual(values = IRS.Pal)
KRT73.LinkagePlot

#-----Saving to tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "KRT73_IRS_PseudotimeGrouping__Early_Late_LinkagePlot.tiff"), res = 600, height = 8, width = 12, units = "in") 
KRT73.LinkagePlot
dev.off()

#--------Plotting GATA3 over Time--------#
HHA.Matrix$MatrixDynamicsOrder <- factor(HHA.Matrix$MatrixAnnotationFine, levels = c("Lower COL17", "Early_IRS_I", "Early_IRS_II",
                                                                                     "IRS_Henle", "IRS_Huxley", "IRS_Cuticle",
                                                                                     "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex",
                                                                                     "LPC", "Medulla", "Early_Cuticle", "Middle_Cuticle", "Late_Cuticle"))
Idents(HHA.Matrix) <- "MatrixDynamicsOrder"
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT36", "EXT1"), stack = T, flip = T)

#-----Plotting Violin Plots Lineage Marker Progressions-----#
ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT36_EXT1.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT36", "EXT1"), stack = T, flip = T)
dev.off()

ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT36_TP63.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT36", "TP63"), stack = T, flip = T)
dev.off()

ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT36.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT36"), stack = T, flip = T)
dev.off()

ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT85.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT85"), stack = T, flip = T)
dev.off()

ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT85_TP63.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT85", "TP63"), stack = T, flip = T)
dev.off()

ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/Dynamics_Vln_COL17A1_KRT85_EXT1.tiff"), res = 600, height = 8, width = 8, units = "in") 
VlnPlot(HHA.Matrix, idents = c("Lower COL17", "Lower COL17", "Early_IRS_I", "Early_IRS_II", "IRS_Cuticle",
                               "Upper COL17", "Early_Cortex", "Middle_Cortex", "Late_Cortex"),
        features = c("COL17A1", "GATA3", "KRT85", "EXT1"), stack = T, flip = T)
dev.off()


#---------Plotting Violin Plots: Upper vs Lower COL17--------#
HHA.COL17 <- HHA.Matrix %>% JoinLayers() %>% subset(MatrixAnnotationFine %in% c("LPC", "Lower COL17", "Upper COL17", "Medulla"))
HHA.COL17$MatrixAnnotationFineOrdered <- factor(HHA.COL17$MatrixAnnotationFine, levels = c("LPC", "Lower COL17", "Upper COL17", "Medulla"))
Idents(HHA.COL17) <- 'MatrixAnnotationFineOrdered'
COL17Pal <- c("#208A42", "#89C75F", "#3BBCA8", "#0C727C")
VlnPlot(HHA.COL17, features = c("COL17A1", "LGR5", "TBX1", "GATA3", "SHH", "FRY", "CPA6", "KRT35", "HOXC13", "FOXQ1", "TCHH", "KRT17"),
        stack = T, flip = T, fill.by = "ident", idents = c("LPC", "Lower COL17", "Upper COL17", "Medulla"), cols = COL17Pal) & theme(aspect.ratio = 0.2) & NoLegend()



