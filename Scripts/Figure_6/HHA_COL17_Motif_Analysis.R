#------Pre-Processing and Analysis of HHA Multiomics Data-----#
#This script contains our workflow for our Matrix Motif Analysis 
#-----Loading Packages-----#
library(Seurat)
library(scCustomize)
library(scDblFinder)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(presto)
library(AUCell)
library(ggpubr)
library(tidyverse)

#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Revisions/Matrix_Motif_Analysis/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#------------------Loading Checkpoint--------------#
load(paste0(SeuratDirectory, "HHA_COL17_Motif_Analysis_2_4_26.RData"))

#---------Loading Objects------#
#ATAC Only
load(paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_Signac_6_16_25.rds"))
Matrix.ATAC
#Multiome RNA/ATAC
load(paste0(SeuratDirectory, "HHA_Multiome_ATAC_Matrix_Integrated_Annotated_6_16_25.rds"))
Matrix.Multiome
#Full RNA
load(paste0(SeuratDirectory, "HHA_Matrix_Integrated_Annotated_5_30_25.rds"))

#-------Setting Default Assay for Multiome Data to RNA-------#
DefaultAssay(Matrix.Multiome) <- "RNA"
#Joining Layers
Matrix.Multiome <- JoinLayers(Matrix.Multiome)

#------Ordering Clusters-------#
MatrixPal <- c("#89C75F", "#3BBCA8",  "#208A42", "#0C727C", "#9ECAE1", "#4292C6", "#08306B", "#E6C2DC", "#C06CAB",  "#89288F", "#D8A767", "#F47D2B", "#F37B7D",  "#7E1416", "#D24B27")
Matrix.Multiome$MatrixAnnotationOrdered <- factor(Matrix.Multiome$MatrixAnnotationFine, levels = c("Lower COL17", "Upper COL17", "LPC", "Medulla", 
                                                                                                   "Early_Cortex", "Middle_Cortex", "Late_Cortex", 
                                                                                                   "Early_Cuticle", "Middle_Cuticle", "Late_Cuticle",
                                                                                                   "Early_IRS_I", "Early_IRS_II", "IRS_Henle", "IRS_Huxley", "IRS_Cuticle"))
#Plotting
DimPlot(Matrix.Multiome, reduction = "UMAP_WNN", label = T, group.by = "MatrixAnnotationOrdered", raster = F, cols = MatrixPal, pt.size = 2)


#------------Pulling JASPAR Motif PFM-------------#
#Pulling PFM from JASPAR 2020
Jaspar.PFM <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
Jaspar.PFM

Motif2TFs
#Collapsing entries with duplicated rows from distinct Uniprot(?) accession IDs
Motif.MetadataList <- map(Jaspar.PFM, ~TFBSTools::tags(.x) %>% as.data.frame() %>% select(!acc) %>% distinct())
Motif.MetadataList
identical(names(Motif.MetadataList), names(Jaspar.PFM))
Motif.Metadata <- bind_rows(Motif.MetadataList)

#-----------Retrieving Family Information----------#
Motif.Metadata <- bind_rows(lapply(Jaspar.PFM, function(m) {
  tg <- TFBSTools::tags(m)
  tibble(MotifID = TFBSTools::ID(m),
    TF  = TFBSTools::name(m), 
    #Grouping motifs for TF complexes with multiple members, or replacing with NA
    TFFamily = if (!is.null(tg$family)) paste(tg$family, collapse = "::") else NA_character_)}))
Motif.Metadata
table(Motif.Metadata$TFFamily) %>% sort()

#---------------Building Peak:Motif Matrix---------------#
DefaultAssay(Matrix.Multiome) <- "peaks"
#Scoring motifs across all peaks
Matrix.Multiome <- AddMotifs(Matrix.Multiome, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = Jaspar.PFM)
#Extracting peaks
Matrix.Peaks <-DietSeurat(Matrix.Multiome, assays = "peaks")
#Resetting Granges names from cell type to peak name
new.granges <- granges(Matrix.Peaks) #This avoids an error in RSE construction during RunChromVAR
names(new.granges) <- paste(seqnames(granges(Matrix.Peaks)), start(granges(Matrix.Peaks)), end(granges(Matrix.Peaks)), sep = "-")
Matrix.Peaks@assays$peaks@ranges <- new.granges

#-----------Running ChromVAR-----------#
#Running ChromVAR
Matrix.Peaks <- RunChromVAR(Matrix.Peaks, genome = BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(Matrix.Peaks) <- 'chromvar'
#Adding to Main Object
Matrix.Multiome[["chromvar"]] <- Matrix.Peaks[["chromvar"]]

#-------Plotting Selected Motifs------#
DefaultAssay(Matrix.Multiome) <- 'chromvar'
#GATA3
FeaturePlot(Matrix.Multiome, features = "MA0037.3", reduction = "UMAP", pt.size = 3, max.cutoff = 8) & 
  scale_color_gradientn(colors = c('#3361A5','#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'))
#HOXC13
FeaturePlot(Matrix.Multiome, features = "MA0907.1", reduction = "UMAP", pt.size = 1, max.cutoff = 5) & 
  scale_color_gradientn(colors = c('#3361A5','#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'))
#NFIB
FeaturePlot(Matrix.Multiome, features = "MA1643.1", reduction = "UMAP", pt.size = 1, max.cutoff = 4) & 
  scale_color_gradientn(colors = c('#3361A5','#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D'))

#------------Saving Seurat Object with chromVAR Scores------------#
#ChromVAR
Matrix.ChromVAR <- Matrix.Peaks
save(Matrix.ChromVAR, file = paste0(SeuratDirectory, file = "HHA_Matrix_ChromVAR_2_4_26.rds"))
#Full Multiome
save(Matrix.Multiome, file = paste0(SeuratDirectory, file = "HHA_Matrix_Multiome_Motifs_2_4_26.rds"))

#-----------------Identifying Differentially Enriched Motifs: ChromVAR--------------------#
Idents(Matrix.Multiome) <- 'MatrixAnnotationFine'
DefaultAssay(Matrix.Multiome) <- 'chromvar'
#Computing Differential Motif Deviations
COL17.DiffCV <- FindMarkers(object = Matrix.Multiome, ident.1 = 'Upper COL17', ident.2 = 'Lower COL17',mean.fxn = rowMeans, fc.name = "avg_diff", logfc.threshold = 0, min.pct = 0)
#Adding TF names
COL17.DiffCV$Motif <- rownames(COL17.DiffCV)
COL17.DiffCV$TF <- plyr::mapvalues(COL17.DiffCV$Motif, from = Motif.Metadata$MotifID, to = Motif.Metadata$TF)
COL17.DiffCV$TFFamily <- plyr::mapvalues(COL17.DiffCV$Motif, from = Motif.Metadata$MotifID, to = Motif.Metadata$TFFamily)
#Adding TF + Motif
COL17.DiffCV <- COL17.DiffCV %>% mutate(TFMotif = paste0(TF, " (", Motif, ")"))
COL17.DiffCV <- COL17.DiffCV %>% mutate(MotifTF = paste0(Motif, " (", TF, ")"))
#Adding Significance Flag
COL17.DiffCV$Significance <- factor(case_when(COL17.DiffCV$p_val_adj < 0.05 & COL17.DiffCV$avg_diff > 0.5 ~ "Upper COL17-Enriched",
                                              COL17.DiffCV$p_val_adj < 0.05 & COL17.DiffCV$avg_diff < -0.5 ~ "Lower COL17-Enriched",
                                              .default = "Non-Significant"), levels = c("Upper COL17-Enriched", "Lower COL17-Enriched", "Non-Significant"))
COL17.DiffCV

#----------Looking at Top hits------------#
#Pulling Lower COL17 motif hits
CV.LowSig <- COL17.DiffCV %>% filter(Significance == "Lower COL17-Enriched") 
CV.LowSig 
#Counting Motifs by Family
LowSig.Enriched <- CV.UpSig %>% group_by(TFFamily) %>% summarize(N_Enriched = n()) %>% arrange(desc(N_Enriched))
LowSig.Enriched
#Plotting Hits by TF Family
ggplot(LowSig.Enriched, aes(x = factor(TFFamily, levels = LowSig.Enriched$TFFamily), y = N_Enriched)) +
  geom_col(fill = "firebrick") + theme_classic() + RotatedAxis() + labs(x = 'TF Family', y = 'N Motif Hits', title = "Total Differentially Enriched Motif Hits: Lower COL17") + theme(aspect.ratio = 1)
#Pulling Upper COL17 motif hits
CV.UpSig <- COL17.DiffCV %>% filter(Significance == "Upper COL17-Enriched") 
CV.UpSig
#Counting Motifs by Family
UpSig.Enriched <-CV.UpSig %>% group_by(TFFamily) %>% summarize(N_Enriched = n()) %>% arrange(desc(N_Enriched))
UpSig.Enriched
#Plotting Hits by TF Family
ggplot(UpSig.Enriched, aes(x = factor(TFFamily, levels = UpSig.Enriched$TFFamily), y = N_Enriched)) +
  geom_col(fill = "cornflowerblue") + theme_classic() + RotatedAxis() + labs(x = 'TF Family', y = 'N Motif Hits', title = "Total Differentially Enriched Motif Hits: Upper COL17") + theme(aspect.ratio = 1)

#-----------Selecting GATA and HOXC13 Motifs------------------#
GATA.Motifs <- COL17.DiffCV[grep("^GATA", COL17.DiffCV$TF),] %>% filter(Significance == "Lower COL17-Enriched")
NFI.Motifs <- COL17.DiffCV[grep("^NFI", COL17.DiffCV$TF),]  %>% filter(Significance == "Lower COL17-Enriched")
HOXC.Motifs <- COL17.DiffCV[grep("^HOX", COL17.DiffCV$TF),] %>% filter(Significance == "Upper COL17-Enriched")
TopMotifs <- COL17.DiffCV %>% filter(MotifTF %in% c(GATA.Motifs$MotifTF, NFI.Motifs$MotifTF, HOXC.Motifs$MotifTF))
TopMotifs 
TopMotifs <- rbind(head(CV.UpSig, 10), head(CV.LowSig, 10))
TopMotifs 

#--------------Plotting VolcanoPlot----------------#
#Labeled 
CV.VolcanoPlot <- ggplot(COL17.DiffCV, aes(x = avg_diff, y = -log10(p_val), color = Significance)) + geom_point(size = 3) + 
  geom_hline(yintercept = -log10(0.05 / nrow(COL17.DiffCV)), linetype = "dashed") + geom_vline(xintercept = 0.5, linetype = "dashed") + geom_vline(xintercept = -0.5, linetype = "dashed") +
  ggrepel::geom_text_repel(data = TopMotifs, aes(x = avg_diff, y = -log10(p_val), label = MotifTF, family = "Arial"), color = "black", force = 1, max.overlaps = Inf) +
  theme_classic(base_family = "Arial") + scale_color_manual(values = c("cornflowerblue", "firebrick", "gray")) + labs(x = "Mean Score Difference", y = "-Log10 P Value") + theme(aspect.ratio = 0.75)
CV.VolcanoPlot
#Unlabeled
CV.VolcanoPlotUnlabeled <- ggplot(COL17.DiffCV, aes(x = avg_diff, y = -log10(p_val), color = Significance)) + geom_point(size = 3) + 
  geom_hline(yintercept = -log10(0.05 / nrow(COL17.DiffCV)), linetype = "dashed") + geom_vline(xintercept = 0.5, linetype = "dashed") + geom_vline(xintercept = -0.5, linetype = "dashed") +
  theme_classic(base_family = "Arial") + scale_color_manual(values = c("cornflowerblue", "firebrick", "gray")) + labs(x = "Mean Score Difference", y = "-Log10 P Value") + theme(aspect.ratio = 0.75)
CV.VolcanoPlotUnlabeled

#----------Saving as PNG--------------#
#Labeled
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_ChromVAR_Volcano_Labeled.png"), res = 300, height = 10, width = 12, units = "in") 
CV.VolcanoPlot
dev.off()
#Unlabeled
ragg::agg_tiff(filename = paste0(PlotsDirectory, "COL17_ChromVAR_Volcano_Unlabeled.tiff"), res = 300, height = 10, width = 12, units = "in") 
CV.VolcanoPlotUnlabeled
dev.off()

cairo_pdf(filename = paste0(PlotsDirectory, "COL17_ChromVAR_Volcano_Labeled.pdf"),width = 14, height = 12)
print(CV.VolcanoPlot)
dev.off()

cairo_pdf(filename = paste0(PlotsDirectory, "COL17_ChromVAR_Volcano_Unlabeled.pdf"),width = 14, height = 12)
print(CV.VolcanoPlotUnlabeled)
dev.off()

#------------Plotting ViolinPlot of GATA3/HOXC13 Motifs--------#
#ChromVAR
DefaultAssay(Matrix.Multiome) <- "chromvar"
CV.VlnPlot <- VlnPlot(Matrix.Multiome, idents = c("Lower COL17", "Upper COL17"), features = c("MA0037.3", "MA0907.1"), cols = c("#89C75F", "#3BBCA8"), stack = T, flip = T, fill.by = "ident") & theme(aspect.ratio = 1) & NoLegend()
CV.VlnPlot
#RNA
DefaultAssay(Matrix.Multiome) <- "RNA"
CVRNA.VlnPlot <- VlnPlot(Matrix.Multiome, idents = c("Lower COL17", "Upper COL17"), features = c("GATA3", "HOXC13"), cols = c("#89C75F", "#3BBCA8"), stack = T, flip = T, fill.by = "ident") & theme(aspect.ratio = 1) & NoLegend()
CVRNA.VlnPlot
#----------Saving as PNG--------------#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_HOXC13_ChromVAR_Violin.png"), res = 300, height = 6, width = 6, units = "in") 
CV.VlnPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_HOXC13_RNA_Violin.png"), res = 300, height = 6, width = 6, units = "in") 
CVRNA.VlnPlot
dev.off()

#-----------Plotting Associated Motif Logos---------#
#GATA3
GATA3.LogoPlot <- MotifPlot(object = Matrix.Multiome, motifs = "MA0037.3", assay = 'peaks')
GATA3.LogoPlot
#HOXC13
HOXC13.LogoPlot <- MotifPlot(object = Matrix.Multiome, motifs = "MA0907.1", assay = 'peaks')
HOXC13.LogoPlot

#----------Saving as PNG--------------#
#GATA3
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_MA0037.3_MotifLogo.png"), res = 300, height = 2, width = 2, units = "in") 
GATA3.LogoPlot
dev.off()
#HOXC13
ragg::agg_tiff(filename = paste0(PlotsDirectory, "HOXC13_MA0907.1_MotifLogo.png"), res = 300, height = 2, width = 2, units = "in") 
HOXC13.LogoPlot
dev.off()

#-----------Writing ChromVAR Data to csv-----------#
write_csv(COL17.DiffCV, file = paste0(PlotsDirectory, "HHA_COL17_Upper_Lower_Differential_ChromVAR_Motif_Deviations.csv"))

DotPlot(Matrix.Multiome, group.by = "MatrixAnnotationOrdered", features = c("GATA2", "GATA3", "GATA4", "GATA6", "HOXC13", "HOXB13", "HOXA13", "HOXD13", "CDX1", "CDX2", "CDX4")) & 
  scale_color_viridis_c(option = "plasma")

DefaultAssay(Matrix.Multiome) <- 'RNA'
VlnPlot(Matrix.Multiome, group.by = "MatrixAnnotationOrdered", features = c("GATA2", "GATA3", "GATA4", "GATA6", "HOXC13", "HOXB13", "HOXA13", "CDX1", "CDX2"), stack = T, flip = T)

#-----------Saving Workspace--------#
save.image(paste0(SeuratDirectory, "HHA_COL17_Motif_Analysis_2_4_26.RData"))










