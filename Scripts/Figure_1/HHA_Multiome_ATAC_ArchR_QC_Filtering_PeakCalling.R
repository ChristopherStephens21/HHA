#-----Loading in Packages-----#
library(readr)
library(dplyr)
library(ggplot2)
library(ArchR)
library(purrr)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

#-----Directories------#
SeuratDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Multiome/Seurat_Data/"
PlotsDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Multiome/Seurat_Data/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
AnnDataDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Scanpy_Data/"
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"
#Setting Working Directory to Store ArchR Outputs
ArchRDirectory <- "/projects/b1217/HHA/Multiome_ArchR/"
setwd(ArchRDirectory)

#-----Setting up Parallel Computing-----#
addArchRThreads(threads = 62) 

#-----Setting up HG38 Genome-----#
addArchRGenome("hg38")

#-----Loading Workspace Image-----#
load(paste0(WorkspaceDirectory, "HHA_Multiome_ATAC_Level_Processing_5_24_25.Rdata"))

#------Palettes-----#
GrayMagma <-c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF") #ArchR
GrayRocket <-c("grey", "#F4875EFF", "#CB1B4FFF", "#611F53FF", "#03051AFF")#ArchR
GrayMako <- c("grey", "#49C1ADFF", "#357BA2FF", "#3E356BFF", "#0B0405FF")
GrayViridis <- c("grey", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
GrayFireworks <- c("grey", "#2488F0", "#7F3F98", "#E22929", "#FCB31A")
#From MetBrewer::met.brewer("Signac", 10) "#92c051" "#2b9b81" "#1f6e9c" "#633372" "#9f5691" "#e87b89" "#de597c" "#d8443c" "#fe9b00" "#f4c40f"
SampleIDPal <- c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f")
Metadata.Pal <- scale_fill_manual(values = c("#92c051", "#1f6e9c", "#633372", "#d8443c", "#f4c40f"))

#-----Paths to Frag Files-----#
FragPaths <- c(EL_C5 = "/projects/b1042/YiLab/Edward/Atlas/Multiomics/EL_C5_CellRangerMultiOut/outs/atac_fragments.tsv.gz",
               EL_C6 = "/projects/b1042/YiLab/Edward/Atlas/Multiomics/EL_C6_CellRangerMultiOut/outs/atac_fragments.tsv.gz",
               EL_E6 = "/projects/b1042/YiLab/Edward/Atlas/Multiomics/EL_E6_CellRangerMultiOut/outs/atac_fragments.tsv.gz",
               EL_E7 = "/projects/b1042/YiLab/Edward/Atlas/Multiomics/EL_E7_CellRangerMultiOut/outs//atac_fragments.tsv.gz",
               EL_E8 = "/projects/b1042/YiLab/Edward/Atlas/Multiomics/EL_E8_CellRangerMultiOut/outs//atac_fragments.tsv.gz")

#------Reading in Barcodes for Cells Passing RNA-level QC-----#
#Reading in csv
FilteredBarcodes <- read_csv(paste0(MetadataDirectory, "HHA_Multiome_PreDF_FilteredBarcodes_11_14_24.csv"))
FilteredBarcodes
#Converting to list of vectors by dataset
FilteredBarcodesList <- FilteredBarcodes %>% filter() %>% group_by(SampleID) %>% group_split() %>% map(~.x %>% pull(OriginalBarcode))
names(FilteredBarcodesList) <- names(FragPaths)
FilteredBarcodesList

#-----Creating Arrow files for each Multiome Sample-----#
#HHA.ArrowFiles <- createArrowFiles(inputFiles = FragPaths, validBarcodes = FilteredBarcodesList,
#                              sampleNames = names(FragPaths), minTSS = 4, minFrags = 1000, maxFrags = 1e5, addTileMat = TRUE,
#                              addGeneScoreMat = TRUE, verbose = T, force = T)
#-----Inferring Doublets from ATAC-Level Data-----#
#HHA.DoubletScores <- addDoubletScores(input = HHA.ArrowFiles, k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#                                      knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
#                                      LSIMethod = 1)
HHA.ArrowFiles <- c("/projects/b1217/HHA/Multiome_ArchR/EL_C5.arrow", "/projects/b1217/HHA/Multiome_ArchR/EL_C6.arrow", "/projects/b1217/HHA/Multiome_ArchR/EL_E6.arrow",
                    "/projects/b1217/HHA/Multiome_ArchR/EL_E7.arrow", "/projects/b1217/HHA/Multiome_ArchR/EL_E8.arrow")

#-----Creating Multiome ArchR Project for ATAC-Level QC and Doublet Filtering-----#
HHA.archr <- ArchRProject(ArrowFiles = HHA.ArrowFiles, outputDirectory = "Multiome_QC_Filtering_5_24_25", copyArrows = TRUE)
HHA.archr #37,915 cells passing RNA level QC

#------Getting ATAC Level Statistics------#
PreQCMetadata.Multiome <- as.data.frame(getCellColData(HHA.archr)) %>% 
  mutate(Log10_Frags = getCellColData(HHA.archr, select = c("log10(nFrags)"))[,1])
PreQCMetadata.Multiome 

#------Plotting Initial ATAC-level QC Statistics-----#
#Violin Plots
plotPDF(
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", pal = SampleIDPal),
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", pal = SampleIDPal),
  #Log10 Fragment Counts
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", pal = SampleIDPal),
  #Blacklist Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin", pal = SampleIDPal),
  #Nucleosome Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin", pal = SampleIDPal),
  #----Ridge Plots-----#
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges", pal = SampleIDPal),
  #Log10 Fragment Counts
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges", pal = SampleIDPal),
  #Blacklist Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "BlacklistRatio", plotAs = "ridges", pal = SampleIDPal),
  #Nucleosome Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges", pal = SampleIDPal),
  #-----TSS Enrichment vs Log10 Fragments-----#
  ggPoint(x = PreQCMetadata.Multiome$Log10_Frags, y = PreQCMetadata.Multiome$TSSEnrichment, 
          colorDensity = TRUE, continuousSet = "horizonExtra", xlabel = "Log10 Unique Fragments",
          ylabel = "TSS Enrichment") + geom_hline(yintercept = 4, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed") + ggtitle("Multiome Data: Log10 Fragment Count vs TSS Enrichment"),
  name = "HHA Multiome Pre-Filtering QC-Sample-Statistics Violins.pdf", ArchRProj = HHA.archr, addDOC = FALSE, width = 4, height = 4)

#------Examining Quantiles for selected statistics------#
quantile(HHA.archr$TSSEnrichment) #These numbers look reasonable
quantile(HHA.archr$nFrags) #These numbers look reasonable
quantile(HHA.archr$BlacklistRatio) #These numbers look reasonable
quantile(HHA.archr$NucleosomeRatio) #These numbers look reasonable
#Computing proportion of cells retained following filtering from RNA Metrics alone. 
PreQCCellProps.Multiome <- PreQCMetadata.Multiome %>% dplyr::group_by(Sample) %>% 
  dplyr::summarize(PropRetained = mean(BlacklistRatio < 0.02 & NucleosomeRatio < 2),
                   CellsRetained = PropRetained * n())
PreQCCellProps.Multiome
PreQCCellProps.Multiome %>% pull(CellsRetained) %>% sum() #37,874 cells pass RNA and ATAC level filtering

#-----Subsetting Project for High Quality Cells and Filtering Doublets-----#
#---Reading in Final QC Barcodes-----#
PostDFFilteringMetadata <- read_csv(paste0(MetadataDirectory, "HHA_Multiome_PostDoubletFilteringMetadata_11_14_24.csv"))
PostDFFilteringMetadata #this contains doublet annotation from scRNAseq

#-----Subsetting Project for High Quality Cells-----#
#Adding in column with barcode name that matches across modaltieis
HHA.archr$UnifiedBarcode <- str_replace_all(HHA.archr$cellNames, pattern = "#", replacement = "_")
HHA.archr <- HHA.archr[which(HHA.archr$UnifiedBarcode %in% PostDFFilteringMetadata$UnifiedBarcode), ]
#Selecting high quality cells and those not annotated as doublets 
Multiome.PassedFilter <- HHA.archr$cellNames[which(HHA.archr$BlacklistRatio < 0.02 & HHA.archr$NucleosomeRatio < 2 & HHA.archr$UnifiedBarcode %in% PostDFFilteringMetadata$UnifiedBarcode)]
Multiome.PassedFilter #33,199 nuclei retained following final filtering
#Subsetting for selected cells 
HHA.archr <- subsetArchRProject(ArchRProj = HHA.archr, cells = Multiome.PassedFilter, outputDirectory = "Multiome_QC_Doublets_Filtered", dropCells = TRUE,
                                threads = getArchRThreads(), force = T)
HHA.archr #33,199 nuclei retained following filtering

#-----Saving ArchR Project------#
saveArchRProject(ArchRProj = HHA.archr, outputDirectory = "Multiome_QC_Doublets_Filtered", load = T)

#------Reading in AMULET Output-----#
Multiome.AmuletResults <- read_csv(paste0(MetadataDirectory, "HHA_Multiome_Amulet_Results_11_29_24.csv"))

#-----Calculating Doublet Formation Rates-----#
Multiome.AmuletResults %>% dplyr::group_by(Sample) %>% 
  dplyr::summarize(N_Singlet = sum(q.value > 0.05), N_Doublet = sum(q.value < 0.1), Doublet_Rate = N_Doublet/n())
Multiome.AmuletResults %>% mutate(DoubletStatus = factor(case_when(q.value < 0.05 ~ "Doublet", .default = "Singlet"), levels = c("Singlet", "Doublet"))) %>% 
  dplyr::group_by(Sample, DoubletStatus) %>% dplyr::summarize(Median_Sites = median(nAbove2))

#-----Adding in AMULET Results to ArchR Project-----#
HHA.archr$UnifiedBarcode <- str_replace_all(HHA.archr$cellNames, pattern = "#", replacement = "_")
#Reordering to match ArchR Object
Amulet.Joined <- left_join(as.data.frame(getCellColData(HHA.archr)), Multiome.AmuletResults, by = "UnifiedBarcode")
sum(Amulet.Joined$UnifiedBarcode == HHA.archr$UnifiedBarcode) == length(HHA.archr$cellNames)
#Adding Results to ArchR Object
HHA.archr$Amulet_nAbove2 <- Multiome.AmuletResults$nAbove2
HHA.archr$Amulet_total_nAbove2 <- Multiome.AmuletResults$total.nAbove2
HHA.archr$Amulet_p_value <-Multiome.AmuletResults$p.value
HHA.archr$Amulet_q_value <-Multiome.AmuletResults$q.value

#------Loading in RNA-Level Data and Adding to ArchR Project------#
Full.Metadata <- read_csv(paste0("/projects/b1217/HHA/Multiome_Seurat/", "HHA_Multiome_Metadata_Final.csv"))
Multiome.Metadata <- Full.Metadata %>% filter(Platform == "Multiome")
MultiomeMetadata.Combined <- left_join(as.data.frame(getCellColData(HHA.archr)), Multiome.Metadata, by = "UnifiedBarcode")
#Reordering Barcodes
MultiomeMetadata.Combined$UnifiedBarcode == HHA.archr$UnifiedBarcode
#Adding in Relevant Metadata Columns 
HHA.archr$nCount_RNA <- MultiomeMetadata.Combined$nCount_RNA
HHA.archr$nFeature_RNA <- MultiomeMetadata.Combined$nFeature_RNA
HHA.archr$percent.mt <- MultiomeMetadata.Combined$percent.mt
HHA.archr$percent.ribo <- MultiomeMetadata.Combined$percent.ribo
HHA.archr$DonorID <- MultiomeMetadata.Combined$DonorID
HHA.archr$Isolation <- MultiomeMetadata.Combined$Isolation
HHA.archr$Sex <- MultiomeMetadata.Combined$Sex
HHA.archr$Age <- MultiomeMetadata.Combined$Age
HHA.archr$DoubletStatus <- MultiomeMetadata.Combined$DoubletStatus
HHA.archr$scDblFinder.score <- MultiomeMetadata.Combined$scDblFinder.score
HHA.archr$Phase <- MultiomeMetadata.Combined$Phase
HHA.archr$MultiomeAnnotationRNA <- MultiomeMetadata.Combined$MultiomeAnnotationRNA 
HHA.archr$BroadAnnotation <- MultiomeMetadata.Combined$BroadAnnotation
HHA.archr$GlobalAnnotation <- MultiomeMetadata.Combined$GlobalAnnotation
HHA.archr$FinalAnnotation <- MultiomeMetadata.Combined$FinalAnnotation
HHA.archr$GeneralAnnotation <- MultiomeMetadata.Combined$GeneralAnnotation

#------Saving ArchR project-----#
HHA.archr <- saveArchRProject(ArchRProj = HHA.archr, outputDirectory = "Multiome_ATAC_Level_Analysis_5_24_25", load = T)

#-----Getting PostQC Metadata-----#
PostQCMetadata.Multiome <- as.data.frame(getCellColData(HHA.archr)) %>% 
  mutate(Log10_Frags = getCellColData(HHA.archr, select = c("log10(nFrags)"))[,1], 
         UnifiedBarcode = str_replace_all(HHA.archr$cellNames, pattern = "#", replacement = "_"),
         OriginalBarcode = str_remove_all(UnifiedBarcode, pattern = "EL_[A-Z][0-9]_"))
PostQCMetadata.Multiome

#-----Saving Multiome ATAC Metadata to csv-----#
write_csv(PostQCMetadata.Multiome, file = paste0(MetadataDirectory, "HHA_Multiome_ATAC_PostfilteringMetadata_11_14_24.csv"))

#-----Plotting post-Filtering QC Statistics-----#
plotPDF(
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", pal = SampleIDPal),
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", pal = SampleIDPal),
  #Log10 Fragment Counts
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", pal = SampleIDPal),
  #Blacklist Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin", pal = SampleIDPal),
  #Nucleosome Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.4,
             colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin", pal = SampleIDPal),
  #----Ridge Plots-----#
  #TSS Enrichment
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges", pal = SampleIDPal),
  #Log10 Fragment Counts
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges", pal = SampleIDPal),
  #Blacklist Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "BlacklistRatio", plotAs = "ridges", pal = SampleIDPal),
  #Nucleosome Ratio
  plotGroups(ArchRProj = HHA.archr, groupBy = "Sample", alpha = 0.6,
             colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges", pal = SampleIDPal),
  #-----TSS Enrichment vs Log10 Fragments-----#
  ggPoint(x = PostQCMetadata.Multiome$Log10_Frags, y = PostQCMetadata.Multiome$TSSEnrichment, 
          colorDensity = TRUE, continuousSet = "horizonExtra", xlabel = "Log10 Unique Fragments",
          ylabel = "TSS Enrichment") + geom_hline(yintercept = 4, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed") + ggtitle("Multiome Data: Log10 Fragment Count vs TSS Enrichment"),
  name = "HHA Multiome Post-Filtering QC-Sample-Statistics.pdf", ArchRProj = HHA.archr, addDOC = FALSE, width = 4, height = 4)

#-----Fragment Size Distribution-----#
MultiomeFragmentSizesPlot <- plotFragmentSizes(ArchRProj = HHA.archr, pal = SampleIDPal)
MultiomeFragmentSizesPlot + ggtitle("Multiome Data: Fragment Size Distribution")
#Saving to pdf
plotPDF(MultiomeFragmentSizesPlot + ggtitle("Multiome Data: Fragment Size Distribution"), name = "HHA_Multiome_Fragment_Size_Distribution.pdf",
        ArchRProj = HHA.archr, addDOC = FALSE, width = 4, height = 4)

#-----TSS Enrichment Distribution-----#
MultiomeTSSEnrichmentPlot <- plotTSSEnrichment(ArchRProj = HHA.archr, pal = SampleIDPal)
MultiomeTSSEnrichmentPlot + ggtitle("Multiome Data: TSS Enrichment")
#Saving to pdf
plotPDF(MultiomeTSSEnrichmentPlot + ggtitle("Multiome Data: TSS Enrichment"), name = "HHA_Multiome_TSS_Enrichment.pdf",
        ArchRProj = HHA.archr, addDOC = FALSE, width = 4, height = 4)

#-----Computing an Initial Low Dimensional Embedding with Tile Matrix-----#
#----Performing DR with Iterative LSI----#
HHA.archr <- addIterativeLSI(ArchRProj = HHA.archr, useMatrix = "TileMatrix", name = "IterativeLSI", 
                             iterations = 2, clusterParams = list(resolution = c(0.1, 0.2), n.start = 10), 
                             varFeatures = 25000, dimsToUse = 1:30, verbose = T, force = T)
#Adding clusters
HHA.archr <- addClusters(input = HHA.archr, reducedDims = "IterativeLSI",
                         method = "Seurat", name = "UnintegratedClusters", resolution = 0.2, force = T)
#Running UMAP
HHA.archr <- addUMAP(ArchRProj = HHA.archr, reducedDims = "IterativeLSI", name = "UMAP_Unintegrated", 
                     nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)

#------Plotting Clustering Results-----#
#Plotting Clusters vs Samples: Clusters are relatively well mixed by Sample
ggAlignPlots(plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "UnintegratedClusters", embedding = "UMAP_Unintegrated", size = 0.1),
             plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Unintegrated", size = 0.1), type = "h")
#Plotting Clusters vs Annotations
#Global Annotation
ggAlignPlots(plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "GlobalAnnotation", embedding = "UMAP_Unintegrated", size = 0.1),
             plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Unintegrated", size = 0.1), type = "h")
#Broad Annotation
ggAlignPlots(plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "BroadAnnotation", embedding = "UMAP_Unintegrated", size = 0.1),
             plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "UnintegratedClusters", embedding = "UMAP_Unintegrated", size = 0.1), type = "h")
#Multiome Annotation
ggAlignPlots(plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "MultiomeAnnotationRNA", embedding = "UMAP_Unintegrated", size = 0.1),
             plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "UnintegratedClusters", embedding = "UMAP_Unintegrated", size = 0.1), type = "h")
#Checking for Doublet Clusters
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "Amulet_q_value", embedding = "UMAP_Unintegrated", size = 0.1, plotAs = "points")
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP_Unintegrated", size = 0.1, plotAs = "points")

#------Performing Integration with Harmony-----#
#integration
HHA.archr <- addHarmony(ArchRProj = HHA.archr, reducedDims = "IterativeLSI",
                        name = "Harmony", groupBy = "Sample", force = T)
#Adding clusters
HHA.archr <- addClusters(input = HHA.archr, reducedDims = "Harmony",
                         method = "Seurat", name = "HarmonyClusters_0.2", resolution = 0.2, force = T)
HHA.archr <- addClusters(input = HHA.archr, reducedDims = "Harmony",
                         method = "Seurat", name = "HarmonyClusters_0.3", resolution = 0.3, force = T)
HHA.archr <- addClusters(input = HHA.archr, reducedDims = "Harmony",
                         method = "Seurat", name = "HarmonyClusters_0.4", resolution = 0.4, force = T)
#Running UMAP
HHA.archr <- addUMAP(ArchRProj = HHA.archr, reducedDims = "Harmony", name = "UMAP_Harmony", 
                     nNeighbors = 30, minDist = 0.5, metric = "cosine", force = T)

#------Plotting Clustering Results-----#
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "HarmonyClusters_0.4", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic()
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic()
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "DonorID", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic()
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "GeneralAnnotation", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic()
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "FinalAnnotation", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic()

#-----Grouping Subclusters for Peak Calling-----#
HHA.archr$PeakCallingAnnotation <- case_when(
                                             HHA.archr$FinalAnnotation %in% c("Cortex") ~ "Cortex",
                                             HHA.archr$FinalAnnotation %in% c("Cuticle") ~ "Cuticle",
                                             HHA.archr$FinalAnnotation %in% c("IRS") ~ "IRS",
                                             HHA.archr$FinalAnnotation %in% c("Proximal ORS") ~ "Proximal_ORS",
                                             HHA.archr$FinalAnnotation %in% c("Proximal CL") ~ "Proximal_CL",
                                             HHA.archr$FinalAnnotation %in% c("ORS_Basal") ~ "ORS_Basal",
                                             HHA.archr$FinalAnnotation %in% c("ORS_Suprabasal") ~ "ORS_Suprabasal",
                                             HHA.archr$FinalAnnotation %in% c("ORS_Suprabasal") ~ "ORS_Suprabasal",
                                             HHA.archr$FinalAnnotation %in% c("Inf_Basal") ~ "Inf_Basal",
                                             HHA.archr$FinalAnnotation %in% c("Inf_Suprabasal") ~ "Inf_Suprabasal",
                                             HHA.archr$FinalAnnotation %in% c("IFE_Basal") ~ "IFE_Basal",
                                             HHA.archr$FinalAnnotation %in% c("IFE_Suprabasal") ~ "IFE_Suprabasal",
                                             HHA.archr$FinalAnnotation %in% c("Isthmus_Suprabasal") ~ "Isthmus_Suprabasal",
                                             HHA.archr$FinalAnnotation %in% c("SG") ~ "SG",
                                             HHA.archr$FinalAnnotation %in% c("Bulge") ~ "Bulge",
                                             HHA.archr$FinalAnnotation %in% c("Bulge") ~ "Bulge",
                                             HHA.archr$FinalAnnotation %in% c("Fibroblasts I", "Fibroblasts II", "Fibroblasts III", "Fibroblasts IV", "DP") ~ "Fibroblasts",
                                             HHA.archr$FinalAnnotation %in% c("SM I", "SM II", "APM") ~ "Muscle",
                                             HHA.archr$FinalAnnotation %in% c("Myeloid I", "Myeloid II", "Myeloid III", "LCs") ~ "Myeloid",
                                             HHA.archr$FinalAnnotation %in% c("NKT", "T_Helper", "T_Reg") ~ "T_Cells",
                                             HHA.archr$FinalAnnotation %in% c("B", "Mast") ~ "B_Cells",
                                             HHA.archr$FinalAnnotation %in% c("Endothelial", "Pericytes", "Lymphatics") ~ "Endo_LV",
                                             HHA.archr$FinalAnnotation %in% c("Melanocytes I", "Melanocytes II", "Schwann") ~ "Schwann",
                                             HHA.archr$FinalAnnotation %in% c("Eccrine_Epithelial",  "Eccrine_Myoepithelial", "Eccrine_Acrosyringium", "Merkel") ~ "Eccrine_Merkel",
                                             HHA.archr$FinalAnnotation %in% c("Eccrine_Acrosyringium") ~ "Eccrine_Acrosyringium")
table(HHA.archr$PeakCallingAnnotation)
plotEmbedding(ArchRProj = HHA.archr, colorBy = "cellColData", name = "PeakCallingAnnotation", embedding = "UMAP_Harmony", size = 0.1, rastr = F,
              legendSize = 5, baseSize = 10) + theme_classic() 

#-----Calling peaks by broad scRNAseq cluster and adding peak matrix------#
#Pseudobulking by supergroup
HHA.archr <- addGroupCoverages(ArchRProj = HHA.archr, groupBy = "PeakCallingAnnotation", minCells = 80, maxCells = 800, force = T)
#Calling peaks with MACS2
HHA.archr <- addReproduciblePeakSet(ArchRProj = HHA.archr, peakMethod = "Macs2", 
                                    groupBy = "PeakCallingAnnotation", pathToMacs2 = "/projects/b1217/Chris/condaenvs/SCENICPlusEnv/bin/macs2", verbose =T, force = T)
#Adding peak matrix for downstream clustering
HHA.archr <- addPeakMatrix(HHA.archr, verbose = T)

#------Saving ArchR project-----#
HHA.archr <- saveArchRProject(ArchRProj = HHA.archr, outputDirectory = "Multiome_ATAC_Level_Analysis_5_24_25", load = T)

#-----Saving Workspace Image-----#
save.image(paste0(WorkspaceDirectory, "HHA_Multiome_ATAC_Level_Processing_5_24_25.Rdata"))

#-------Printing session info------#
sessionInfo()



