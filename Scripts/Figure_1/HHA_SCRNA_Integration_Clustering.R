#-------Clustering and Integration of HHA scRNA-seq Data-----#
#-----Loading in Packages----#
library(Seurat) #scRNA-seq
library(EnsDb.Hsapiens.v86) #Genome for ClusterProfiler
library(clusterProfiler) #Gene set overrepresentation analysis 
library(tidyverse) #Data handling
library(scCustomize) #Reading in CellBender output

#-----Wrapper Functions-----# 
#Plots in viridis color schemes
PlotGenesViridis <- function(SeuratObj, reduction = "UMAP", genes, viridis_option = "viridis", ...) {
  FeaturePlot(SeuratObj, reduction = reduction, features = genes, ...) & scale_color_viridis_c(option = viridis_option) & 
    labs(x = "UMAP 1", y = "UMAP 2")
}

#-----Directories-----#
SeuratDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/SeuratData/"
PlotsDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Plots/"
AnnDataDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Scanpy_Data/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/"

#-----Loading in Seurat Object-----#
load(file = paste0(SeuratDirectory, "HHA_PreIntegration_10_10_24.RData"))

#-----Performing Log Normalization and Cell Cycle Regression-----#
#Performing Log Normalization 
HHA <- NormalizeData(HHA, normalization.method = "LogNormalize", scale.factor = 10000)
#Splitting our Seurat object into a list of our original datasets
HHA <- SplitObject(HHA, split.by = 'orig.ident')
#Performing cell cycle scoring on each dataset individually. 
HHA <- lapply(X = HHA, FUN = function(x){ #Using 2019 updated list of cell cycle genes from Seurat
  x <- CellCycleScoring(x, s.features= cc.genes.updated.2019$s.genes, g2m.features= cc.genes.updated.2019$g2m.genes, set.ident=TRUE)
})
#Merging our datasets back together
HHA <- merge(x = HHA[[1]], y = HHA[2:length(HHA)])
head(HHA[[]])

#-----Unintegrated Clustering/DR-----#
#Selecting 2000 variable features for dimensionality reduction 
HHA <- FindVariableFeatures(HHA, selection.method = "vst", nfeatures = 2000)
HHA.TopFeaturePlot <- VariableFeaturePlot(HHA) %>% LabelPoints(points = head(VariableFeatures(HHA), 10), repel = TRUE) + ggtitle("Highly Variable Genes")
HHA.TopFeaturePlot
HHA <- ScaleData(HHA)
HHA <- RunPCA(HHA)
DimHeatmap(HHA, dims = 1:15, cells = 1000, balanced = T)
#Performing dimensionality reduction and clustering without integration
HHA <-FindNeighbors(HHA, dims = 1:50, reduction = "pca")
HHA <- RunUMAP(HHA, dims = 1:50, reduction = "pca", reduction.name = "UMAP_Unintegrated")
HHA <- FindClusters(HHA, resolution = 1, cluster.name = "unintegrated_clusters")

#--------Plotting Unintegrated Results-----#
#Plotting unintegrated UMAP
UnintegratedUMAP <- DimPlot(HHA, reduction = "UMAP_Unintegrated", raster = F) + 
  labs(title = "Unintegrated Datasets: Louvain Clusters", subtitle = "50 PCs", x = "UMAP 1", y = "UMAP 2") + 
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15), plot.title = element_text(size = 20, hjust = 0.5)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAP
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "UMAPs/UnintegratedUMAP_10_10_24.tiff"), res = 300, height = 12, width = 12, units = "in") 
UnintegratedUMAP
dev.off()
#Plotting unintegrated UMAP
UnintegratedUMAPBySample <-DimPlot(HHA, reduction = "UMAP_Unintegrated", group.by = "orig.ident", raster = F, cols = MetBrewer::met.brewer("Signac", 11), shuffle = T) + 
  labs(title = "Unintegrated Datasets: Sample ID", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Sample") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAPBySample
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "UMAPs/UnintegratedUMAPBySample_10_10_24.tiff"), res = 300, height = 12, width = 12, units = "in") 
UnintegratedUMAPBySample
dev.off()
UnintegratedUMAPByPlatform <- DimPlot(HHA, reduction = "UMAP_Unintegrated", group.by = "Platform", raster = F, shuffle = T, cols = c("steelblue3", "firebrick")) + 
  labs(title = "Unintegrated Datasets: Sequencing Platform", subtitle = "50 PCs", x = "Unintegrated UMAP 1", y = "Unintegrated UMAP 2", color = "Platform") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15), 
        axis.line= element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends = "last"), linewidth = 0.75),
        axis.title = element_text(face = "bold", size = 15), aspect.ratio = 1)
UnintegratedUMAPByPlatform
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "UMAPs/UnintegratedUMAPByPlatform_10_10_24.tiff"), res = 300, height = 12, width = 12, units = "in") 
UnintegratedUMAPByPlatform
dev.off()
#Plotting QC Metrics
UnintegratedQCMetrics <- PlotGenesViridis(HHA, reduction = "UMAP_Unintegrated", genes = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), viridis_option = "inferno") & NoAxes()
UnintegratedQCMetrics
#Saving as tiff
ragg::agg_tiff(filename = paste0(PlotsDirectory, "UMAPs/UnintegratedQCMetrics_10_10_24.tiff"), res = 300, height = 12, width = 12, units = "in") 
UnintegratedQCMetrics
dev.off()

#-------Performing Scanpy Conversion for SCVI Integration-----#
AnnDataDirectory <- "/projects/p31989/HHA_Spatial_SC_Atlas/Scanpy_Data/"
#Subsetting Atlas for Variable Genes 
Diet.HHA <- HHA[VariableFeatures(HHA)]

#----Writing Metadata----#
HHA.metadata <- Diet.HHA[[]] %>% select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, SampleID, DonorID, Platform, Isolation, Sex, Age, Source, Phase)
HHA.metadata$barcode <- rownames(HHA.metadata)
write.csv(HHA.metadata, file = paste0(AnnDataDirectory, "HHA_Metadata_PreIntegration_10_10_24.csv"), quote = F, row.names = F)
#----Writing Count Matrix and Gene Names----#
#Joining layers and extracting count matrix
DefaultAssay(Diet.HHA) <- "RNA"
DietHHA.counts <- Diet.HHA %>% JoinLayers() %>% SeuratObject::LayerData(assay = 'RNA', layer = 'counts')
#Writing to mtx file 
Matrix::writeMM(DietHHA.counts , file = paste0(AnnDataDirectory, "HHA_RawCounts_10_10_24.mtx"))
#Writing gene names
write.table(data.frame('gene'=rownames(DietHHA.counts)), file = paste0(AnnDataDirectory, "HHA_GeneNames_10_10_24.csv"),
            quote = F, row.names = F, col.names = T)

#-----Saving Results------#
save(HHA, file = paste0(SeuratDirectory, "HHA_PreSCVI_10_10_24.RData"))
rm(DietHHA, HHA.metadata, DietHHA.counts)
load(paste0(SeuratDirectory, "HHA_PreSCVI_10_10_24.RData"))

#------Loading in SCVI Integration Output------#
HHA.scvi <- read_csv(paste0(AnnDataDirectory, "HHA_SCVI_LatentRep_10_10_24.csv")) %>% column_to_rownames("Barcode")
head(HHA.scvi)
# put it back in our original Seurat object
HHA[["scvi"]] <- CreateDimReducObject(embeddings = as.matrix(HHA.scvi), key = "scvi_")
rm(HHA.scvi)

#-----Performing Final DR and Clustering-----#
HHA <-FindNeighbors(HHA, dims = 1:50, reduction = "scvi")
HHA <- RunUMAP(HHA, dims = 1:50, reduction = "scvi", reduction.name = "UMAP")
HHA<- FindClusters(HHA, resolution = 2.0, cluster.name = "scvi_clusters_2.0")
HHA<- FindClusters(HHA, resolution = 1.5, cluster.name = "scvi_clusters_1.5")
HHA<- FindClusters(HHA, resolution = 1.2, cluster.name = "scvi_clusters_1.2")
HHA<- FindClusters(HHA, resolution = 1, cluster.name = "scvi_clusters_1.0")
HHA<- FindClusters(HHA, resolution = 0.8, cluster.name = "scvi_clusters_0.8")
HHA<- FindClusters(HHA, resolution = 0.75, cluster.name = "scvi_clusters_0.75")
HHA<- FindClusters(HHA, resolution = 0.6, cluster.name = "scvi_clusters_0.6")
HHA<- FindClusters(HHA, resolution = 0.5, cluster.name = "scvi_clusters_0.5")
HHA<- FindClusters(HHA, resolution = 0.4, cluster.name = "scvi_clusters_0.4")
HHA<- FindClusters(HHA, resolution = 0.3, cluster.name = "scvi_clusters_0.3")

#Plotting Clusters
DimPlot(HHA, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.8", raster = F, pt.size = 0.01) + 
  NoAxes() + theme(aspect.ratio = 1) 
#Plotting by Sample ID
DimPlot(HHA, reduction = "UMAP", label = F, shuffle = T, group.by = "SampleID", raster = F,
        cols = MetBrewer::met.brewer("Signac", 11), pt.size = 0.1) + theme(aspect.ratio = 1) + NoAxes()
#Plotting by Isolation Method
HHA$Isolation <- factor(HHA$Isolation, levels = c("4mm", "FUE", "Bulb"))
DimPlot(HHA, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", split.by = "Isolation", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()
#Plotting by Platform
DimPlot(HHA, reduction = "UMAP", label = T, group.by = "scvi_clusters_0.75", split.by = "Platform", raster = F, pt.size = 0.01) & NoLegend() & NoAxes()
DimPlot(HHA, reduction = "UMAP", group.by = "Platform", raster = F, pt.size = 0.01, shuffle = T) + 
  NoAxes() + theme(aspect.ratio = 1)
#Plotting by Sex
DimPlot(HHA, reduction = "UMAP", group.by = "Sex", raster = F, pt.size = 0.01, shuffle = T) + 
  NoAxes() + theme(aspect.ratio = 1)

#-----FeaturePlots for General Clusters-----#
DefaultAssay(HHA) <- "RNA"
#QC Metrics
PlotGenesViridis(HHA, genes = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                 viridis_option = "inferno", raster = F) & NoAxes() & theme(aspect.ratio = 1)
#General Lineage Markers
PlotGenesViridis(HHA, genes = c("KRT5", "VIM", "PTPRC", "CD3G", "PDGFRA",
                                "PECAM1", "ACTA2", "DCT", "MSX2"), ncol = 3,
                 viridis_option = "inferno", raster = F)
#Melanocytes
PlotGenesViridis(HHA, genes = c("DCT", "PMEL", "TYRP1", "SOX10"), viridis_option = "inferno")
#Professional Antigen Presenting (Dendritic Cells, Langerhans Cells, Monocytes/Macrophages)
PlotGenesViridis(HHA, genes = c("CD74", "CD207", "CD14", "CD68"), viridis_option = "inferno", raster = T)
#Endothelial Cells
PlotGenesViridis(HHA, genes = c("PECAM1", "VWF", "CDH5", "PROX1"), raster = F, viridis_option = "inferno")
#Dermal Fibroblasts and Dermal Papilla
PlotGenesViridis(HHA, genes = c("COL1A1", "PDGFRA", "ITGA9", "CRABP1", "SOX2", "DCN"), ncol = 3, raster = F, viridis_option = "inferno")
#T Cells
PlotGenesViridis(HHA, genes = c("CD3G", "NKG7", "FOXP3", "CD4"),viridis_option = "turbo", raster = F)
#B Cells and Mast Cells
PlotGenesViridis(HHA, genes = c("IGHG1", "CD79A", "FCER1G", "TPSB2"),viridis_option = "plasma", raster = "F")
#Smooth Muscle
PlotGenesViridis(HHA, genes = c("RGS5", 'DES', "ACTA2", "ITGA8", 'ACTA1'), viridis_option = "plasma", raster = F)
#Schwann Cell Markers
PlotGenesViridis(HHA, genes = c("MPZ"), raster = F)
#Eccrine
PlotGenesViridis(HHA, genes = c("KRT7", "KRT18", "SAA1", "AQP5"), raster = F)
#Melanocytes
PlotGenesViridis(HHA, genes = c("MSX2", "KRT35", "GATA3", "KRT75"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("MKI67"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("CST6"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("S100A8", "KRT6A", "MGST1", "CST6"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("IGFBP3", "S100A8", "KRT1", "MGST1"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("KRT75"), raster = F, viridis_option = "inferno")
PlotGenesViridis(HHA, genes = c("RGS5"), raster = F, viridis_option = "inferno")

#-----Finding Markers for Each Cluster-----#
GeneralMarkers.TFIDF <- HHA %>% JoinLayers() %>% GetAssayData(assay = "RNA", slot = "counts") %>%
  SoupX::quickMarkers(clusters = HHA$scvi_clusters_0.8, N = 20, FDR = 0.01)
GeneralMarkers.TFIDF
GeneralMarkers.TFIDFList <- GeneralMarkers.TFIDF %>% mutate(cluster = as.numeric(cluster)) %>% group_by(cluster) %>% arrange(desc(tfidf))%>%
  select(cluster, gene, geneFrequency, geneFrequencySecondBest, geneFrequencyOutsideCluster, geneFrequencyGlobal, secondBestClusterName, tfidf,   idf,  qval) %>%  group_split() 
names(GeneralMarkers.TFIDFList) <- str_c("Cluster ", 0:(length(unique(HHA$scvi_clusters_0.8))-1))
GeneralMarkers.TFIDFList[15]

#------Assigning General ClusterIDs-----#
InitialAnnotation <- c(
  "Bulge", #0
  "Cortex/Cuticle", #1
  "Suprabasal ORS", #2
  "Basal Infundibulum", #3
  "Spinous/Granular I", #4
  "Basal IFE", #5
  "Fibroblasts I", #6
  "T Helper Cells", #7
  "SM/Dermal Sheath", #8
  "Fibroblasts II", #9
  "IRS/LPC", #10
  "Basal ORS I", #11
  "Pre-Sebaceous", #12
  "Transitional Infundibulum", #13
  "Quiescent ORS", #14
  "Antigen Presenting I", #15
  "Basal ORS II", #16
  "Upper Companion", #17
  "Endothelial", #18
  "Spinous/Granular II", #19
  "Melanocytes I", #20
  "Proliferating Basal", #21
  "T Reg", #22
  "NK Cells", #23
  "Proliferating ORS", #24
  "Plasma B Cells", #25
  "Pericytes", #26
  "Langerhans Cells", #27
  "Antigen Presenting II", #28
  "Schwann I", #29
  "Eccrine Epithelial", #30
  "Melanocytes II", #31
  "Eccrine Myoepithelial", #32
  "Mast Cells", #33
  "Antigen Presenting III", #34
  "Pro-B Cells", #35
  "Fibroblasts III", #36
  "Middle Companion", #37
  "Arrector Pili", #38
  "Lymph Vasculature", #39
  "Schwann II", #40
  "Dermal Papilla", #41
  "Merkel Cells" #42
)
#Renaming Idents
Idents(HHA) <- "scvi_clusters_0.8"
names(InitialAnnotation) <- levels(HHA)
HHA <- RenameIdents(HHA, InitialAnnotation)
HHA$InitialAnnotation<- HHA@active.ident
DimPlot(HHA, reduction = "UMAP", label = T, repel = T) & NoAxes() & theme(aspect.ratio = 1) + NoLegend()
DimPlot(HHA, reduction = "UMAP", label = T, repel = T, split.by = "Isolation")  & NoLegend() & NoAxes() 
table(HHA$InitialAnnotation)

#------Plotting Sample Characteristics by Cluster-----#
#Plotting by Sample
HHA[[]] %>% ggplot(aes(x = "", fill = SampleID)) + geom_bar(position = "fill") + facet_wrap(~InitialAnnotation) + 
  theme_bw() + scale_fill_manual(values = MetBrewer::met.brewer("Signac", 12)) + labs(x = "", y = "Proportion")
#Plotting by Donor
HHA[[]] %>% ggplot(aes(x = "", fill = DonorID)) + geom_bar(position = "fill") + facet_wrap(~InitialAnnotation) + 
  theme_bw() + scale_fill_manual(values = MetBrewer::met.brewer("Signac", 8))
#Plotting by Sex
HHA[[]] %>% ggplot(aes(x = "", fill = Sex)) + geom_bar(position = "fill") + 
  facet_wrap(~InitialAnnotation) + theme_bw()
#Plotting by Isolation Method
HHA[[]] %>% ggplot(aes(x = "", fill = Isolation)) + geom_bar(position = "fill") + facet_wrap(~InitialAnnotation) + theme_bw()
#Plotting by Platform
HHA[[]] %>% ggplot(aes(x = "", fill = Platform)) + geom_bar(position = "fill") + facet_wrap(~InitialAnnotation) + theme_bw()

HHA$InitialAnnotationOrdered <- factor(HHA$InitialAnnotation, levels = c("Cortex/Cuticle", "IRS/LPC", "Bulge", 'Basal ORS I', "Basal ORS II", "Proliferating ORS", "Quiescent ORS", "Suprabasal ORS", "Middle Companion", "Upper Companion",
                              "Basal Infundibulum", "Transitional Infundibulum", "Pre-Sebaceous", "Basal IFE", "Proliferating Basal", "Spinous/Granular I", "Spinous/Granular II",
                              "Merkel Cells", "Eccrine Epithelial",  "Eccrine Myoepithelial", "Arrector Pili", "SM/Dermal Sheath", "Pericytes", "Endothelial", "Lymph Vasculature",
                              "Fibroblasts I", "Fibroblasts II", "Fibroblasts III", "Dermal Papilla", "Schwann I", "Schwann II",  "Melanocytes I",  "Melanocytes II",
                              "Antigen Presenting I", "Antigen Presenting II", "Antigen Presenting III", "Langerhans Cells",
                              'Mast Cells', "T Helper Cells", "T Reg", "NK Cells", "Pro-B Cells", "Plasma B Cells"))
DimPlot(HHA, reduction = "UMAP", group.by = "InitialAnnotationOrdered", label = T, repel = T)  & NoLegend() & NoAxes() & theme(aspect.ratio = 1)
#Plotting Dot Plot
DotPlot(HHA,  group.by = "InitialAnnotationOrdered", features = c("MSX2", "KRT35", "KRT25", "KRT28", 'KRT5', "SOX9", "DIO2", "TCEAL2","LGR5", "CTNND2", "BARX2", "KRT6A", "KRT75", "S100A8", "S100A9", "SAA1", "IGFBP3", "WNT3", "KRT1", 'KRT20', "KRT18", "KRT7",
                                                          "CHRM3", "DES", "ACTA1", "ACTA2", "ITGA8",  "RGS5", "PECAM1", "CDH5", "PROX1", "CCL21", "DCN", "PDGFRA", "RSPO3", "CDH19", "NRXN1", "MPZ", "DCT", "PMEL",
                                                          "PTPRC", "FCER1G", "CD14", 'CD68', "CD207", "TPSB2", "CD3D", "CD3G", "FOXP3", "NKG7", "MS4A1", "CD79A", "MZB1",
                                                          "IGHG1", 'VIM')) + RotatedAxis() + scale_color_gradientn(colors = c("black","#2488F0","#7F3F98","#E22929", "#FCB31A")) + 
  theme(aspect.ratio = 1, panel.background = element_rect(color = "black"))

#Saving image
save.image(paste0(WorkspaceDirectory, "HHA_Integration_Initial_Annotation_11_1_24.Rdata"))
#Saving object
save(HHA, file = paste0(SeuratDirectory, "HHA_SCVI_Annotated_10_26_24.rds"))
#-----Recording session information-----#
sessionInfo()
#writing to text file
writeLines(capture.output(sessionInfo()), "/home/cms2775/HHA_Spatial_scRNAseq_2024/sessionInfo/HHA_PreIntegration_DR_Clustering_11_1_24_sessionInfo.txt")
