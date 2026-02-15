#-----Loading in Packages-----#
library(readr)
library(dplyr)
library(ggplot2)
library(ArchR)
library(purrr)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86) 
library(ArchRtoSignac) 
library(Signac)


#-----Directories------#
SeuratDirectory <- "/projects/b1217/HHA/Multiome_Seurat/"
PlotsDirectory <- "/projects/b1217/HHA/Bulb_Seurat_Plots/Figure_6/"
WorkspaceDirectory <- "/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Workspaces/" 
MetadataDirectory <- "/projects/b1217/HHA/Multiome_Metadata/"

#---------Loading in archR object----------#
Matrix.archr <- loadArchRProject("/projects/b1217/HHA/Multiome_ArchR/Multiome_Matrix_5_30_25")
Matrix.archr

#-------Loading in RNA Seurat Object------#
load(paste0(SeuratDirectory, "HHA_Multiome_Matrix_Integrated_Annotated_5_30_25.rds"))
Matrix.Multiome

#-----Getting Annotations-------#
##ENSEMBL v107
annotations <- rtracklayer::import('/projects/b1217/HHA/Dictys/Inputs/gene.gtf')
#Removing Non-Standard Chromosomes
annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
#Restandardizing seqnames
seqlevelsStyle(annotations) <- "UCSC"
annotations$tx_id <- annotations$transcript_id

#-----Getting Peak Matrix------#
pkm <- getPeakMatrix(Matrix.archr)

#---------Fragment Files for Each Sample-------#
samplelist <- list("EL_C5", "EL_C6", "EL_E6", "EL_E7", "EL_E8")
fragments_dir <- list("/projects/b1217/HHA/Multiomics_frag/EL_C5/",
                      "/projects/b1217/HHA/Multiomics_frag/EL_C6/",
                      "/projects/b1217/HHA/Multiomics_frag/EL_E6/",
                      "/projects/b1217/HHA/Multiomics_frag/EL_E7/",
                      "/projects/b1217/HHA/Multiomics_frag/EL_E8/")

#-----Converting to Signac------#
Matrix.ATAC <- ArchR2Signac(
  ArchRProject = Matrix.archr,
  refversion = "hg38",
  samples = samplelist,
  fragments_dir = fragments_dir,
  pm = pkm,
  fragments_fromcellranger = "NO",
  fragments_file_extension = 'atac_fragments.tsv.gz',
  annotation = annotations)

#------Adding Chromatin Data to Multiome Seurat Object-------$
#Resetting colnames to UnifiedBarcode
colnames(Matrix.Multiome) <- Matrix.Multiome$UnifiedBarcode
#Adding Assay
Matrix.Multiome[["peaks"]] <- Matrix.ATAC[["peaks"]]
Matrix.Multiome

#--------Testing Workflow------#
DefaultAssay(Matrix.Multiome) <- "peaks"
Idents(Matrix.Multiome) <- "MatrixAnnotationBroad"
CoveragePlot(
  object = Matrix.Multiome,
  region = c("KRT35"),
  idents = c("Cortex", "Cuticle", "IRS"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1)

#-----Loading in Significant Peak2Gene Links-------#
Matrix.P2G <- read_csv("/projects/b1217/HHA/Matrix_ATAC_SEACells/Matrix_HRG_PG_Links_Stress_CC_Filtered_6_5_25.csv")
#Converting Peak names to match signac format
Matrix.P2G$peak <- str_replace(Matrix.P2G$peak, pattern = ":", replacement = "-")
Matrix.P2G

#------Extracting TSS Positions-------#
Genes.TSS <- GetTSSPositions(annotations, biotypes = NULL)
Genes.TSS
P2G.TSS <- Matrix.P2G %>% dplyr::filter(gene %in% Genes.TSS$gene_name) %>% left_join(data.frame(gene = Genes.TSS$gene_name, TSS = start(Genes.TSS@ranges)))
P2G.TSS

#---------Building GRanges Object--------#
#Only using protein coding genes
P2G.use <- Matrix.P2G %>% dplyr::filter(gene %in% Genes.TSS$gene_name)
#get midpoint of each peak
P2G.PeakRanges <- StringToGRanges(
  regions = P2G.use$peak,
  sep = c("-", "-"))
P2G.midpoints <- start(x = P2G.PeakRanges) + (width(x = P2G.PeakRanges) / 2)
# create dataframe
P2G.df <- data.frame(
  chromosome = as.character(x = seqnames(P2G.PeakRanges)),
  tss = P2G.TSS$TSS,
  pk = P2G.midpoints,
  score = P2G.use$cor,
  gene = P2G.use$gene,
  peak = P2G.use$peak)
#work out start and end coords
P2G.df$start <- ifelse(test = P2G.df$tss < P2G.df$pk, yes = P2G.df$tss, no = P2G.df$pk)
P2G.df$end <- ifelse(test = P2G.df$tss < P2G.df$pk, yes = P2G.df$pk, no = P2G.df$tss)
P2G.df$tss <- NULL
P2G.df$pk <- NULL
#convert to granges
P2G.gr <- sort(x = makeGRangesFromDataFrame(df = P2G.df, keep.extra.columns = TRUE))
P2G.gr

#-------Adding Links to Object---------#
DefaultAssay(Matrix.Multiome) <- "peaks"
Links(Matrix.Multiome) <- P2G.gr

#--------Testing Workflow------#
DefaultAssay(Matrix.Multiome) <- "peaks"
Idents(Matrix.Multiome) <- "MatrixAnnotationFine"
CoveragePlot(
  object = Matrix.Multiome,
  region = c("KRT36"),
  idents = c("Early_Cortex", "Middle_Cortex", "Late_Cortex"),
  extend.upstream = 10000,
  extend.downstream = 10000,
  ncol = 1)

#-------Saving Objects---------#
save(Matrix.ATAC, file = paste0(SeuratDirectory, "HHA_Matrix_Multiome_ATAC_Signac_6_16_25.rds"))
save(Matrix.Multiome, file = paste0(SeuratDirectory, "HHA_Multiome_ATAC_Matrix_Integrated_Annotated_6_16_25.rds"))

