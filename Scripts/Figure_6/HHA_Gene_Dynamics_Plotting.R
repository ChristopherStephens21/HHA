#-----Loading Packages-----#
library(Seurat)
library(scCustomize)
library(tidyverse)

#-------Loading in Motif Dynamics Data------#
#DORCs
Combined.NormDORCs <- read_csv("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/DORC_Trends.csv")
#Genes
Combined.NormGenes <- read_csv("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/Gene_Trends.csv")
#Cortex Lags
Cortex.Plotting <- read_csv("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/Cortex_Gene_Dynamics_Trends.csv")
#IRS Cuticle Lags
IRSCuticle.Plotting <- read_csv("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/Figure_6/Gene_Trends_Data/IRS_Cuticle_Gene_Dynamics_Trends.csv")

#--------Plotting Combined Results-------#
#---------GATA3---------#
#Plotting DORC Accesibility
GATA3.DORCPlot <- ggplot() + 
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.5) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.8) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  theme_bw() +
  ggtitle("DORC Accessibility") +
  xlab("Pseudotime") +
  ylab("GATA3")
GATA3.DORCPlot
#Plotting Gene Expression 
GATA3.GEXPlot <- ggplot() + 
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.5) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.8) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = GATA3),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  theme_bw() +
  ggtitle("Gene Expression") +
  xlab("Pseudotime") +
  ylab("GATA3")
GATA3.GEXPlot
#Plotting Chromatin-RNA latency
GATA3.LagPlot <- ggplot(IRSCuticle.Plotting, aes(x = Pseudotime, y = GATA3_Lag)) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "darkseagreen4", alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "darkseagreen4", size = 2) +
  theme_bw() +
  labs(title = "Chromatin-RNA Latency", y = "GATA3", x = "Pseudotime")
GATA3.LagPlot

#---------KRT36---------#
#Plotting DORC Accessibility
KRT36.DORCPlot <- ggplot() + 
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.5) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.8) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  theme_bw() +
  ggtitle("DORC Accessibility") +
  xlab("Pseudotime") +
  ylab("KRT36")
KRT36.DORCPlot
#Plotting Gene Expression
KRT36.GEXPlot <- ggplot() + 
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.5) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.8) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = KRT36),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  theme_bw() +
  ggtitle("Gene Expression") +
  xlab("Pseudotime") +
  ylab("KRT36")
KRT36.GEXPlot
#Plotting Chromatin-RNA Latency
KRT36.LagPlot <- ggplot(Cortex.Plotting, aes(x = Pseudotime, y = KRT36_Lag)) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "darkseagreen4", alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "darkseagreen4", size = 2) +
  theme_bw() +
  labs(title = "Chromatin-RNA Latency", y = "KRT36", x = "Pseudotime")
KRT36.LagPlot

#---------TP63---------#
#Plotting DORC Accessibility
TP63.DORCPlot <- ggplot() + 
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.5) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.8) +
  stat_smooth(data = Combined.NormDORCs %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  theme_bw() +
  ggtitle("DORC Accessibility") +
  xlab("Pseudotime") +
  ylab("TP63")
TP63.DORCPlot
#Plotting Gene Expression
TP63.GEXPlot <- ggplot() + 
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "firebrick", alpha = 0.5) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "IRS Cuticle"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "firebrick", size = 2) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "cornflowerblue", alpha = 0.8) +
  stat_smooth(data = Combined.NormGenes %>% filter(Trajectory == "Cortex"),
              aes(x = Pseudotime, y = TP63),
              method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "cornflowerblue", size = 2) +
  theme_bw() +
  ggtitle("Gene Expression") +
  xlab("Pseudotime") +
  ylab("TP63")
TP63.GEXPlot
#Plotting Chromatin-RNA Latency: Cortex
TP63.LagPlot <- ggplot(Cortex.Plotting, aes(x = Pseudotime, y = TP63_Lag)) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "darkseagreen4", alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "darkseagreen4", size = 2) +
  theme_bw() +
  labs(title = "Chromatin-RNA Latency", y = "TP63", x = "Pseudotime")
TP63.LagPlot
#Plotting Chromatin-RNA Latency: IRS
TP63.IRSLagPlot <- ggplot(IRSCuticle.Plotting, aes(x = Pseudotime, y = TP63_Lag)) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "area", fill = "darkseagreen4", alpha = 0.5) +
  stat_smooth(method = "loess", span = 0.1, se = FALSE,
              geom = "line", color = "darkseagreen4", size = 2) +
  theme_bw() +
  labs(title = "Chromatin-RNA Latency", y = "TP63", x = "Pseudotime")
TP63.IRSLagPlot

#-------Saving as tiff-----#
PlotsDirectory <- '/projects/b1217/Edward/R_Projects/HHA/Figures/Figure_06/'
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_DORC.tiff"), res = 300, height = 6, width = 6, units = "in") 
GATA3.DORCPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_GEX.tiff"), res = 300, height = 6, width = 6, units = "in") 
GATA3.GEXPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_Lag.tiff"), res = 300, height = 6, width = 6, units = "in") 
GATA3.LagPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "GATA3_Motif.tiff"), res = 300, height = 6, width = 6, units = "in") 
GATA3.MotifPlot
dev.off()

#-------Saving as tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "TP63_DORC.tiff"), res = 300, height = 6, width = 6, units = "in") 
TP63.DORCPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "TP63_GEX.tiff"), res = 300, height = 6, width = 6, units = "in") 
TP63.GEXPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "TP63_Lag.tiff"), res = 300, height = 6, width = 6, units = "in") 
TP63.LagPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "TP63_Motif.tiff"), res = 300, height = 6, width = 6, units = "in") 
TP63.MotifPlot
dev.off()
ragg::agg_tiff(filename = paste0("/projects/b1042/YiLab/HHA_scRNA_Spatial_Shared/TP63_Plots/", "TP63_Lag.tiff"), res = 300, height = 6, width = 6, units = "in") 
TP63.IRSLagPlot
dev.off()

#-------Saving as tiff-----#
ragg::agg_tiff(filename = paste0(PlotsDirectory, "KRT36_DORC.tiff"), res = 300, height = 6, width = 6, units = "in") 
KRT36.DORCPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "KRT36_GEX.tiff"), res = 300, height = 6, width = 6, units = "in") 
KRT36.GEXPlot
dev.off()
ragg::agg_tiff(filename = paste0(PlotsDirectory, "KRT36_Lag.tiff"), res = 300, height = 6, width = 6, units = "in") 
KRT36.LagPlot
dev.off()








