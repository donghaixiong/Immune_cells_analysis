

############################################################################################################################################################################################################################################

# R v4.0.2

rm(list=ls())

library(Seurat)

library(BiocParallel)

library(Biobase)
library(GenomicRanges)

library(SummarizedExperiment)

library(DESeq2)


library(rocc)

library(Biobase)

library(edgeR)

library(limma)


library("RColorBrewer")

library(biomaRt)

library(org.Hs.eg.db)

library(annotate)

library(plyr)

library(ggplot2)

library(cowplot)

library(dplyr)

library(cancerclass)


setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes')

### the following datafile 'ImmuneCells23clusters.rds' is big and deposited in GitLab instead of GitHub 
### The download address of the 'ImmuneCells23clusters.rds' file is: https://gitlab.com/donghaixiong/immune.cells.big.datafiles/-/blob/master/ImmuneCells23clusters.rds

immune.combined.more.clusters <- readRDS('ImmuneCells23clusters.rds')

head(immune.combined.more.clusters@meta.data)

table(Idents(immune.combined.more.clusters))

sum(table(Idents(immune.combined.more.clusters)))

MM_cluster6_12_23 <- subset(immune.combined.more.clusters, idents = c("6","12","23"))

str(MM_cluster6_12_23@meta.data)

dim(MM_cluster6_12_23@meta.data)

table(Idents(MM_cluster6_12_23))

sort(unique(MM_cluster6_12_23@meta.data$Cell.Type))


setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes/output_test')

# Figure 2a

featured.markers <- c("IDO1", "FCER1A", "NAPSB", "S100A12", "APOBEC3A", "CD1C", "SELL", "CXCL10", "LGALS2",         
"TREM2", "SPP1", "RNASE1", "MT1G", "SEPP1", "FOLR2", "NUPR1", "KLHDC8B", "CCL18", "MMP12", "APOC2", "C3", "C1QA", "C1QB", "C1QC",
"LCK", "TRAC", "CD3D", "CD3G", "TIGIT", "PTPRCAP", "KLRK1", "CD8A", "PRF1", "CD2", "LAT", "IL32", "CST7", "NKG7", "IFITM1", "CCL5", "CD69", "CD8B", "CD96")

tiff("Fig2a.Macrophage.3clusters.withLeg.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = featured.markers, hjust = 0.5, angle = 0, slot = "scale.data")

dev.off()


# Figure 2b


### VlnPlot for IDO1.markers of cluster6

IDO1.markers <- c('IDO1', 'APOBEC3A', 'CXCL10')

tiff("MM_3clusters.VlnPlot.IDO1.markers.tiff", width = 5000, height = 2500, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23, features = IDO1.markers, pt.size = 0.02, log = T)

dev.off()



### VlnPlot for TREM2.markers of cluster12

TREM2.markers <- c('TREM2', 'RNASE1', 'MT1G', 'SEPP1', 'FOLR2', 'NUPR1')

tiff("MM_3clusters.VlnPlot.TREM2.markers.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23, features = TREM2.markers, pt.size = 0.02, log = T) 

dev.off()



### VlnPlot for Immunoregulatory.markers of cluster6

Immunoregulatory.markers <- c('TRAC', 'KLRK1', 'CD96')

tiff("MM_3clusters.VlnPlot.Immunoregulatory.markers.tiff", width = 5000, height = 2500, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23, features = Immunoregulatory.markers, pt.size = 0.02, log = T) 

dev.off()

