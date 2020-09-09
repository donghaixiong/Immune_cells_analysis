
############################################################################################################################################################################################################################################

# R v4.0.2

rm(list=ls())

####### o.rocc

#library(DEsingle)

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

#library(RUVSeq)

#library(EDASeq)

library(org.Hs.eg.db)

library(annotate)

library(plyr)

library(ggplot2)

library(Seurat)

library(cowplot)

library(dplyr)

library(cancerclass)

human = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "oct2016.archive.ensembl.org")
mouse = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "oct2016.archive.ensembl.org")


setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes')

ImmuneCells23clusters <- readRDS('ImmuneCells23clusters.rds')

table(Idents(ImmuneCells23clusters))


tiff(file = "Fig1a_ImmuneCells_Clustering_condition_sidebyside.tiff", width = 8000, height = 8000, units = "px", res = 800)

DimPlot(object = ImmuneCells23clusters, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 4.5)

dev.off()


pdf(file = "Figure1a.pdf", width = 8, height = 8)

DimPlot(object = ImmuneCells23clusters, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 4.5)

dev.off()



new.ids2 <- c("CD8+ T cells", "Regulatory T cells", "CD4+ T cells", "CD8+ T cells", "CD8+ T cells", "Macrophages", 
"CD8+ T cells", "Tgd cells", "MKI67hi", "CD8+ T cells", "CD8+ T cells", "Macrophages", "B cells", "B cells", 
"NK cells", "MKI67hi", "B cells", "Plasma cells", "Dendritic cells", "CD8+ T cells", "Tgd cells", "B cells", "Macrophages")

names(new.ids2) <- as.character(1:23)

new.ids2



ImmuneCells23clusters <- RenameIdents(ImmuneCells23clusters, new.ids2)

ImmuneCells23clusters[['Cell.Type']] <- Idents(ImmuneCells23clusters)

table(Idents(ImmuneCells23clusters))

str(ImmuneCells23clusters@meta.data)


pdf(file = "Figure1b.pdf", width = 8, height = 8)

DimPlot(object = ImmuneCells23clusters, reduction = "umap", label = TRUE, label.size = 4.5)

dev.off()





