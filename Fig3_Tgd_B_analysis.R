
# using R x64 3.5.3

rm(list=ls())

####### o.rocc

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

library(RUVSeq)

library(EDASeq)

library(cancerclass)

library(org.Hs.eg.db)

library(annotate)

library(plyr)

library(ggplot2)

library(Seurat)

library(cowplot)

library(dplyr)

library(fgsea)

library(data.table) 

library(ggplot2)

library('org.Hs.eg.db')

library(tidyverse)


setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes')

immune.combined.more.clusters <- readRDS('ImmuneCells23clusters.rds')

table(Idents(immune.combined.more.clusters))

head(rownames(immune.combined.more.clusters))

length(rownames(immune.combined.more.clusters))


### Tgd cells data analyses

Tgd_clusters <- subset(immune.combined.more.clusters, idents = c("8","21"))

# Run the standard workflow for visualization and clustering

Tgd_clusters <- ScaleData(object = Tgd_clusters, features = rownames(Tgd_clusters))

DefaultAssay(object = Tgd_clusters)


# UMAP and Clustering

Tgd_clusters <- RunUMAP(object = Tgd_clusters, reduction = "pca", dims = 1:20)

Tgd_clusters <- FindNeighbors(object = Tgd_clusters, reduction = "pca", dims = 1:20)


str(Tgd_clusters@meta.data)

dim(Tgd_clusters@meta.data)

sort(unique(Tgd_clusters@meta.data$Cell.Type))


setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes/Tgd_B')

DefaultAssay(object = Tgd_clusters) <- "RNA"

table(Idents(Tgd_clusters))


# find markers for every cluster compared to all remaining cells

Tgd_clusters.markers <- FindAllMarkers(Tgd_clusters, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.4)  

top20 <- Tgd_clusters.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

top20$gene



Tgd_21vs8.DE <- FindMarkers(object = Tgd_clusters, ident.1 = "21", ident.2 = "8", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

Tgd_21vs8.DE <- cbind(rownames(Tgd_21vs8.DE),Tgd_21vs8.DE)

colnames(Tgd_21vs8.DE)[1] <- 'Gene'

Tgd_21vs8.DE <- Tgd_21vs8.DE[order(Tgd_21vs8.DE$p_val_adj),]



Tgd_21vs8.DE$FoldChange <- Tgd_21vs8.DE$pct.1/Tgd_21vs8.DE$pct.2

Tgd_21vs8.DE <- Tgd_21vs8.DE[order(Tgd_21vs8.DE$p_val_adj, -Tgd_21vs8.DE$FoldChange),]

range(Tgd_21vs8.DE$p_val)


Tgd_21vs8.DE.topFCgenes <- rownames(Tgd_21vs8.DE[Tgd_21vs8.DE$FoldChange > 10 & Tgd_21vs8.DE$p_val_adj < 0.01 & !is.infinite(Tgd_21vs8.DE$FoldChange) & !is.na(Tgd_21vs8.DE$FoldChange),])

length(Tgd_21vs8.DE.topFCgenes)

head(Tgd_21vs8.DE.topFCgenes)

selected_Tgd_markers_forplot_V2 <- c(top20$gene[1:20],Tgd_21vs8.DE.topFCgenes[1:20])


tiff("Tgd_clusters.heatmap.BWR_V2.tiff", width = 4500, height = 4000, units = "px", res = 600)
  
DoHeatmap(Tgd_clusters, features = selected_Tgd_markers_forplot_V2, hjust =1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

dev.off()


#############################################################################################################################################


### B cells data analyses

B_clusters <- subset(immune.combined.more.clusters, idents = c("13","14","17","22"))

# Run the standard workflow for visualization and clustering

B_clusters <- ScaleData(object = B_clusters, features = rownames(B_clusters))

DefaultAssay(object = B_clusters)


B_clusters <- RunUMAP(object = B_clusters, reduction = "pca", dims = 1:20)

B_clusters <- FindNeighbors(object = B_clusters, reduction = "pca", dims = 1:20)


str(B_clusters@meta.data)

dim(B_clusters@meta.data)

sort(unique(B_clusters@meta.data$Cell.Type))


DefaultAssay(object = B_clusters) <- "RNA"

table(Idents(B_clusters))


# find markers for every cluster compared to all remaining cells

B_clusters.markers <- FindAllMarkers(B_clusters, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.4)  

top20 <- B_clusters.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

top20$gene

selected_B_markers_forplot <- top20$gene


tiff("B_clusters.heatmap.BWR.tiff", width = 6000, height = 6000, units = "px", res = 600)
  
DoHeatmap(B_clusters, features = selected_B_markers_forplot, hjust =1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) 

dev.off()


#############################################################################################################################################


### Tgd_B significant pathways

Tidy.Tgd.c21.sigOnly.short.DP.rev <- read.table('Tidy.Tgd.c21.sigOnly.short.DP.rev.txt',sep='\t',header=T)

dim(Tidy.Tgd.c21.sigOnly.short.DP.rev)

head(Tidy.Tgd.c21.sigOnly.short.DP.rev)


Tidy.B.c22.sigOnly.short.DP.rev <- read.table('Tidy.B.c22.sigOnly.short.DP.rev.txt',sep='\t',header=T)

dim(Tidy.B.c22.sigOnly.short.DP.rev)

head(Tidy.B.c22.sigOnly.short.DP.rev)


Tgd.B.2clusters.selected.pathways <- unique(c(as.character(Tidy.Tgd.c21.sigOnly.short.DP.rev$pathway),as.character(Tidy.B.c22.sigOnly.short.DP.rev$pathway)))

length(Tgd.B.2clusters.selected.pathways)

Tgd.B.2clusters.selected.pathways


Tidy.Tgd.c21.sigOnly.short.DP <- read.table('Tidy.Tgd.c21.sigOnly.short.DP.xls',sep='\t',header=T)
Tidy.B.c22.sigOnly.short.DP <- read.table('Tidy.B.c22.sigOnly.short.DP.xls',sep='\t',header=T)

Tidy.Tgd.c21.B.c22.sigOnly.short.DP <- rbind(Tidy.Tgd.c21.sigOnly.short.DP, Tidy.B.c22.sigOnly.short.DP)


Tgd_B_combined <- Tidy.Tgd.c21.B.c22.sigOnly.short.DP

Tgd_B_combined <- Tgd_B_combined[Tgd_B_combined$pathway %in% Tgd.B.2clusters.selected.pathways,]


Tgd_B_combined$cluster <- factor(Tgd_B_combined$cluster, levels=c("Tgd_c21", "B_c22"))


tiff(file = "Tgd_B_sig_pathways.tiff", width = 8000, height = 3500, units = "px", res = 600)

p <- ggplot(Tgd_B_combined, aes(cluster, pathway))

p + geom_point(aes(color=FDR, size=NES.abs, shape=Direction)) +
		scale_color_gradientn(colours=rainbow(4), limits=c(0, 0.05)) +
		scale_shape_manual(values=c(19,17))+
    theme(panel.background=element_rect(fill="gray95", colour="gray95"),
          panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
          axis.text.x = element_text(colour="black",angle=0,hjust=0.5,size=12),
        	axis.text.y = element_text(colour="black",size=12), 
        	axis.title.x = element_text(size=8),
          legend.title = element_text(size=8),
          legend.text = element_text(size = 8),
          axis.title.y=element_blank()) 
    
dev.off()
 

Tgd_B_combined.out <- data.frame(Tgd_B_combined)

class(Tgd_B_combined.out)

dim(Tgd_B_combined.out)

head(Tgd_B_combined.out)

tail(Tgd_B_combined.out)

write.table(Tgd_B_combined.out, 'Tgd_B_combined.out.xls', sep='\t', row.names=F, quote=F)


###################################################################################################################################################################################################
###################################################################################################################################################################################################
###################################################################################################################################################################################################


