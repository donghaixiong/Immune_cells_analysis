
# using R x64 3.5.3

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

library(cancerclass)

library(org.Hs.eg.db)

library(annotate)

library(plyr)

library(ggplot2)

library(Seurat)

library(cowplot)

library(dplyr)


human = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "oct2016.archive.ensembl.org")
mouse = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "oct2016.archive.ensembl.org")


### note of 5/3/2020: cancerclass failed

#setwd("/rcc/temp_dxiong/DefineT_CellPaper_Analysis")


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

source("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/heatmap3.R")


# All CD45ImmuneCells

CD45ImmuneCells_tpm <- read.table('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt',sep='\t',fill=T,header=T)

dim(CD45ImmuneCells_tpm)

head(CD45ImmuneCells_tpm[,1:30])

tail(CD45ImmuneCells_tpm[,1:30])

head(CD45ImmuneCells_tpm[,16283:16292])

tail(CD45ImmuneCells_tpm[,16283:16292])


class(CD45ImmuneCells_tpm[1,])

colnames(CD45ImmuneCells_tpm)[1:16291] <- colnames(CD45ImmuneCells_tpm)[2:16292]

CD45ImmuneCells_tpm <- CD45ImmuneCells_tpm[-1,-16292]


class(CD45ImmuneCells_tpm)

#range(CD45ImmuneCells_tpm) # tab-delimited text file containing log2(TPM+1) values (transcript-per-million reads) for 55,737 transcripts (rows) and 16,291 cells (columns)



DefineT_PhenoInfo <- read.table("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/DefineT_PhenoInfo.txt",sep="\t",header=T)

DefineT_PhenoInfo <- DefineT_PhenoInfo[,-c(1,3,4)]

colnames(DefineT_PhenoInfo)[1] <- "Cell.Name"


dim(DefineT_PhenoInfo)

head(DefineT_PhenoInfo)

tail(DefineT_PhenoInfo)

unique(DefineT_PhenoInfo[[2]])

unique(DefineT_PhenoInfo[[4]])


DefineT_PhenoInfo$Cell.Name <- as.character(DefineT_PhenoInfo$Cell.Name)

class(DefineT_PhenoInfo$Cell.Name)


all(colnames(CD45ImmuneCells_tpm)==DefineT_PhenoInfo$Cell.Name)

length(intersect(colnames(CD45ImmuneCells_tpm),DefineT_PhenoInfo$Cell.Name))

setdiff(colnames(CD45ImmuneCells_tpm),DefineT_PhenoInfo$Cell.Name)

setdiff(DefineT_PhenoInfo$Cell.Name,colnames(CD45ImmuneCells_tpm))

unmatchedID <- data.frame(unmatchedExpVsPheno=setdiff(colnames(CD45ImmuneCells_tpm),DefineT_PhenoInfo$Cell.Name),
unmatchedPhenoVsExp=setdiff(DefineT_PhenoInfo$Cell.Name,colnames(CD45ImmuneCells_tpm)))

dim(unmatchedID)

head(unmatchedID)

tail(unmatchedID)



colnames(CD45ImmuneCells_tpm)[colnames(CD45ImmuneCells_tpm) %in% unmatchedID$unmatchedExpVsPheno] <- 
gsub("\\.","-",colnames(CD45ImmuneCells_tpm)[colnames(CD45ImmuneCells_tpm) %in% unmatchedID$unmatchedExpVsPheno])

all(colnames(CD45ImmuneCells_tpm)==DefineT_PhenoInfo$Cell.Name)

dim(CD45ImmuneCells_tpm)

head(CD45ImmuneCells_tpm[,1:30])

head(CD45ImmuneCells_tpm[,16282:16291])

tail(CD45ImmuneCells_tpm[,1:30])

tail(CD45ImmuneCells_tpm[,16282:16291])

#save(CD45ImmuneCells_tpm, file="CD45ImmuneCells_tpm.RData")



DefineT_ClusterInfo <- read.table("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/DefineT_ClusterInfo.txt",sep="\t",header=T)

DefineT_ClusterInfo <- DefineT_ClusterInfo[,1:2]

DefineT_ClusterInfo$Cell.Name <- gsub("_DN1","_T_enriched",DefineT_ClusterInfo$Cell.Name)

DefineT_ClusterInfo$Cell.Name <- gsub("_DN","_T_enriched",DefineT_ClusterInfo$Cell.Name)

DefineT_ClusterInfo$Cell.Name <- gsub("_DP1","_T_enriched",DefineT_ClusterInfo$Cell.Name)

dim(DefineT_ClusterInfo)

head(DefineT_ClusterInfo)

tail(DefineT_ClusterInfo)


length(intersect(DefineT_PhenoInfo$Cell.Name, DefineT_ClusterInfo$Cell.Name))

setdiff(DefineT_PhenoInfo$Cell.Name, DefineT_ClusterInfo$Cell.Name)

remainIDs <- DefineT_PhenoInfo$Cell.Name[c(15301:15418,16106:16291)]

all(remainIDs==setdiff(DefineT_PhenoInfo$Cell.Name, DefineT_ClusterInfo$Cell.Name))


remainIDs_trim <- gsub("_myeloid_enriched","",remainIDs)

remainIDs_trim <- gsub("_T_enriched","",remainIDs_trim)

all(remainIDs_trim==DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

remainIDs_trim; DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)]

class(DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

class(remainIDs)

remainIDs <- as.character(remainIDs)

DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)] <- remainIDs

all(remainIDs==DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

length(DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

head(DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

tail(DefineT_ClusterInfo$Cell.Name[c(15301:15418,16106:16291)])

length(intersect(DefineT_PhenoInfo$Cell.Name, DefineT_ClusterInfo$Cell.Name))


setdiff(DefineT_ClusterInfo$Cell.Name, DefineT_PhenoInfo$Cell.Name)

any(duplicated(DefineT_PhenoInfo$Cell.Name))

any(duplicated(DefineT_ClusterInfo$Cell.Name))


DefineT_PhenoInfo_V2 <- merge(DefineT_PhenoInfo,DefineT_ClusterInfo,by="Cell.Name")

colnames(DefineT_PhenoInfo_V2)[2:4] <- c("SampleID","Response","Therapy")



sort(unique(DefineT_PhenoInfo_V2$Cluster.number))

class(DefineT_PhenoInfo_V2$Cluster.number)

DefineT_PhenoInfo_V2$TreatTime <- DefineT_PhenoInfo_V2$SampleID

DefineT_PhenoInfo_V2$TreatTime <- gsub('_.*$','',DefineT_PhenoInfo_V2$TreatTime)

table(DefineT_PhenoInfo_V2[[2]])

table(DefineT_PhenoInfo_V2[[3]])

table(DefineT_PhenoInfo_V2[[4]])

table(DefineT_PhenoInfo_V2$TreatTime)


rownames(DefineT_PhenoInfo_V2) <- DefineT_PhenoInfo_V2$Cell.Name

dim(DefineT_PhenoInfo_V2)

head(DefineT_PhenoInfo_V2)

tail(DefineT_PhenoInfo_V2)

save(DefineT_PhenoInfo_V2, file="DefineT_PhenoInfo_V2.RData")


clustered.cells <- list()

cluster.names <- c('B.cells','Plasma.cells','Monocytes/Macrophages','Dendritic.cells','Lymphocytes','Exhausted.CD8.T.cells','Regulatory.T.cells',
'CytotoxicityLymphocytes','Exhausted/HS.CD8.T.cells','Memory.T.cells','Lymphocytes.Exhausted/Cell.cycle')

for (i in 1:11) {

		clustered.cells[[i]] <- as.character(DefineT_PhenoInfo_V2$Cell.Name[DefineT_PhenoInfo_V2$Cluster.number == i])
		
		names(clustered.cells)[i] <- cluster.names[i]
		
		print(length(clustered.cells[[i]]))

}

str(clustered.cells)



##################### use SingleR to assign cell types ##################### 

dim(CD45ImmuneCells_tpm)

head(CD45ImmuneCells_tpm[,1:30])

(cluster.names.V2 <- gsub("\\/",".",cluster.names))


### divide and output into 11 groups for local Seurat 2.3.4 analyses

CD45ImmuneCells_tpm_11groups <- list()

for (i in 1:11) {

		CD45ImmuneCells_tpm_11groups[[i]] <- CD45ImmuneCells_tpm[,clustered.cells[[i]]]
		
		colnumber <- ncol(CD45ImmuneCells_tpm_11groups[[i]])
		
		CD45ImmuneCells_tpm_11groups[[i]]$Gene <- rownames(CD45ImmuneCells_tpm_11groups[[i]])

		
		CD45ImmuneCells_tpm_11groups[[i]] <- CD45ImmuneCells_tpm_11groups[[i]][,c((colnumber+1),1:colnumber)]
		
		filename <- paste0("G",i,"_",cluster.names.V2[i],"_CD45ImmuneCells_tpm.txt")

		names(CD45ImmuneCells_tpm_11groups)[i] <- gsub("\\.txt","",filename)
		
		print(dim(CD45ImmuneCells_tpm_11groups[[i]]))
		
		print(head(CD45ImmuneCells_tpm_11groups[[i]][,1:4]))
		
		write.table(CD45ImmuneCells_tpm_11groups[[i]], filename, sep='\t', row.names=F, quote=F)

}

names(CD45ImmuneCells_tpm_11groups)


########### Line 284 to Line 322 omitted now

### divide and output into 20 groups for local Seurat 2.3.4 analyses


CD45ImmuneCells_tpm_20groups <- list()

for (i in 1:20) {

	if (i == 20) {
	
		CD45ImmuneCells_tpm_20groups[[i]] <- CD45ImmuneCells_tpm[,15201:16291]
		
	} else {
	
		CD45ImmuneCells_tpm_20groups[[i]] <- CD45ImmuneCells_tpm[,(800*i-799):(800*i)]
		
	}	
		colnumber <- ncol(CD45ImmuneCells_tpm_20groups[[i]])
		
		CD45ImmuneCells_tpm_20groups[[i]]$Gene <- rownames(CD45ImmuneCells_tpm_20groups[[i]])

		
		CD45ImmuneCells_tpm_20groups[[i]] <- CD45ImmuneCells_tpm_20groups[[i]][,c((colnumber+1),1:colnumber)]
		
		filename <- paste0("Group",i,"_CD45ImmuneCells_tpm.txt")

		names(CD45ImmuneCells_tpm_20groups)[i] <- gsub("\\.txt","",filename)
		
		print(dim(CD45ImmuneCells_tpm_20groups[[i]]))
		
		print(head(CD45ImmuneCells_tpm_20groups[[i]][,1:4]))
		
		write.table(CD45ImmuneCells_tpm_20groups[[i]], filename, sep='\t', row.names=F, quote=F)

}

names(CD45ImmuneCells_tpm_20groups)



######################################################################

# separate NR (NonResponder) & R (Responder) cells


NR.cells <- as.character(DefineT_PhenoInfo_V2$Cell.Name[DefineT_PhenoInfo_V2$Response=='Non-responder'])

length(NR.cells)

head(NR.cells)

tail(NR.cells)


R.cells <- as.character(DefineT_PhenoInfo_V2$Cell.Name[DefineT_PhenoInfo_V2$Response=='Responder'])

length(R.cells)

head(R.cells)

tail(R.cells)


length(NR.cells)+length(R.cells)


DefineT_PhenoInfo_V2$Response <- factor(DefineT_PhenoInfo_V2$Response,levels=c("Responder","Non-responder"))

DefineT_PhenoInfo_V2 <- DefineT_PhenoInfo_V2[order(DefineT_PhenoInfo_V2$Response,DefineT_PhenoInfo_V2$SampleID,DefineT_PhenoInfo_V2$Therapy),]

dim(DefineT_PhenoInfo_V2)

head(DefineT_PhenoInfo_V2)

tail(DefineT_PhenoInfo_V2)

table(DefineT_PhenoInfo_V2$Cluster.number)

sum(table(DefineT_PhenoInfo_V2$Cluster.number))


DefineT_PhenoInfo_V2_Responder <- DefineT_PhenoInfo_V2[DefineT_PhenoInfo_V2$Cell.Name %in% R.cells,]

dim(DefineT_PhenoInfo_V2_Responder)

head(DefineT_PhenoInfo_V2_Responder)


DefineT_PhenoInfo_V2_NonResponder <- DefineT_PhenoInfo_V2[DefineT_PhenoInfo_V2$Cell.Name %in% NR.cells,]

dim(DefineT_PhenoInfo_V2_NonResponder)

head(DefineT_PhenoInfo_V2_NonResponder)


############################## Follow "Seurat1_Tutorial Integrating stimulated vs. control PBMC datasets to learn cell-type specific responses.pdf" ################################ 

######## split into R and NR groups

CD45ImmuneCells_tpm_R <- CD45ImmuneCells_tpm[,R.cells]

dim(CD45ImmuneCells_tpm_R)

head(CD45ImmuneCells_tpm_R[,1:6])


CD45ImmuneCells_tpm_NR <- CD45ImmuneCells_tpm[,NR.cells]

dim(CD45ImmuneCells_tpm_NR)


# Set up Responder object

Responder <- CreateSeuratObject(counts = CD45ImmuneCells_tpm_R, project = "IMMUNE_Responder", min.cells = 5, meta.data = DefineT_PhenoInfo_V2_Responder) # min.cells=5: 0.1% of the data of 5564 Responder

Responder$Pheno <- "Responder"

Responder <- subset(x = Responder, subset = nFeature_RNA > 200)

Responder <- NormalizeData(object = Responder, verbose = FALSE)

Responder <- FindVariableFeatures(object = Responder, selection.method = "vst", nfeatures = 2000)

str(Responder)

dim(Responder@meta.data)

head(Responder@meta.data)

tail(Responder@meta.data)


# Set up NonResponder object

NonResponder <- CreateSeuratObject(counts = CD45ImmuneCells_tpm_NR, project = "IMMUNE_NonResponder", min.cells = 10, meta.data = DefineT_PhenoInfo_V2_NonResponder)  # min.cells=10: 0.1% of the data of 10727 NonResponder

NonResponder$Pheno <- "NonResponder"

NonResponder <- subset(x = NonResponder, subset = nFeature_RNA > 200)

NonResponder <- NormalizeData(object = NonResponder, verbose = FALSE)

NonResponder <- FindVariableFeatures(object = NonResponder, selection.method = "vst", nfeatures = 2000)

str(NonResponder)

dim(NonResponder@meta.data)

head(NonResponder@meta.data)

tail(NonResponder@meta.data)


immune.anchors <- FindIntegrationAnchors(object.list = list(Responder, NonResponder), dims = 1:20)

immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


DefaultAssay(object = immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering

immune.combined <- ScaleData(object = immune.combined, verbose = FALSE)

immune.combined <- RunPCA(object = immune.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering

immune.combined <- RunUMAP(object = immune.combined, reduction = "pca", dims = 1:20)

immune.combined <- FindNeighbors(object = immune.combined, reduction = "pca", dims = 1:20)

immune.combined <- FindClusters(immune.combined, resolution = 0.5)


immune.combined@meta.data$seurat_clusters_V2 <- immune.combined@meta.data$seurat_clusters

immune.combined@meta.data$seurat_clusters_V2 <- as.character(immune.combined@meta.data$seurat_clusters_V2)

immune.combined@meta.data$seurat_clusters_V2[immune.combined@meta.data$seurat_clusters_V2=='0'] <- '13'

Idents(immune.combined) <- "seurat_clusters_V2"

table(Idents(immune.combined))

str(immune.combined@meta.data)


# Visualization

tiff(file = "Fig1_ImmuneCells_Clustering.tiff", width = 12000, height = 7000, units = "px", res = 600)

p1 <- DimPlot(object = immune.combined, reduction = "umap", group.by = "Pheno")
p2 <- DimPlot(object = immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

dev.off()


tiff(file = "Fig2_ImmuneCells_Clustering_condition_sidebyside.tiff", width = 12000, height = 7000, units = "px", res = 600)

DimPlot(object = immune.combined, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 6)

dev.off()


DefaultAssay(object = immune.combined) <- "RNA"

Cluster6.markers <- FindConservedMarkers(object = immune.combined, ident.1 = 7, grouping.var = "Pheno", verbose = FALSE)

head(x = Cluster6.markers)


table(Idents(immune.combined))

str(immune.combined@meta.data)

head(immune.combined@meta.data)

tail(immune.combined@meta.data)


unique(immune.combined@meta.data$Pheno)

Idents(immune.combined) <- immune.combined@meta.data$Pheno

immune.combined.Res <- subset(immune.combined, idents = "Responder")

immune.combined.NonRes <- subset(immune.combined, idents = "NonResponder")

str(immune.combined.Res@meta.data)

str(immune.combined.NonRes@meta.data)

Idents(immune.combined) <- "seurat_clusters_V2"

Idents(immune.combined.Res) <- "seurat_clusters_V2"

Idents(immune.combined.NonRes) <- "seurat_clusters_V2"


PanglaoDB_markers <- read.table('cellmarkers.short.txt',sep='\t',header=T,fill=T)

dim(PanglaoDB_markers)

head(PanglaoDB_markers)

tail(PanglaoDB_markers)

sort(unique(PanglaoDB_markers$Organ))


length(unique(PanglaoDB_markers$Cell.type))

unique(PanglaoDB_markers$Cell.type)


PanglaoDB_markers_immune <- PanglaoDB_markers[PanglaoDB_markers$Organ == 'Immune system',]

dim(PanglaoDB_markers_immune)

head(PanglaoDB_markers_immune)


length(unique(PanglaoDB_markers_immune$Cell.type))

sort(unique(PanglaoDB_markers_immune$Cell.type))


immune.type <- as.character(unique(PanglaoDB_markers_immune$Cell.type))

immune.type

immune.type.rename <- gsub('\\s+','.',immune.type)

immune.type.rename <- gsub('-','_',immune.type.rename)

immune.type.rename


celltype.list <- list()

for (i in 1:length(immune.type)) {
		
		celltype.list[[i]] <- as.character(PanglaoDB_markers_immune$Official.gene.symbol[PanglaoDB_markers_immune$Cell.type==immune.type[i]])
				
		names(celltype.list)[i] <- immune.type.rename[i]
	
}

str(celltype.list)

class(celltype.list$Monocytes)

length(celltype.list$Monocytes)

celltype.list$Monocytes

length(celltype.list$Dendritic.cells)

celltype.list$Dendritic.cells

celltype.list[20:25]


#tiff(file = "Fp1_T_NK_markers_R.tiff", width = 12000, height = 7000, units = "px", res = 600)

#FeaturePlot(object = immune.combined.Res, features = c("CD3D", "CD3E", "CD3G", "CD8A",
#"CD4", "FOXP3", "MKI67", "NCR1""SELL", "CREM", "CD8A", "GNLY", "CD79A", 
#    "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")

#dev.off()


# cell type specific markers

tiff('Fp1_T_cycling_NK_FeaturePlot_Res.tiff', width = 7000, height = 7000, units = "px", res = 600)

FeaturePlot(immune.combined.Res, features = c("CD3D", "CD3E", "CD3G","CD8A","CD4","FOXP3","MKI67","NCR1","NCAM1"), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('Fp2_T_cycling_NK_FeaturePlot_NonRes.tiff', width = 7000, height = 7000, units = "px", res = 600)

FeaturePlot(immune.combined.NonRes, features = c("CD3D", "CD3E", "CD3G","CD8A","CD4","FOXP3","MKI67","NCR1","NCAM1"), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('Fp3_B_Plasma_MM_DC_FeaturePlot_Res.tiff', width = 7000, height = 7000, units = "px", res = 600)

FeaturePlot(immune.combined.Res, features = c("CD19","MZB1","MARCO","MERTK","FCER1A","NAPSA"), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('Fp4_B_Plasma_MM_DC_FeaturePlot_NonRes.tiff', width = 7000, height = 7000, units = "px", res = 600)

FeaturePlot(immune.combined.NonRes, features = c("CD19","MZB1","MARCO","MERTK","FCER1A","NAPSA"), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Res

tiff('FeaturePlot_Monocytes_Res.tiff', width = 16000, height = 18000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = celltype.list$Monocytes, min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()




tiff('FeaturePlot_DC_pt1_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = celltype.list$Dendritic.cells[1:35], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt2_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = celltype.list$Dendritic.cells[36:70], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt3_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = celltype.list$Dendritic.cells[71:105], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt4_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = celltype.list$Dendritic.cells[106:133], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# NonRes

tiff('FeaturePlot_Monocytes_NonRes.tiff', width = 16000, height = 18000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = celltype.list$Monocytes, min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()




tiff('FeaturePlot_DC_pt1_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = celltype.list$Dendritic.cells[1:35], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt2_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = celltype.list$Dendritic.cells[36:70], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt3_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = celltype.list$Dendritic.cells[71:105], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_DC_pt4_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = celltype.list$Dendritic.cells[106:133], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Replot clustering of R and NR based on canonical markers determined cell types

sort(levels(immune.combined.Res))

levels(immune.combined.NonRes)


new.cluster.ids <- c("CD8+ T cells", "CD8+ T cells", "B cells", "CD8+ T cells", "Macrophages/Monocytes", "Tgd cells", "NK cells", 
"Macrophages/Monocytes", "CD8+ T cells", "Plasma cells", "Dendritic cells", "B cells", "CD4+ T cells")

names(new.cluster.ids) <- as.character(1:13)

new.cluster.ids



# immune.combined.Res

immune.combined.Res <- RenameIdents(immune.combined.Res, new.cluster.ids)

tiff(file = "Fig3_ImmuneCells_Res_Clustering_MainTypes.tiff", width = 5000, height = 5000, units = "px", res = 600)

DimPlot(immune.combined.Res, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()

dev.off()


#immune.combined.NonRes

immune.combined.NonRes <- RenameIdents(immune.combined.NonRes, new.cluster.ids)

tiff(file = "Fig4_ImmuneCells_NonRes_Clustering_MainTypes.tiff", width = 5000, height = 5000, units = "px", res = 600)

DimPlot(immune.combined.NonRes, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()

dev.off()


saveRDS(immune.combined, file="immune.combined.Jan2020.rds")

saveRDS(immune.combined.Res, file="immune.combined.Res.Jan2020.rds")

saveRDS(immune.combined.NonRes, file="immune.combined.NonRes.Jan2020.rds")


################################################################# More clusters here from 717 ##############################################################

# Want much more clusters

str(immune.combined@meta.data)


# Run the standard workflow for visualization and clustering

immune.combined.more.clusters <- immune.combined

immune.combined.more.clusters <- ScaleData(object = immune.combined.more.clusters, verbose = FALSE)

immune.combined.more.clusters <- RunPCA(object = immune.combined.more.clusters, npcs = 30, verbose = FALSE)


# t-SNE and Clustering

immune.combined.more.clusters <- RunUMAP(object = immune.combined.more.clusters, reduction = "pca", dims = 1:20)

immune.combined.more.clusters <- FindNeighbors(object = immune.combined.more.clusters, reduction = "pca", dims = 1:20)

immune.combined.more.clusters <- FindClusters(immune.combined.more.clusters, resolution = 1)

str(immune.combined.more.clusters@meta.data)


# Visualization

new.ids <- as.character(1:23)

names(new.ids) <- as.character(0:22)

new.ids

immune.combined.more.clusters <- RenameIdents(immune.combined.more.clusters, new.ids)

table(Idents(immune.combined.more.clusters))

immune.combined.more.clusters[['RNA_snn_res.1.V2']] <- Idents(immune.combined.more.clusters)

colnames(immune.combined.more.clusters@meta.data)[7] <- 'integrated_snn_res.0.5.V2'

str(immune.combined.more.clusters@meta.data)

head(rownames(immune.combined.more.clusters@meta.data))



tiff(file = "FigM1_ImmuneCells_Clustering.tiff", width = 12000, height = 7000, units = "px", res = 600)

p1 <- DimPlot(object = immune.combined.more.clusters, reduction = "umap", group.by = "Pheno")
p2 <- DimPlot(object = immune.combined.more.clusters, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

dev.off()


tiff(file = "FigM2_ImmuneCells_Clustering_condition_sidebyside.tiff", width = 12000, height = 9000, units = "px", res = 600)

DimPlot(object = immune.combined.more.clusters, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 5)

dev.off()


tiff(file = "FigM2_V2_ImmuneCells_Clustering_condition_sidebyside.tiff", width = 8500, height = 5000, units = "px", res = 1000)

DimPlot(object = immune.combined.more.clusters, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 4.5)

dev.off()



new.ids2 <- c("CD8+ T cells", "Regulatory T cells", "CD4+ T cells", "CD8+ T cells", "CD8+ T cells", "Macrophages/Monocytes", 
"CD8+ T cells", "Tgd cells", "MKI67hi Lymph.", "CD8+ T cells", "CD8+ T cells", "Macrophages/Monocytes", "B cells", "B cells", 
"NK cells", "MKI67hi Lymph.", "B cells", "Plasma cells", "Dendritic cells", "CD8+ T cells", "Tgd cells", "B cells", "Macrophages/Monocytes")

names(new.ids2) <- as.character(1:23)

new.ids2




immune.combined.more.clusters <- RenameIdents(immune.combined.more.clusters, new.ids2)

immune.combined.more.clusters[['Cell.Type']] <- Idents(immune.combined.more.clusters)

table(Idents(immune.combined.more.clusters))

str(immune.combined.more.clusters@meta.data)


Idents(immune.combined.more.clusters) <- immune.combined.more.clusters@meta.data$Pheno

immune.combined.more.clusters.Res <- subset(immune.combined.more.clusters, idents = "Responder")

immune.combined.more.clusters.NonRes <- subset(immune.combined.more.clusters, idents = "NonResponder")

str(immune.combined.more.clusters.Res@meta.data)

str(immune.combined.more.clusters.NonRes@meta.data)

Idents(immune.combined.more.clusters) <- "Cell.Type"

Idents(immune.combined.more.clusters.Res) <- "Cell.Type"

Idents(immune.combined.more.clusters.NonRes) <- "Cell.Type"


immune.combined.more.clusters@meta.data$Cell.Type.V2 <- as.character(immune.combined.more.clusters@meta.data$Cell.Type)


immune.combined.more.clusters@meta.data$Cell.Type.V2[immune.combined.more.clusters@meta.data$Cell.Type.V2=="Macrophages/Monocytes"] <- "Macrophages"

immune.combined.more.clusters@meta.data$Cell.Type.V2[immune.combined.more.clusters@meta.data$Cell.Type.V2=="MKI67hi Lymph."] <- "MKI67hi"


table(immune.combined.more.clusters@meta.data$Cell.Type.V2)



Idents(immune.combined.more.clusters) <- "Cell.Type.V2"

table(Idents(immune.combined.more.clusters))


tiff(file = "PaperFig1C_ImmuneCells_Clustering_Overall.tiff", width = 5500, height = 5500, units = "px", res = 1000)

DimPlot(object = immune.combined.more.clusters, reduction = "umap", label = TRUE, label.size = 4) + NoLegend()

dev.off()


CellMarkers <- c('CD3D','CD8A','CD4','FOXP3','MKI67','NCR1','NCAM1','CD19','MZB1','MARCO','MERTK','FCER1A') # omit 'CD3E','CD3G','NAPSA'

tiff('PaperFig1B_FeaturePlot_immune.combined.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.more.clusters, features = CellMarkers, min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()



tiff('PaperFig1B_pt1_6genes_FeaturePlot_immune.combined.tiff', width = 6000, height = 6000, units = "px", res = 1000)

FeaturePlot(immune.combined.more.clusters, features = CellMarkers[c(1,2,5,6,9,10)], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "darkblue"))

dev.off()


tiff('PaperFig1B_pt2_6genes_FeaturePlot_immune.combined.tiff', width = 6000, height = 6000, units = "px", res = 1000)

FeaturePlot(immune.combined.more.clusters, features = CellMarkers[c(3,4,7,8,11,12)], min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "darkblue"))

dev.off()


# immune.combined.more.clusters.Res

tiff(file = "FigM3_ImmuneCells_Res_MoreCellClusters.tiff", width = 7000, height = 7000, units = "px", res = 600)

DimPlot(immune.combined.more.clusters.Res, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()

dev.off()


#immune.combined.more.clusters.NonRes

tiff(file = "FigM4_ImmuneCells_NonRes_MoreCellClusters.tiff", width = 7000, height = 7000, units = "px", res = 600)

DimPlot(immune.combined.more.clusters.NonRes, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()

dev.off()


############# Find marker genes for each of the 23 clusters ##############

Idents(immune.combined.more.clusters) <- "RNA_snn_res.1.V2"

table(Idents(immune.combined.more.clusters))


# find markers for every cluster compared to all remaining cells, report only the positive ones

immune.combined.more.clusters.markers <- FindAllMarkers(immune.combined.more.clusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

immune.combined.more.clusters.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


top20 <- immune.combined.more.clusters.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

tiff("immune23clusters.heatmap.tiff", width = 15000, height = 15000, units = "px", res = 600)
  
DoHeatmap(immune.combined.more.clusters, features = top20$gene, hjust = 0, angle = 90) + NoLegend()

dev.off()

table(immune.combined.more.clusters@meta.data$Pheno)


########### calculate proportion of each cluster in R and NR


###

Idents(immune.combined.more.clusters.Res) <- 'RNA_snn_res.1.V2'

table(Idents(immune.combined.more.clusters.Res))

sum(table(Idents(immune.combined.more.clusters.Res)))

table(Idents(immune.combined.more.clusters.Res))/sum(table(Idents(immune.combined.more.clusters.Res)))

table(Idents(immune.combined.more.clusters.Res))/nrow(immune.combined@meta.data)

str(immune.combined.more.clusters.Res@meta.data)

immune.combined.more.clusters.Res@meta.data$SampleID_Drug <- paste0(immune.combined.more.clusters.Res@meta.data$SampleID, '_', immune.combined.more.clusters.Res@meta.data$integrated_snn_res.0.5.V2)

table(immune.combined.more.clusters.Res@meta.data$SampleID_Drug)

R.stratified.group.names <- names(table(immune.combined.more.clusters.Res@meta.data$SampleID_Drug))

R.stratified.group.names <- gsub('_.*_anti','',R.stratified.group.names)

R.stratified.group.names

table(R.stratified.group.names)

immune.combined.more.clusters.Res@meta.data$Timing_Drug <- paste0(immune.combined.more.clusters.Res@meta.data$TreatTime, '_', immune.combined.more.clusters.Res@meta.data$integrated_snn_res.0.5.V2)

table(immune.combined.more.clusters.Res@meta.data$Timing_Drug)


df.table1 <- table(Idents(immune.combined.more.clusters.Res))/sum(table(Idents(immune.combined.more.clusters.Res)))

Res.proportions.AllClusters <- as.data.frame(df.table1)

Res.proportions.AllClusters$Phenotype <- "Responder/Regression"

colnames(Res.proportions.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(Res.proportions.AllClusters)

head(Res.proportions.AllClusters)

tail(Res.proportions.AllClusters)


###

Idents(immune.combined.more.clusters.NonRes) <- 'RNA_snn_res.1.V2'

table(Idents(immune.combined.more.clusters.NonRes))

sum(table(Idents(immune.combined.more.clusters.NonRes)))

table(Idents(immune.combined.more.clusters.NonRes))/sum(table(Idents(immune.combined.more.clusters.NonRes)))

table(Idents(immune.combined.more.clusters.NonRes))/nrow(immune.combined@meta.data)


str(immune.combined.more.clusters.NonRes@meta.data)

immune.combined.more.clusters.NonRes@meta.data$SampleID_Drug <- paste0(immune.combined.more.clusters.NonRes@meta.data$SampleID, '_', immune.combined.more.clusters.NonRes@meta.data$integrated_snn_res.0.5.V2)

table(immune.combined.more.clusters.NonRes@meta.data$SampleID_Drug)

NR.stratified.group.names <- names(table(immune.combined.more.clusters.NonRes@meta.data$SampleID_Drug))

NR.stratified.group.names <- gsub('_.*_anti','',NR.stratified.group.names)

NR.stratified.group.names

table(NR.stratified.group.names)

immune.combined.more.clusters.NonRes@meta.data$Timing_Drug <- paste0(immune.combined.more.clusters.NonRes@meta.data$TreatTime, '_', immune.combined.more.clusters.NonRes@meta.data$integrated_snn_res.0.5.V2)

table(immune.combined.more.clusters.NonRes@meta.data$Timing_Drug)


df.table2 <- table(Idents(immune.combined.more.clusters.NonRes))/sum(table(Idents(immune.combined.more.clusters.NonRes)))

NonRes.proportions.AllClusters <- as.data.frame(df.table2)

NonRes.proportions.AllClusters$Phenotype <- "Non-responder/Progression"

colnames(NonRes.proportions.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(NonRes.proportions.AllClusters)

head(NonRes.proportions.AllClusters)

tail(NonRes.proportions.AllClusters)



# Proportion of all clusters

RvsNR.proportions.AllClusters <- rbind(Res.proportions.AllClusters, NonRes.proportions.AllClusters)

dim(RvsNR.proportions.AllClusters)

head(RvsNR.proportions.AllClusters)

tail(RvsNR.proportions.AllClusters)


tiff('Barplot.RvsNR.proportions.AllClusters.tiff', width = 6000, height = 6000, units = "px", res = 300)

ggplot(RvsNR.proportions.AllClusters,aes(x=factor(Phenotype), y = Percentage)) + 
    facet_wrap(~Cluster) + theme(strip.text.x = element_text(size = 18, colour = "red", angle = 0)) +
    geom_bar(aes(fill = factor(Phenotype)),stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

write.table(RvsNR.proportions.AllClusters, 'RvsNR.proportions.AllClusters.xls', sep='\t', row.names=F, quote=F)



##### NR vs R: before-anti-PD-1

table(immune.combined.more.clusters.Res@meta.data$Timing_Drug)

table(immune.combined.more.clusters.NonRes@meta.data$Timing_Drug)


Idents(immune.combined.more.clusters.Res) <- 'Timing_Drug'

Idents(immune.combined.more.clusters.NonRes) <- 'Timing_Drug'


table(Idents(immune.combined.more.clusters.Res))

table(Idents(immune.combined.more.clusters.NonRes))



# R_before_PD1

R_before_PD1 <- subset(immune.combined.more.clusters.Res, idents='Pre_anti-PD1')

str(R_before_PD1@meta.data)

Idents(R_before_PD1) <- 'RNA_snn_res.1.V2'

table(Idents(R_before_PD1))


R_before_PD1.table <- table(Idents(R_before_PD1))/sum(table(Idents(R_before_PD1)))

R_before_PD1.fraction.AllClusters <- as.data.frame(R_before_PD1.table)

R_before_PD1.fraction.AllClusters$Phenotype <- "Responder/Regression"

colnames(R_before_PD1.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(R_before_PD1.fraction.AllClusters)

head(R_before_PD1.fraction.AllClusters)

tail(R_before_PD1.fraction.AllClusters)


# NR_before_PD1

NR_before_PD1 <- subset(immune.combined.more.clusters.NonRes, idents='Pre_anti-PD1')

str(NR_before_PD1@meta.data)

Idents(NR_before_PD1) <- 'RNA_snn_res.1.V2'

table(Idents(NR_before_PD1))


NR_before_PD1.table <- table(Idents(NR_before_PD1))/sum(table(Idents(NR_before_PD1)))

NR_before_PD1.fraction.AllClusters <- as.data.frame(NR_before_PD1.table)

NR_before_PD1.fraction.AllClusters$Phenotype <- "Non-responder/Progression"

colnames(NR_before_PD1.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(NR_before_PD1.fraction.AllClusters)

head(NR_before_PD1.fraction.AllClusters)

tail(NR_before_PD1.fraction.AllClusters)



# fractions of all clusters

NRvsR.before_PD1.fractions.AllClusters <- rbind(R_before_PD1.fraction.AllClusters, NR_before_PD1.fraction.AllClusters)

NRvsR.before_PD1.fractions.AllClusters <- rbind(NRvsR.before_PD1.fractions.AllClusters, data.frame(Cluster='21',Percentage=0,Phenotype='Responder/Regression'))

NRvsR.before_PD1.fractions.AllClusters <- rbind(NRvsR.before_PD1.fractions.AllClusters, data.frame(Cluster='23',Percentage=0,Phenotype='Non-responder/Progression'))

NRvsR.before_PD1.fractions.AllClusters$Phenotype[NRvsR.before_PD1.fractions.AllClusters$Phenotype=='Responder/Regression'] <- 'R'

NRvsR.before_PD1.fractions.AllClusters$Phenotype[NRvsR.before_PD1.fractions.AllClusters$Phenotype=='Non-responder/Progression'] <- 'NR'

NRvsR.before_PD1.fractions.AllClusters$Phenotype <- factor(NRvsR.before_PD1.fractions.AllClusters$Phenotype, levels=c('R','NR'))

NRvsR.before_PD1.fractions.AllClusters$Cluster <- as.integer(NRvsR.before_PD1.fractions.AllClusters$Cluster)

NRvsR.before_PD1.fractions.AllClusters <- NRvsR.before_PD1.fractions.AllClusters[order(NRvsR.before_PD1.fractions.AllClusters$Phenotype, NRvsR.before_PD1.fractions.AllClusters$Cluster),]

NRvsR.before_PD1.fractions.AllClusters[35,2] <- 0.052273425

NRvsR.before_PD1.fractions.AllClusters[21,2] <- 0.000768049

NRvsR.before_PD1.fractions.AllClusters[44,2] <- 0.058774139

NRvsR.before_PD1.fractions.AllClusters[45,2] <- 0.002039102


NRvsR.before_PD1.fractions.AllClusters


tiff('Barplot.NRvsR.before_PD1.fractions.AllClusters.tiff', width = 4000, height = 6000, units = "px", res = 300)

ggplot(NRvsR.before_PD1.fractions.AllClusters,aes(x=factor(Phenotype), y = Percentage)) + 
    facet_wrap(~Cluster) + theme(strip.text.x = element_text(size = 18, colour = "red", angle = 0)) +
    geom_bar(aes(fill = factor(Phenotype)),stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

write.table(NRvsR.before_PD1.fractions.AllClusters, 'NRvsR.before_PD1.fractions.AllClusters.xls', sep='\t', row.names=F, quote=F)



# plot fold changes of different clusters between NR and R groups


data <-read.table("BeforePD1.fractions.AllClusters.FoldChange.NRvsR.txt",sep='\t',header=T) 

tiff('BeforePD1.fractions.AllClusters.FoldChange.NRvsR.tiff',width = 6000, height = 4000, units = "px", res = 600)

bp <- barplot(data$NRvsR_FoldChange, main = "", xlab = "Single-cell Cluster", 
ylab = "Fold change of percentage of CD45+ immune cells (NR/R)", ylim = c(-10,16),names.arg = data$Cluster, axes=F,
    col = rainbow(nrow(data)))

xval = seq(-10,16,by=2)

axis(side = 2, at = xval, labels = FALSE, xpd=T)
axis(side = 2, at = xval, tick = FALSE, labels = xval, xpd=T)

abline(h=0)
abline(h=6,lty=2)
abline(h=-6,lty=2)

dev.off()




##### NR vs R: after-anti-PD-1

table(immune.combined.more.clusters.Res@meta.data$Timing_Drug)

table(immune.combined.more.clusters.NonRes@meta.data$Timing_Drug)


Idents(immune.combined.more.clusters.Res) <- 'Timing_Drug'

Idents(immune.combined.more.clusters.NonRes) <- 'Timing_Drug'


table(Idents(immune.combined.more.clusters.Res))

table(Idents(immune.combined.more.clusters.NonRes))



# R_after_PD1

R_after_PD1 <- subset(immune.combined.more.clusters.Res, idents='Post_anti-PD1')

str(R_after_PD1@meta.data)

Idents(R_after_PD1) <- 'RNA_snn_res.1.V2'

table(Idents(R_after_PD1))


R_after_PD1.table <- table(Idents(R_after_PD1))/sum(table(Idents(R_after_PD1)))

R_after_PD1.fraction.AllClusters <- as.data.frame(R_after_PD1.table)

R_after_PD1.fraction.AllClusters$Phenotype <- "Responder/Regression"

colnames(R_after_PD1.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(R_after_PD1.fraction.AllClusters)

head(R_after_PD1.fraction.AllClusters)

tail(R_after_PD1.fraction.AllClusters)


# NR_after_PD1

NR_after_PD1 <- subset(immune.combined.more.clusters.NonRes, idents='Post_anti-PD1')

str(NR_after_PD1@meta.data)

Idents(NR_after_PD1) <- 'RNA_snn_res.1.V2'

table(Idents(NR_after_PD1))


NR_after_PD1.table <- table(Idents(NR_after_PD1))/sum(table(Idents(NR_after_PD1)))

NR_after_PD1.fraction.AllClusters <- as.data.frame(NR_after_PD1.table)

NR_after_PD1.fraction.AllClusters$Phenotype <- "Non-responder/Progression"

colnames(NR_after_PD1.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(NR_after_PD1.fraction.AllClusters)

head(NR_after_PD1.fraction.AllClusters)

tail(NR_after_PD1.fraction.AllClusters)



# fractions of all clusters

NRvsR.after_PD1.fractions.AllClusters <- rbind(R_after_PD1.fraction.AllClusters, NR_after_PD1.fraction.AllClusters)

#NRvsR.after_PD1.fractions.AllClusters <- rbind(NRvsR.after_PD1.fractions.AllClusters, data.frame(Cluster='21',Percentage=0,Phenotype='Responder/Regression'))

#NRvsR.after_PD1.fractions.AllClusters <- rbind(NRvsR.after_PD1.fractions.AllClusters, data.frame(Cluster='23',Percentage=0,Phenotype='Non-responder/Progression'))

NRvsR.after_PD1.fractions.AllClusters$Phenotype[NRvsR.after_PD1.fractions.AllClusters$Phenotype=='Responder/Regression'] <- 'R'

NRvsR.after_PD1.fractions.AllClusters$Phenotype[NRvsR.after_PD1.fractions.AllClusters$Phenotype=='Non-responder/Progression'] <- 'NR'

NRvsR.after_PD1.fractions.AllClusters$Phenotype <- factor(NRvsR.after_PD1.fractions.AllClusters$Phenotype, levels=c('R','NR'))

NRvsR.after_PD1.fractions.AllClusters$Cluster <- as.integer(NRvsR.after_PD1.fractions.AllClusters$Cluster)

NRvsR.after_PD1.fractions.AllClusters <- NRvsR.after_PD1.fractions.AllClusters[order(NRvsR.after_PD1.fractions.AllClusters$Phenotype, NRvsR.after_PD1.fractions.AllClusters$Cluster),]

NRvsR.after_PD1.fractions.AllClusters[35,2] <- 0.052273425

#NRvsR.after_PD1.fractions.AllClusters[21,2] <- 0.000768049

#NRvsR.after_PD1.fractions.AllClusters[44,2] <- 0.058774139

NRvsR.after_PD1.fractions.AllClusters[22,2] <- 0.01268504


NRvsR.after_PD1.fractions.AllClusters


tiff('Barplot.NRvsR.after_PD1.fractions.AllClusters.tiff', width = 4000, height = 6000, units = "px", res = 300)

ggplot(NRvsR.after_PD1.fractions.AllClusters,aes(x=factor(Phenotype), y = Percentage)) + 
    facet_wrap(~Cluster) + theme(strip.text.x = element_text(size = 18, colour = "red", angle = 0)) +
    geom_bar(aes(fill = factor(Phenotype)),stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

write.table(NRvsR.after_PD1.fractions.AllClusters, 'NRvsR.after_PD1.fractions.AllClusters.xls', sep='\t', row.names=F, quote=F)



# plot fold changes of different clusters between NR and R groups


data <-read.table("AfterPD1.fractions.AllClusters.FoldChange.NRvsR.txt",sep='\t',header=T) 

tiff('AfterPD1.fractions.AllClusters.FoldChange.NRvsR.tiff',width = 6000, height = 4000, units = "px", res = 600)

bp <- barplot(data$NRvsR_FoldChange, main = "", xlab = "Single-cell Cluster", 
ylab = "Fold change of percentage of CD45+ immune cells (NR/R)", ylim = c(-10,16),names.arg = data$Cluster, axes=F,
    col = rainbow(nrow(data)))

xval = seq(-10,16,by=2)

axis(side = 2, at = xval, labels = FALSE, xpd=T)
axis(side = 2, at = xval, tick = FALSE, labels = xval, xpd=T)

abline(h=0)
abline(h=6,lty=2)
abline(h=-6,lty=2)

dev.off()




##### NR vs R: after-anti-both

table(immune.combined.more.clusters.Res@meta.data$Timing_Drug)

table(immune.combined.more.clusters.NonRes@meta.data$Timing_Drug)


Idents(immune.combined.more.clusters.Res) <- 'Timing_Drug'

Idents(immune.combined.more.clusters.NonRes) <- 'Timing_Drug'


table(Idents(immune.combined.more.clusters.Res))

table(Idents(immune.combined.more.clusters.NonRes))



# R_after_both

R_after_both <- subset(immune.combined.more.clusters.Res, idents='Post_anti-CTLA4+PD1')

str(R_after_both@meta.data)

Idents(R_after_both) <- 'RNA_snn_res.1.V2'

table(Idents(R_after_both))


R_after_both.table <- table(Idents(R_after_both))/sum(table(Idents(R_after_both)))

R_after_both.fraction.AllClusters <- as.data.frame(R_after_both.table)

R_after_both.fraction.AllClusters$Phenotype <- "Responder/Regression"

colnames(R_after_both.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(R_after_both.fraction.AllClusters)

head(R_after_both.fraction.AllClusters)

tail(R_after_both.fraction.AllClusters)


# NR_after_both

NR_after_both <- subset(immune.combined.more.clusters.NonRes, idents='Post_anti-CTLA4+PD1')

str(NR_after_both@meta.data)

Idents(NR_after_both) <- 'RNA_snn_res.1.V2'

table(Idents(NR_after_both))


NR_after_both.table <- table(Idents(NR_after_both))/sum(table(Idents(NR_after_both)))

NR_after_both.fraction.AllClusters <- as.data.frame(NR_after_both.table)

NR_after_both.fraction.AllClusters$Phenotype <- "Non-responder/Progression"

colnames(NR_after_both.fraction.AllClusters)[1:2] <- c("Cluster","Percentage")

dim(NR_after_both.fraction.AllClusters)

head(NR_after_both.fraction.AllClusters)

tail(NR_after_both.fraction.AllClusters)



# fractions of all clusters

NRvsR.after_both.fractions.AllClusters <- rbind(R_after_both.fraction.AllClusters, NR_after_both.fraction.AllClusters)

NRvsR.after_both.fractions.AllClusters <- rbind(NRvsR.after_both.fractions.AllClusters, data.frame(Cluster='23',Percentage=0,Phenotype='Responder/Regression'))

NRvsR.after_both.fractions.AllClusters <- rbind(NRvsR.after_both.fractions.AllClusters, data.frame(Cluster='4',Percentage=0,Phenotype='Non-responder/Progression'))

NRvsR.after_both.fractions.AllClusters <- rbind(NRvsR.after_both.fractions.AllClusters, data.frame(Cluster='23',Percentage=0,Phenotype='Non-responder/Progression'))


NRvsR.after_both.fractions.AllClusters$Phenotype[NRvsR.after_both.fractions.AllClusters$Phenotype=='Responder/Regression'] <- 'R'

NRvsR.after_both.fractions.AllClusters$Phenotype[NRvsR.after_both.fractions.AllClusters$Phenotype=='Non-responder/Progression'] <- 'NR'

NRvsR.after_both.fractions.AllClusters$Phenotype <- factor(NRvsR.after_both.fractions.AllClusters$Phenotype, levels=c('R','NR'))

NRvsR.after_both.fractions.AllClusters$Cluster <- as.integer(NRvsR.after_both.fractions.AllClusters$Cluster)

NRvsR.after_both.fractions.AllClusters <- NRvsR.after_both.fractions.AllClusters[order(NRvsR.after_both.fractions.AllClusters$Phenotype, NRvsR.after_both.fractions.AllClusters$Cluster),]

NRvsR.after_both.fractions.AllClusters[35,2] <- 0.01907563025

#NRvsR.after_both.fractions.AllClusters[21,2] <- 0.000768049

NRvsR.after_both.fractions.AllClusters[44,2] <- 0.008403361

NRvsR.after_both.fractions.AllClusters[22,2] <- 0.011041825


NRvsR.after_both.fractions.AllClusters


tiff('Barplot.NRvsR.after_both.fractions.AllClusters.tiff', width = 4000, height = 6000, units = "px", res = 300)

ggplot(NRvsR.after_both.fractions.AllClusters,aes(x=factor(Phenotype), y = Percentage)) + 
    facet_wrap(~Cluster) + theme(strip.text.x = element_text(size = 18, colour = "red", angle = 0)) +
    geom_bar(aes(fill = factor(Phenotype)),stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()

write.table(NRvsR.after_both.fractions.AllClusters, 'NRvsR.after_both.fractions.AllClusters.xls', sep='\t', row.names=F, quote=F)



# plot fold changes of different clusters between NR and R groups


data <-read.table("AfterBoth.fractions.AllClusters.FoldChange.NRvsR.txt",sep='\t',header=T) 

tiff('AfterBoth.fractions.AllClusters.FoldChange.NRvsR.tiff',width = 6000, height = 4000, units = "px", res = 600)

bp <- barplot(data$NRvsR_FoldChange, main = "", xlab = "Single-cell Cluster", 
ylab = "Fold change of percentage of CD45+ immune cells (NR/R)", ylim = c(-10,16),names.arg = data$Cluster, axes=F,
    col = rainbow(nrow(data)))

xval = seq(-10,16,by=2)

axis(side = 2, at = xval, labels = FALSE, xpd=T)
axis(side = 2, at = xval, tick = FALSE, labels = xval, xpd=T)

abline(h=0)
abline(h=6,lty=2)
abline(h=-6,lty=2)

dev.off()





# Macrophage analyses

table(Idents(immune.combined.more.clusters))

sum(table(Idents(immune.combined.more.clusters)))

MM_cluster6_12_23 <- subset(immune.combined.more.clusters, idents = c("6","12","23"))

str(MM_cluster6_12_23@meta.data)

dim(MM_cluster6_12_23@meta.data)

table(Idents(MM_cluster6_12_23))

sort(unique(MM_cluster6_12_23@meta.data$Cell.Type))



# find markers for every cluster compared to all remaining cells, report only the positive ones in MM population - 3 clusters of 6, 12, 23

MM_cluster6_12_23.markers <- FindAllMarkers(MM_cluster6_12_23, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

MM_cluster6_12_23.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


top20 <- MM_cluster6_12_23.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

tiff("MM_3clusters.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = top20$gene, hjust = 0.5, angle = 0) + NoLegend()

dev.off()

dim(MM_cluster6_12_23.markers)

levels(MM_cluster6_12_23)


### DE analysis of each of the 3 MM clusters

# MM.cluster6 DE analysis

MM.cluster6.de.markers <- FindMarkers(MM_cluster6_12_23, ident.1 = "6", ident.2 = NULL, only.pos = TRUE)

dim(MM.cluster6.de.markers)

MM.cluster6.de.markers

MM.cluster6.de.markers.V2 <- MM.cluster6.de.markers

MM.cluster6.de.markers.V2$FoldChange <- MM.cluster6.de.markers.V2$pct.1/MM.cluster6.de.markers.V2$pct.2

MM.cluster6.de.markers.V2 <- MM.cluster6.de.markers.V2[order(MM.cluster6.de.markers.V2$FoldChange, decreasing=T),]

range(MM.cluster6.de.markers.V2$p_val)

MM.cluster6.de.markers.V2

MM.cluster6.de.markers.V2.topFCgenes <- rownames(MM.cluster6.de.markers.V2[MM.cluster6.de.markers.V2$FoldChange > 1.5,])

MM.cluster6.de.markers.V2.topFCgenes <- MM.cluster6.de.markers.V2.topFCgenes[!MM.cluster6.de.markers.V2.topFCgenes %in% c('ISG20','BIRC3','CLEC10A','FCN1','CD52','CXCL10')]

MM.cluster6.de.markers.V2.topFCgenes



MM.cluster6.DE <- FindMarkers(object = MM_cluster6_12_23, ident.1 = "6", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

dim(MM.cluster6.DE)

head(MM.cluster6.DE)

tail(MM.cluster6.DE)


MM.cluster6.de.markers <- MM.cluster6.de.markers[order(MM.cluster6.de.markers$p_val),]

head(MM.cluster6.de.markers)


# MM.cluster12 DE analysis

MM.cluster12.de.markers <- FindMarkers(MM_cluster6_12_23, ident.1 = "12", ident.2 = NULL, only.pos = TRUE)

dim(MM.cluster12.de.markers)

MM.cluster12.DE <- FindMarkers(object = MM_cluster6_12_23, ident.1 = "12", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

dim(MM.cluster12.DE)

head(MM.cluster12.DE)

tail(MM.cluster12.DE)


MM.cluster12.de.markers <- MM.cluster12.de.markers[order(MM.cluster12.de.markers$p_val),]

head(MM.cluster12.de.markers)



# MM.cluster23 DE analysis

MM.cluster23.de.markers <- FindMarkers(MM_cluster6_12_23, ident.1 = "23", ident.2 = NULL, only.pos = TRUE)

MM.cluster23.de.markers <- MM.cluster23.de.markers[MM.cluster23.de.markers$p_val < 0.05,]

dim(MM.cluster23.de.markers)

range(MM.cluster23.de.markers$p_val)


MM.cluster23.de.markers.V2 <- MM.cluster23.de.markers

MM.cluster23.de.markers.V2$FoldChange <- MM.cluster23.de.markers.V2$pct.1/MM.cluster23.de.markers.V2$pct.2

MM.cluster23.de.markers.V2 <- MM.cluster23.de.markers.V2[order(MM.cluster23.de.markers.V2$FoldChange, decreasing=T),]

range(MM.cluster23.de.markers.V2$p_val)

MM.cluster23.de.markers.V2

MM.cluster23.de.markers.V2.topFCgenes <- rownames(MM.cluster23.de.markers.V2[MM.cluster23.de.markers.V2$FoldChange > 2 & MM.cluster23.de.markers.V2$pct.1 > 0.7,])

MM.cluster23.de.markers.V2.topFCgenes <- MM.cluster23.de.markers.V2.topFCgenes[!MM.cluster23.de.markers.V2.topFCgenes %in% c('ISG20','BIRC3','CLEC10A','FCN1','CD52','CXCL10')]

MM.cluster23.de.markers.V2.topFCgenes




MM.cluster23.DE <- FindMarkers(object = MM_cluster6_12_23, ident.1 = "23", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

dim(MM.cluster23.DE)

head(MM.cluster23.DE)

tail(MM.cluster23.DE)




###### MM_cluster12 NR vs R analysis

MM_cluster12 <- subset(MM_cluster6_12_23, idents='12')

table(Idents(MM_cluster12))

str(MM_cluster12@meta.data)

Idents(MM_cluster12) <- MM_cluster12@meta.data$Pheno
		
table(Idents(MM_cluster12))

tiff(file = "MM_cluster12_Heatmap_NRvR_deGenes.tiff", width = 3000, height = 2500, units = "px", res = 600)

DoHeatmap(object = MM_cluster12, features = rownames(MM.cluster12.de.markers), slot = "scale.data", hjust = 0, angle=0, size=4)

dev.off()


MM_cluster12_Responder_ID <- rownames(MM_cluster12@meta.data[MM_cluster12@meta.data$Pheno=='Responder',])

length(MM_cluster12_Responder_ID)

MM_cluster12_Responder_ID


MM_cluster12_NonResponder_ID <- rownames(MM_cluster12@meta.data[MM_cluster12@meta.data$Pheno=='NonResponder',])

length(MM_cluster12_NonResponder_ID)

MM_cluster12_NonResponder_ID



MM.cluster12.NRvsR.DE <- FindMarkers(object = MM_cluster12, ident.1 = "Responder", ident.2 = "NonResponder", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster12.NRvsR.DE <- MM.cluster12.NRvsR.DE[order(MM.cluster12.NRvsR.DE$avg_logFC,decreasing=T),]

dim(MM.cluster12.NRvsR.DE)

head(MM.cluster12.NRvsR.DE)

tail(MM.cluster12.NRvsR.DE)

MM.cluster12.NRvsR.DE[rownames(MM.cluster12.NRvsR.DE)=='TREM2',]


MM.cluster12.NRvsR.DE.sig <- MM.cluster12.NRvsR.DE[MM.cluster12.NRvsR.DE$p_val_adj < 0.05,]

dim(MM.cluster12.NRvsR.DE.sig)

head(MM.cluster12.NRvsR.DE.sig)

tail(MM.cluster12.NRvsR.DE.sig)


tiff(file = "MM_cluster12_Heatmap_NRvR_deGenes_DE.sig.tiff", width = 3000, height = 2500, units = "px", res = 600)

DoHeatmap(object = MM_cluster12, features = rownames(MM.cluster12.NRvsR.DE.sig)[1:20], slot = "scale.data", hjust = 0, angle=90, size=4)

dev.off()


MM.cluster12.overlap.cluster.pheno.markers <- intersect(rownames(MM.cluster12.de.markers),rownames(MM.cluster12.NRvsR.DE.sig))

length(MM.cluster12.overlap.cluster.pheno.markers)


tiff(file = "MM_cluster12_Heatmap_cluster.pheno.markers.tiff", width = 3000, height = 2000, units = "px", res = 600)

DoHeatmap(object = MM_cluster12, features = MM.cluster12.overlap.cluster.pheno.markers, slot = "scale.data", hjust = 0, angle=90, size=3)

dev.off()


# Feature plot of MM_cluster12 genes

tiff('FeaturePlot_MM_cluster12_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = rownames(MM.cluster12.de.markers), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Feature plot of MM_cluster12 genes

tiff('FeaturePlot_MM_cluster12_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = rownames(MM.cluster12.de.markers), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Feature plot of MM_cluster12_selected_markers genes

MM.cluster12.de.markers.selected <- c("TREM2","SPP1","RNASE1","MT1G","SEPP1","FOLR2","NUPR1","KLHDC8B","CCL18","MMP12","APOC2")

length(MM.cluster12.de.markers.selected)


tiff('FeaturePlot_MM_cluster12_selected_markers_Res.tiff', width = 6000, height = 6000, units = "px", res = 300)

FeaturePlot(immune.combined.more.clusters.Res, features = MM.cluster12.de.markers.selected, min.cutoff = "q99", label=T, label.size = 4, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_MM_cluster12_selected_markers_NonRes.tiff', width = 6000, height = 6000, units = "px", res = 300)

FeaturePlot(immune.combined.more.clusters.NonRes, features = MM.cluster12.de.markers.selected, min.cutoff = "q9", label=T, label.size = 4, cols = c("lightgrey", "red"))

dev.off()


DefaultAssay(object = MM_cluster6_12_23) <- "RNA"

MM_cluster6_12_23_rescale <- ScaleData(object = MM_cluster6_12_23, features = rownames(MM_cluster6_12_23))

DefaultAssay(object = MM_cluster6_12_23_rescale)

# DoHeatmap for re-selected genes for cluster 6 and 12

topMMgenes.reselected <- unique(c(MM.cluster6.de.markers.V2.topFCgenes,MM.cluster12.de.markers.selected,MM.cluster23.de.markers.V2.topFCgenes))

topMMgenes.reselected <- topMMgenes.reselected[!topMMgenes.reselected %in% c('RPS4Y1','TRBC2')]

topMMgenes.reselected <- c(topMMgenes.reselected,'CXCL10')

topMMgenes.reselected <- topMMgenes.reselected[c(1:7,37,8:36)]

length(topMMgenes.reselected)

topMMgenes.reselected

#MM.3clusters.topMMgenes.reselected <- MM_cluster6_12_23


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.topMMgenes.reselected.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = topMMgenes.reselected, hjust = 0.5, angle = 0) + NoLegend()

dev.off()


## add some markers for cluster12

topMMgenes.reselected.V2 <- c(topMMgenes.reselected,'C3', 'C1QA', 'C1QB','C1QC')

topMMgenes.reselected.V2 <- topMMgenes.reselected.V2[c(2,1,3:20,38:41,21:37)]

topMMgenes.reselected.V2


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.topMMgenes.reselected.V2.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = topMMgenes.reselected.V2, hjust = 0.5, angle = 0) + NoLegend()

dev.off()



## add some markers for cluster23

topMMgenes.reselected.V3 <- unique(c(topMMgenes.reselected.V2,
'CD8A', 'CD8B', 'IFITM1', 'CD3G', 'CD3D', 'KLRK1', 'CD96', 'TRAC', 'C5', 'C2', 'C4A', 'C4B'))

setdiff(c('CD8A', 'CD8B', 'IFITM1', 'CD3G', 'CD3D', 'KLRK1', 'CD96', 'TRAC'),topMMgenes.reselected.V2)


topMMgenes.reselected.V3

intersect(topMMgenes.reselected.V3,c('CD8A', 'IGKV4-1', 'IGKV2D-30', 'IGHV4-39', 'TRAV19', 'KLRD1', 'IGHV3-11', 'CD8B', 'IGLC1', 'IFITM1', 'CD3G', 'CD3D', 'IGHV2-5', 'KLRK1', 'CD96', 'TRAC', 'IGHV1-46'))

write.table(topMMgenes.reselected.V3, 'Fig3B_topMMgenes.reselected.V3.txt', sep='\t', row.names=F, col.names=F, quote=F)


### To draw heatmap for the 3 MM clusters

DefaultAssay(object = MM_cluster6_12_23) <- "integrated"

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.topMMgenes.reselected.V3.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = topMMgenes.reselected.V3, hjust = 0.5, angle = 0) + NoLegend()

dev.off()


TREM2.correlated.markers <- topMMgenes.reselected.V3[10:24]

TREM2.correlated.markers



######### relationship of 3 clusters with M1 and M2 polarization cells

M1.markers <- c('Azin1','Cd38','Cd86','Cxcl10','Cxcl9','Fpr2','Gpr18','Il12b','Il18','Il1b','Irf5','Nfkbiz','Nos2','Ptgs2','Socs3','Tlr4','Tnf')

M2.markers <- c('Alox15','Arg1','Chil3','Chil4','Clec7a','Egr2','Il10','Irf4','Klf4','Mrc1','Myc','Socs2','Tgm2')

genesMatch.M1 = getLDS(attributes = c('mgi_symbol'), filters = 'mgi_symbol', values = M1.markers, 
mart = mouse, attributesL = 'hgnc_symbol', martL = human, uniqueRows=T)

dim(genesMatch.M1)

head(genesMatch.M1)

genesMatch.M2 = getLDS(attributes = c('mgi_symbol'), filters = 'mgi_symbol', values = M2.markers, 
mart = mouse, attributesL = 'hgnc_symbol', martL = human, uniqueRows=T)

dim(genesMatch.M2)

head(genesMatch.M2)

load('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/CombinedSignatures.RData')

str(CombinedSignatures)


M1.genes.combined <- as.character(unique(c(genesMatch.M1$HGNC.symbol, CombinedSignatures$M1.Macrophage.Polarization, CombinedSignatures$M1.Macrophage)))

length(M1.genes.combined)

M1.genes.combined

M2.genes.combined <- as.character(unique(c(genesMatch.M2$HGNC.symbol, CombinedSignatures$M2.Macrophage.Polarization, CombinedSignatures$M2.Macrophage)))

length(M2.genes.combined)

M2.genes.combined

commonMacGenes <- intersect(M1.genes.combined, M2.genes.combined)


M1.genes.combined <- setdiff(M1.genes.combined,commonMacGenes)

length(M1.genes.combined)

M1.genes.combined


M2.genes.combined <- setdiff(M2.genes.combined,commonMacGenes)

length(M2.genes.combined)

M2.genes.combined


### To draw heatmap for the M1.genes.combined

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.M1.genes.combined.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, features = M1.genes.combined, hjust = 0.5, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red"))

dev.off()


### To draw heatmap for the M2.genes.combined

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.M2.genes.combined.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = M2.genes.combined, hjust = 0.5, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red"))

dev.off()


### DE genes for cluster 6, 12, 23 for the object MM_cluster6_12_23_rescale

table(Idents(MM_cluster6_12_23_rescale))

saveRDS(MM_cluster6_12_23_rescale, file='MM_cluster6_12_23_rescale.rds')


MM.cluster12.rescale.noRev.DE <- FindMarkers(object = MM_cluster6_12_23_rescale, ident.1 = "12", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster12.rescale.noRev.DE <- cbind(rownames(MM.cluster12.rescale.noRev.DE),MM.cluster12.rescale.noRev.DE)

colnames(MM.cluster12.rescale.noRev.DE)[1] <- 'Gene'

dim(MM.cluster12.rescale.noRev.DE)

head(MM.cluster12.rescale.noRev.DE)

tail(MM.cluster12.rescale.noRev.DE)

write.table(MM.cluster12.rescale.noRev.DE,'MM.cluster12.rescale.noRev.vsOthers.xls',sep='\t',row.names=F,quote=F)

#MM.cluster12.rescale.noRev.DE <- read.table('MM.cluster12.rescale.noRev.vsOthers.xls',sep='\t',header=T)



MM.cluster6.rescale.noRev.DE <- FindMarkers(object = MM_cluster6_12_23_rescale, ident.1 = "6", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster6.rescale.noRev.DE <- cbind(rownames(MM.cluster6.rescale.noRev.DE),MM.cluster6.rescale.noRev.DE)

colnames(MM.cluster6.rescale.noRev.DE)[1] <- 'Gene'

dim(MM.cluster6.rescale.noRev.DE)

head(MM.cluster6.rescale.noRev.DE)

tail(MM.cluster6.rescale.noRev.DE)

write.table(MM.cluster6.rescale.noRev.DE,'MM.cluster6.rescale.noRev.vsOthers.xls',sep='\t',row.names=F,quote=F)



MM.cluster23.rescale.noRev.DE <- FindMarkers(object = MM_cluster6_12_23_rescale, ident.1 = "23", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster23.rescale.noRev.DE <- cbind(rownames(MM.cluster23.rescale.noRev.DE),MM.cluster23.rescale.noRev.DE)

colnames(MM.cluster23.rescale.noRev.DE)[1] <- 'Gene'

dim(MM.cluster23.rescale.noRev.DE)

head(MM.cluster23.rescale.noRev.DE)

tail(MM.cluster23.rescale.noRev.DE)

write.table(MM.cluster23.rescale.noRev.DE,'MM.cluster23.rescale.noRev.vsOthers.xls',sep='\t',row.names=F,quote=F)


###### Developing MM.cluster12 signatures

dim(MM.cluster12.rescale.noRev.DE)

head(MM.cluster12.rescale.noRev.DE)


MM.cluster12.rescale.noRev.DE.sig <- MM.cluster12.rescale.noRev.DE[MM.cluster12.rescale.noRev.DE$p_val_adj < 0.05,]

MM.cluster12.rescale.noRev.DE.sig <- MM.cluster12.rescale.noRev.DE.sig[!is.na(MM.cluster12.rescale.noRev.DE.sig$Gene),]

dim(MM.cluster12.rescale.noRev.DE.sig)

head(MM.cluster12.rescale.noRev.DE.sig)

tail(MM.cluster12.rescale.noRev.DE.sig)

write.table(MM.cluster12.rescale.noRev.DE.sig,'MM.cluster12.rescale.noRev.DE.sig.txt',sep='\t',row.names=F,quote=F)



######### M1 genes for rescale.noRev.DE

MM.cluster6.rescale.noRev.DE.M1.genes <- MM.cluster6.rescale.noRev.DE[MM.cluster6.rescale.noRev.DE$Gene %in% M1.genes.combined,]

dim(MM.cluster6.rescale.noRev.DE.M1.genes)

head(MM.cluster6.rescale.noRev.DE.M1.genes)


MM.cluster12.rescale.noRev.DE.M1.genes <- MM.cluster12.rescale.noRev.DE[MM.cluster12.rescale.noRev.DE$Gene %in% M1.genes.combined,]

dim(MM.cluster12.rescale.noRev.DE.M1.genes)

head(MM.cluster12.rescale.noRev.DE.M1.genes)


MM.cluster23.rescale.noRev.DE.M1.genes <- MM.cluster23.rescale.noRev.DE[MM.cluster23.rescale.noRev.DE$Gene %in% M1.genes.combined,]

dim(MM.cluster23.rescale.noRev.DE.M1.genes)

head(MM.cluster23.rescale.noRev.DE.M1.genes)


MM.cluster6.rescale.noRev.DE.M1.genes[MM.cluster6.rescale.noRev.DE.M1.genes$avg_logFC > 0 & MM.cluster6.rescale.noRev.DE.M1.genes$p_val_adj < 0.05,]

MM.cluster23.rescale.noRev.DE.M1.genes[MM.cluster23.rescale.noRev.DE.M1.genes$avg_logFC > 0 & MM.cluster23.rescale.noRev.DE.M1.genes$p_val_adj < 0.05,]


M1.genes.selected <- c('IL1B','IDO1','IRF1','PTGS2','FPR2','CXCL11')



######### M2 genes for MM.cluster12.rescale.noRev.DE


MM.cluster12.rescale.noRev.DE.M2.genes <- MM.cluster12.rescale.noRev.DE[MM.cluster12.rescale.noRev.DE$Gene %in% M2.genes.combined,]

dim(MM.cluster12.rescale.noRev.DE.M2.genes)

head(MM.cluster12.rescale.noRev.DE.M2.genes)


MM.cluster12.rescale.noRev.DE.M2.genes[MM.cluster12.rescale.noRev.DE.M2.genes$avg_logFC > 0 & MM.cluster12.rescale.noRev.DE.M2.genes$p_val_adj < 0.05,]

M2.genes.selected <- as.character(MM.cluster12.rescale.noRev.DE.M2.genes$Gene)

M2.genes.selected

intersect(M1.genes.selected,M2.genes.selected)


M1.M2.genes.selected <- c(M1.genes.selected, M2.genes.selected)

M1.M2.genes.selected


M1.M2.genes.selected.V2 <- M1.M2.genes.selected[c(1:15,22,25)]

M1.M2.genes.selected.V2


### HeatMap for M1 and M2 polarization genes

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.SigDE.M1.M2.polarization.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = M1.M2.genes.selected.V2, hjust = 1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()

dev.off()


### HeatMap for TREM2 macrophages for Complement, M2 polarization genes

#TREM2.Complement.M2.genes <- c(M1.M2.genes.selected.V2[1:6], 'TREM2','C1QA','C1QB','C1QC','C2','C3','C4A','C4B','C5',M1.M2.genes.selected.V2[7:17])

TREM2.Complement.M2.genes <- c('TREM2','C1QA','C1QB','C1QC','C2','C3',M1.M2.genes.selected.V2[7:17])

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.SigDE.TREM2.Complement.M2.tiff", width = 3000, height = 3000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = TREM2.Complement.M2.genes, hjust = 1, angle = 0)  + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()

dev.off()

all(TREM2.Complement.M2.genes %in% MM.cluster12.rescale.noRev.DE.sig$Gene)

TREM2.related.genes <- as.character(unique(c(TREM2.correlated.markers, TREM2.Complement.M2.genes)))

length(TREM2.related.genes)

TREM2.related.genes

write.table(TREM2.related.genes,'TREM2.related.genes.txt',sep='\t',quote=F,row.names=F,col.names=F)


### HeatMap for TREM2sig.overlap.BMS038DE.GSE78220DE

#TREM2sig.overlap.BMS038DE.GSE78220DE <- read.table('TREM2sig.overlap.BMS038DE.GSE78220DE.txt',sep='\t',header=F)

#TREM2sig.overlap.BMS038DE.GSE78220DE <- as.character(unlist(TREM2sig.overlap.BMS038DE.GSE78220DE[[1]]))


TREM2sig.overlap.BMS038DE.GSE78220DE <- c("TREM2", "SPP1", "RNASE1", "MT1G", "SEPP1", "FOLR2", "NUPR1", "KLHDC8B", "CCL18", "MMP12", "APOC2", "C3", "C1QA", "C1QB", "C1QC", "C2",             
"MMP14", "CD276", "FN1", "MRC1", "CCL13", "LYVE1", "PDCD1LG2", "MMP9", "TGFB2", "ARG2", "CYTL1", "CLDN7", "TMEM171", "PRUNE2", "ITGA3", "STC1", "LY86-AS1", "TM4SF19", "ALDH1L2",           
"TRAF3IP2", "MAP2K5", "ZNF219", "TSHZ3", "STC2")   


length(TREM2sig.overlap.BMS038DE.GSE78220DE)

TREM2sig.overlap.BMS038DE.GSE78220DE


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.TREM2sig.overlap.BMS038DE.GSE78220DE.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = TREM2sig.overlap.BMS038DE.GSE78220DE, hjust = 1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()

dev.off()


#tiff("MM_3clusters.TREM2sig.overlap.BMS038DE.GSE78220DE.V2.tiff", width = 5000, height = 5000, units = "px", res = 600)
#DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = TREM2sig.overlap.BMS038DE.GSE78220DE, hjust = 1, angle = 0) + NoLegend()
#dev.off()


# TREM2sig.overlap.BMS038DE.GSE78220DE in MM.cluster12.rescale.noRev.DE

MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE <- MM.cluster12.rescale.noRev.DE[MM.cluster12.rescale.noRev.DE$Gene %in% TREM2sig.overlap.BMS038DE.GSE78220DE,]

MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE$Gene <- factor(MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE$Gene, levels=TREM2sig.overlap.BMS038DE.GSE78220DE)

MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE <- MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE[order(MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE$Gene),]

dim(MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE)

MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE

write.table(MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE, 'MM.cluster12.rescale.noRev.DE_TREM2sig.overlap.BMS038DE.GSE78220DE.xls',sep='\t',row.names=F,quote=F)


##################

BMS038_top51_100_DEgenes <- read.table('BMS038_top51_100_DEgenes.xls', sep='\t', header=F)

BMS038_top51_100_DEgenes <- as.character(unlist(BMS038_top51_100_DEgenes[[1]]))

length(BMS038_top51_100_DEgenes)

BMS038_top51_100_DEgenes

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.BMS038_top51_100_DEgenes.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = BMS038_top51_100_DEgenes, hjust = 1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()

dev.off()


##################

BMS038_top144_DEgenes_diff_TREM2 <- read.table('BMS038_top144_DEgenes_diff_TREM2.xls', sep='\t', header=F)

BMS038_top144_DEgenes_diff_TREM2 <- as.character(unlist(BMS038_top144_DEgenes_diff_TREM2[[1]]))

length(BMS038_top144_DEgenes_diff_TREM2)

BMS038_top144_DEgenes_diff_TREM2


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.BMS038_top144_DEgenes_diff_TREM2.tiff", width = 9000, height = 9000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = BMS038_top144_DEgenes_diff_TREM2, hjust = 1, angle = 0) + scale_fill_gradientn(colors = c("blue", "white", "red")) + NoLegend()

dev.off()



##################

BMS038_Select_V2 <- as.character(c('LINC00221','LGI4','RMRP','IL21','PGA3','PGA4','FCRL3','ALK','SYT6','PLA2G2D', 'RGS9'))

length(BMS038_Select_V2)

BMS038_Select_V2


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.BMS038_Select_V2.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23_rescale, assay = "RNA", features = BMS038_Select_V2, hjust = 0, angle = 0) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red")) 

dev.off()



########################################################################################################################################################

### VlnPlot for M2 markers


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

tiff("MM_3clusters.VlnPlot.M2.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23_rescale, features = M1.M2.genes.selected.V2[7:17]) 

dev.off()


### VlnPlot for TREM2.markers of cluster12

TREM2.markers <- c('TREM2', 'RNASE1', 'MT1G', 'SEPP1', 'FOLR2', 'NUPR1')

tiff("MM_3clusters.VlnPlot.TREM2.markers.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23_rescale, features = TREM2.markers, pt.size = 0.02, log = T) 

dev.off()


### VlnPlot for IDO1.markers of cluster6

TREM2.markers <- c('TREM2', 'RNASE1', 'MT1G', 'SEPP1', 'FOLR2', 'NUPR1')

tiff("MM_3clusters.VlnPlot.TREM2.markers.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23_rescale, features = TREM2.markers, pt.size = 0.02, log = T) 

dev.off()


tiff("MM_3clusters.VlnPlot.TREM2.no.rescale.markers.tiff", width = 5000, height = 5000, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23, features = TREM2.markers, pt.size = 0.02, log = T) 

dev.off()


### VlnPlot for IDO1.markers of cluster6

IDO1.markers <- c('IDO1', 'APOBEC3A', 'CXCL10')

tiff("MM_3clusters.VlnPlot.IDO1.markers.tiff", width = 5000, height = 2500, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23_rescale, features = IDO1.markers, pt.size = 0.02, log = T) 

dev.off()


### VlnPlot for Immunoregulatory.markers of cluster6

Immunoregulatory.markers <- c('TRAC', 'KLRK1', 'CD96')

tiff("MM_3clusters.VlnPlot.Immunoregulatory.markers.tiff", width = 5000, height = 2500, units = "px", res = 600)
  
VlnPlot(MM_cluster6_12_23_rescale, features = Immunoregulatory.markers, pt.size = 0.02, log = T) 

dev.off()



######### TREM positive cell percentage in each patient and in the scale of CD45 and CD68 (or MARCO): Start #########

Each.Responder.CellIDs <- readRDS('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Each.Responder.CellIDs.rds')

Each.NonResponder.CellIDs <- readRDS('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Each.NonResponder.CellIDs.rds')

str(Each.Responder.CellIDs)

sum(sapply(Each.Responder.CellIDs,length))

str(Each.NonResponder.CellIDs)

sum(sapply(Each.NonResponder.CellIDs,length))

sum(sapply(Each.Responder.CellIDs,length))+sum(sapply(Each.NonResponder.CellIDs,length))




############################## Redo clustering of immune.combined.more.clusters using revised cluster 6 and 12

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis/cluster12revised')

TREM2.cutoff <- 0

TREM2fetched.cluster12 <- as.data.frame(FetchData(MM_cluster12, vars='TREM2') > TREM2.cutoff)

TREM2fetched.cluster12$cellID <- rownames(TREM2fetched.cluster12)

class(TREM2fetched.cluster12)

dim(TREM2fetched.cluster12)

length(TREM2fetched.cluster12$cellID[TREM2fetched.cluster12$TREM2])

head(TREM2fetched.cluster12$cellID[TREM2fetched.cluster12$TREM2])

tail(TREM2fetched.cluster12$cellID[TREM2fetched.cluster12$TREM2])

TREM2positive <- TREM2fetched.cluster12$cellID[TREM2fetched.cluster12$TREM2]

length(TREM2positive)


table(Idents(MM_cluster12))

length(intersect(MM_cluster12_Responder_ID, TREM2positive))

length(intersect(MM_cluster12_NonResponder_ID, TREM2positive))

length(intersect(MM_cluster12_Responder_ID, TREM2positive))/length(MM_cluster12_Responder_ID)

length(intersect(MM_cluster12_NonResponder_ID, TREM2positive))/length(MM_cluster12_NonResponder_ID)


MM_cluster12_Responder_ID_TREM2pos <- intersect(MM_cluster12_Responder_ID, TREM2positive)

MM_cluster12_Responder_ID_TREM2pos


dim(DefineT_PhenoInfo_V2[DefineT_PhenoInfo_V2$Cell.Name %in% MM_cluster12_Responder_ID_TREM2pos,])

DefineT_PhenoInfo_V2[DefineT_PhenoInfo_V2$Cell.Name %in% MM_cluster12_Responder_ID_TREM2pos,]

MM_cluster12_Responder_TREM2pos_Info <- DefineT_PhenoInfo_V2[DefineT_PhenoInfo_V2$Cell.Name %in% MM_cluster12_Responder_ID_TREM2pos,]

MM_cluster12_Responder_TREM2pos_Info$SampleID <- as.character(MM_cluster12_Responder_TREM2pos_Info$SampleID)

table(MM_cluster12_Responder_TREM2pos_Info$SampleID)

MM_cluster12_Responder_ID_TREM2pos_change_to6 <- 
MM_cluster12_Responder_TREM2pos_Info$Cell.Name[MM_cluster12_Responder_TREM2pos_Info$SampleID %in% c('Post_P17','Post_P19','Pre_P26','Pre_P29')]

length(MM_cluster12_Responder_ID_TREM2pos_change_to6)

MM_cluster12_Responder_ID_TREM2pos_change_to6


immune.combined.more.clusters.V2 <- immune.combined.more.clusters

immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters <- immune.combined.more.clusters.V2@meta.data$RNA_snn_res.1.V2


dim(immune.combined.more.clusters.V2@meta.data)

head(immune.combined.more.clusters.V2@meta.data)

tail(immune.combined.more.clusters.V2@meta.data)

table(immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters)



immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters <- as.character(immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters)

immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters[rownames(immune.combined.more.clusters.V2@meta.data) %in% MM_cluster12_Responder_ID_TREM2pos_change_to6] <- '6'

class(immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters[rownames(immune.combined.more.clusters.V2@meta.data) %in% MM_cluster12_Responder_ID_TREM2pos_change_to6])

unique(immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters[rownames(immune.combined.more.clusters.V2@meta.data) %in% MM_cluster12_Responder_ID_TREM2pos_change_to6])


Idents(immune.combined.more.clusters.V2) <- immune.combined.more.clusters.V2@meta.data$cluster12revised_clusters

Idents(immune.combined.more.clusters.V2) <- factor(Idents(immune.combined.more.clusters.V2),levels=as.character(1:23))


class(Idents(immune.combined.more.clusters.V2))

table(Idents(immune.combined.more.clusters.V2))

table(Idents(immune.combined.more.clusters))


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis/cluster12revised')
 
tiff(file = "ImmuneCells_Clustering_cluster12_revised.tiff", width = 12000, height = 7000, units = "px", res = 600)

p1 <- DimPlot(object = immune.combined.more.clusters.V2, reduction = "umap", group.by = "Pheno", label.size = 8)
p2 <- DimPlot(object = immune.combined.more.clusters.V2, reduction = "umap", label = TRUE, label.size = 8)
plot_grid(p1, p2)

dev.off()


tiff(file = "ImmuneCells_Clustering_sidebyside_cluster12_revised.tiff", width = 12000, height = 9000, units = "px", res = 600)

DimPlot(object = immune.combined.more.clusters.V2, reduction = "umap", split.by = "Pheno", label = TRUE, label.size = 9)

dev.off()




########### calculate proportion of each cluster in R and NR for cluster 12 revised clustering

# immune.combined.more.clusters cell type plot

#CD8.exp <- expression(paste("CD",8^"+" ," T cells"))

immune.combined.more.clusters.V2@meta.data$Cell.Type.V2 <- immune.combined.more.clusters.V2@meta.data$Cell.Type

immune.combined.more.clusters.V2@meta.data$Cell.Type.V2 <- as.character(immune.combined.more.clusters.V2@meta.data$Cell.Type.V2)

immune.combined.more.clusters.V2@meta.data$Cell.Type.V2[immune.combined.more.clusters.V2@meta.data$Cell.Type.V2=='Macrophages/Monocytes'] <- 'Macrophages'

immune.combined.more.clusters.V2@meta.data$Cell.Type.V2[immune.combined.more.clusters.V2@meta.data$Cell.Type.V2=='MKI67hi Lymph.'] <- 'MKI67hi'



sort(unique(immune.combined.more.clusters.V2@meta.data$Cell.Type.V2))


Idents(immune.combined.more.clusters.V2) <- immune.combined.more.clusters.V2@meta.data$Pheno

immune.combined.more.clusters.V2.Res <- subset(immune.combined.more.clusters.V2, idents = "Responder")

immune.combined.more.clusters.V2.NonRes <- subset(immune.combined.more.clusters.V2, idents = "NonResponder")

str(immune.combined.more.clusters.V2.Res@meta.data)

str(immune.combined.more.clusters.V2.NonRes@meta.data)

Idents(immune.combined.more.clusters.V2) <- "Cell.Type.V2"

Idents(immune.combined.more.clusters.V2.Res) <- "Cell.Type.V2"

Idents(immune.combined.more.clusters.V2.NonRes) <- "Cell.Type.V2"


tiff(file = "ImmuneCells_MoreCellClusters_V2_CellType_DimPlot.tiff", width = 7000, height = 7000, units = "px", res = 600)

DimPlot(immune.combined.more.clusters.V2, reduction = "umap", label = TRUE, label.size = 9) + NoLegend()

dev.off()



tiff(file = "ImmuneCells_Res_MoreCellClusters_V2_CellType_DimPlot.tiff", width = 7000, height = 7000, units = "px", res = 600)

DimPlot(immune.combined.more.clusters.V2.Res, reduction = "umap", label = TRUE, label.size = 9) + NoLegend()

dev.off()


tiff(file = "ImmuneCells_NonRes_MoreCellClusters_V2_CellType_DimPlot.tiff", width = 7000, height = 7000, units = "px", res = 600)

DimPlot(immune.combined.more.clusters.V2.NonRes, reduction = "umap", label = TRUE, label.size = 9) + NoLegend()

dev.off()


DefaultAssay(object = immune.combined.more.clusters.V2) <- "RNA"


# Feature plot of canonical markers

tiff('immune.combined.more.clusters.V2.canonical.tiff', width = 7000, height = 7000, units = "px", res = 600)

FeaturePlot(immune.combined.more.clusters.V2, features = c("CD3D", "CD3E", "CD3G","CD8A","CD4","FOXP3","MKI67","NCR1","NCAM1","CD19","MZB1","MARCO","MERTK","FCER1A","NAPSA"), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


###

Idents(immune.combined.more.clusters.V2.Res) <- 'cluster12revised_clusters'

table(Idents(immune.combined.more.clusters.V2.Res))

sum(table(Idents(immune.combined.more.clusters.V2.Res)))

table(Idents(immune.combined.more.clusters.V2.Res))/sum(table(Idents(immune.combined.more.clusters.V2.Res)))

table(Idents(immune.combined.more.clusters.V2.Res))/nrow(immune.combined@meta.data)


df.table1 <- table(Idents(immune.combined.more.clusters.V2.Res))/sum(table(Idents(immune.combined.more.clusters.V2.Res)))

Res.proportions.AllClusters.cluster12.revised <- as.data.frame(df.table1)

Res.proportions.AllClusters.cluster12.revised$Phenotype <- "Responder/Regression"

colnames(Res.proportions.AllClusters.cluster12.revised)[1:2] <- c("Cluster","Percentage")

dim(Res.proportions.AllClusters.cluster12.revised)

head(Res.proportions.AllClusters.cluster12.revised)

tail(Res.proportions.AllClusters.cluster12.revised)


###

Idents(immune.combined.more.clusters.V2.NonRes) <- 'cluster12revised_clusters'

table(Idents(immune.combined.more.clusters.V2.NonRes))

sum(table(Idents(immune.combined.more.clusters.V2.NonRes)))

table(Idents(immune.combined.more.clusters.V2.NonRes))/sum(table(Idents(immune.combined.more.clusters.V2.NonRes)))

table(Idents(immune.combined.more.clusters.V2.NonRes))/nrow(immune.combined@meta.data)


df.table2 <- table(Idents(immune.combined.more.clusters.V2.NonRes))/sum(table(Idents(immune.combined.more.clusters.V2.NonRes)))

NonRes.proportions.AllClusters.cluster12.revised <- as.data.frame(df.table2)

NonRes.proportions.AllClusters.cluster12.revised$Phenotype <- "Non-responder/Progression"

colnames(NonRes.proportions.AllClusters.cluster12.revised)[1:2] <- c("Cluster","Percentage")

dim(NonRes.proportions.AllClusters.cluster12.revised)

head(NonRes.proportions.AllClusters.cluster12.revised)

tail(NonRes.proportions.AllClusters.cluster12.revised)


# Proportion of all clusters

RvsNR.proportions.AllClusters.cluster12.revised <- rbind(Res.proportions.AllClusters.cluster12.revised, NonRes.proportions.AllClusters.cluster12.revised)

dim(RvsNR.proportions.AllClusters.cluster12.revised)

head(RvsNR.proportions.AllClusters.cluster12.revised)

tail(RvsNR.proportions.AllClusters.cluster12.revised)


RvsNR.proportions.AllClusters.cluster12.revised$Cluster <- factor(RvsNR.proportions.AllClusters.cluster12.revised$Cluster,levels=as.character(1:23))

levels(RvsNR.proportions.AllClusters.cluster12.revised$Cluster)


sort(unique(RvsNR.proportions.AllClusters.cluster12.revised$Phenotype))

RvsNR.proportions.AllClusters.cluster12.revised$Phenotype[RvsNR.proportions.AllClusters.cluster12.revised$Phenotype=='Non-responder/Progression'] <- 'NR'

RvsNR.proportions.AllClusters.cluster12.revised$Phenotype[RvsNR.proportions.AllClusters.cluster12.revised$Phenotype=='Responder/Regression'] <- 'R'

RvsNR.proportions.AllClusters.cluster12.revised$Phenotype <- factor(RvsNR.proportions.AllClusters.cluster12.revised$Phenotype,levels=c('R','NR'))


tiff('Barplot.RvsNR.proportions.cluster12.revised.AllClusters.cluster12.revised.tiff', width = 5000, height = 4000, units = "px", res = 600)

ggplot(RvsNR.proportions.AllClusters.cluster12.revised,aes(x=factor(Phenotype), y = Percentage)) + 
    facet_wrap(~Cluster) + theme(strip.text.x = element_text(size = 20, colour = "red", angle = 0)) +
    labs(y = "Percentage of CD45+ immune cells") +
    geom_bar(aes(fill = factor(Phenotype)),stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

dev.off()


write.table(RvsNR.proportions.AllClusters.cluster12.revised, 'RvsNR.proportions.cluster12.revised.AllClusters.cluster12.revised.xls', sep='\t', row.names=F, quote=F)


# plot fold changes of different clusters between NR and R groups

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis/cluster12revised')

data <-read.table("RvsNR.proportions.cluster12.revised.FoldChange.NRvsR.txt",sep='\t',header=T) 

tiff('Barplot.AllClusters.cluster12.revised.NRvsR.proportions.FoldChange.tiff',width = 6000, height = 4000, units = "px", res = 600)

bp <- barplot(data$NRvsR_FoldChange, main = "", xlab = "Single-cell Cluster", 
ylab = "Fold change of percentage of CD45+ immune cells (NR/R)", ylim = c(-10,16),names.arg = data$Cluster, axes=F,
    col = rainbow(nrow(data)))

xval = seq(-10,16,by=2)

axis(side = 2, at = xval, labels = FALSE, xpd=T)
axis(side = 2, at = xval, tick = FALSE, labels = xval, xpd=T)

abline(h=0)
abline(h=6,lty=2)
abline(h=-6,lty=2)

dev.off()


# Macrophage analyses with cluster 12 revised

Idents(immune.combined.more.clusters.V2) <- 'cluster12revised_clusters'

table(Idents(immune.combined.more.clusters.V2))

sum(table(Idents(immune.combined.more.clusters.V2)))


MM_3clusters_cluster12.revised <- subset(immune.combined.more.clusters.V2, idents = c("6","12","23"))

str(MM_3clusters_cluster12.revised@meta.data)

dim(MM_3clusters_cluster12.revised@meta.data)

table(Idents(MM_3clusters_cluster12.revised))

sort(unique(MM_3clusters_cluster12.revised@meta.data$Cell.Type))

Idents(MM_3clusters_cluster12.revised) <- factor(Idents(MM_3clusters_cluster12.revised),levels=c('6','12','23'))


# find markers for every cluster compared to all remaining cells, report only the positive ones in MM population - 3 clusters of 6, 12, 23

MM_3clusters_cluster12.revised.markers <- FindAllMarkers(MM_3clusters_cluster12.revised, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

MM_3clusters_cluster12.revised.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)


top20 <- MM_3clusters_cluster12.revised.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

tiff("MM_3clusters_cluster12.revised.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_3clusters_cluster12.revised, features = top20$gene, hjust = 0.5, angle = 0) + NoLegend()

dev.off()

dim(MM_3clusters_cluster12.revised.markers)

levels(MM_3clusters_cluster12.revised)


### DE analysis of each of the 3 MM clusters

# MM.cluster6.cluster12.revised DE analysis

MM.cluster6.cluster12.revised.de.markers <- FindMarkers(MM_3clusters_cluster12.revised, ident.1 = "6", ident.2 = NULL, only.pos = TRUE)

dim(MM.cluster6.cluster12.revised.de.markers)

MM.cluster6.cluster12.revised.de.markers

MM.cluster6.cluster12.revised.de.markers.V2 <- MM.cluster6.cluster12.revised.de.markers

MM.cluster6.cluster12.revised.de.markers.V2$FoldChange <- MM.cluster6.cluster12.revised.de.markers.V2$pct.1/MM.cluster6.cluster12.revised.de.markers.V2$pct.2

MM.cluster6.cluster12.revised.de.markers.V2 <- MM.cluster6.cluster12.revised.de.markers.V2[order(MM.cluster6.cluster12.revised.de.markers.V2$FoldChange, decreasing=T),]

range(MM.cluster6.cluster12.revised.de.markers.V2$p_val)

MM.cluster6.cluster12.revised.de.markers.V2

MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes <- rownames(MM.cluster6.cluster12.revised.de.markers.V2[MM.cluster6.cluster12.revised.de.markers.V2$FoldChange > 1.5,])

MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes <- MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes[!MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes %in% c('ISG20','BIRC3','CLEC10A','FCN1','CD52','CXCL10')]

MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes



MM.cluster6.cluster12.revised.DE <- FindMarkers(object = MM_3clusters_cluster12.revised, ident.1 = "6", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster6.cluster12.revised.DE <- cbind(rownames(MM.cluster6.cluster12.revised.DE),MM.cluster6.cluster12.revised.DE)

colnames(MM.cluster6.cluster12.revised.DE)[1] <- 'Gene'

dim(MM.cluster6.cluster12.revised.DE)

head(MM.cluster6.cluster12.revised.DE)

tail(MM.cluster6.cluster12.revised.DE)


write.table(MM.cluster6.cluster12.revised.DE,'MM.cluster6vsOthers.DEresults.xls',sep='\t',row.names=F,quote=F)


MM.cluster6.cluster12.revised.de.markers <- MM.cluster6.cluster12.revised.de.markers[order(MM.cluster6.cluster12.revised.de.markers$p_val),]

head(MM.cluster6.cluster12.revised.de.markers)


# Enriched cluster6 genes

range(MM.cluster6.cluster12.revised.DE$avg_logFC[MM.cluster6.cluster12.revised.DE$avg_logFC>0])

MM.cluster6.enriched <- MM.cluster6.cluster12.revised.DE[MM.cluster6.cluster12.revised.DE$avg_logFC > 0 & MM.cluster6.cluster12.revised.DE$p_val < 0.05,]

MM.cluster6.enriched <- MM.cluster6.enriched[order(-MM.cluster6.enriched$avg_logFC),]

dim(MM.cluster6.enriched)

head(MM.cluster6.enriched)

tail(MM.cluster6.enriched)

write.table(MM.cluster6.enriched,'MM.cluster6.enriched.xls',sep='\t',row.names=F,quote=F)


############################################################# fgsea analyses for MM cluster 6 first as an example #############################################################

library(fgsea)

library(data.table) 

library(ggplot2)

library('org.Hs.eg.db')

library(tidyverse)

######## fgsea based on gene symbol ########


# Load the pathways into a named list

pathways.hallmark <- gmtPathways("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Seurat3_Analysis/msigdb/h.all.v6.2.symbols.gmt")
pathways.kegg <- gmtPathways("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Seurat3_Analysis/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
pathways.reactome <- gmtPathways("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Seurat3_Analysis/msigdb/c2.cp.reactome.v6.2.symbols.gmt")
pathways.cp <- gmtPathways("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Seurat3_Analysis/msigdb/c2.cp.v6.2.symbols.gmt")
pathways.go <- gmtPathways("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/Seurat3_Analysis/msigdb/c5.all.v6.2.symbols.gmt")

length(pathways.hallmark)
length(pathways.kegg)
length(pathways.reactome)
length(pathways.cp)
length(pathways.go)


pathways.combined <- c(pathways.hallmark, pathways.kegg, pathways.reactome, pathways.cp, pathways.go)

length(pathways.combined)

head(pathways.combined)

tail(pathways.combined)


MM.c6.DE <- MM.cluster6.cluster12.revised.DE

MM.c6.DE$stat <- (-log10(MM.c6.DE$p_val)*MM.c6.DE$avg_logFC)
#MM.c6.DE$Gene <- rownames(MM.c6.DE)

dim(MM.c6.DE)
head(MM.c6.DE)


res2.cluster <- MM.c6.DE %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
  
dim(res2.cluster)

class(res2.cluster)

head(res2.cluster)

names(res2.cluster)


ranks.cluster <- deframe(res2.cluster)
head(ranks.cluster, 20)


###################################### combined.cluster geneset that were significant for cluster #####################################

combined.cluster <- fgsea(pathways=pathways.combined, stats=ranks.cluster, minSize=15, maxSize=500, nperm=100000)

combined.cluster <- combined.cluster[order(combined.cluster$pval),]

class(combined.cluster)

dim(combined.cluster)

head(combined.cluster)

tail(combined.cluster)

range(combined.cluster$pval)



combined.cluster.sigPathwaysUp <- unique(combined.cluster[ES > 0][pval <= 0.05][order(pval),pathway])

length(combined.cluster.sigPathwaysUp)

combined.cluster.sigPathwaysUp

combined.cluster.sigPathwaysDown <- unique(combined.cluster[ES < 0][pval <= 0.05][order(pval),pathway])

length(combined.cluster.sigPathwaysDown)

combined.cluster.sigPathwaysDown


length(pathways.combined)

combined.cluster.sigPathways <- c(combined.cluster.sigPathwaysUp,combined.cluster.sigPathwaysDown)

pathways.combined.cluster.sig <- pathways.combined[names(pathways.combined) %in% combined.cluster.sigPathways]

pathways.combined.cluster.sig <- pathways.combined.cluster.sig[!duplicated(names(pathways.combined.cluster.sig))]

any(duplicated(names(pathways.combined.cluster.sig)))

length(pathways.combined.cluster.sig)

head(pathways.combined.cluster.sig)

tail(pathways.combined.cluster.sig)

names(pathways.combined.cluster.sig)


# fgsea only sig pathways above

MM.c6.DE.sigOnly <- fgsea(pathways=pathways.combined.cluster.sig, stats=ranks.cluster, minSize=15, maxSize=500, nperm=100000)


###

Tidy.cluster.sigOnly <- MM.c6.DE.sigOnly %>%
  as_tibble() %>%
  arrange(desc(NES))

outname1 <- 'MM.c6.DE'

tiff(paste0(outname1,'.sigOnly.NES.tiff'), width = 4000, height = 1800, units = "px", res = 300)
  
ggplot(Tidy.cluster.sigOnly, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES<0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Combined pathways NES from GSEA") + 
  theme_minimal()

dev.off()

cluster.sigOnly.genes <- pathways.combined.cluster.sig %>% 
  enframe("pathway", "Gene") %>% 
  unnest() %>% 
  inner_join(res2.cluster, by="Gene")

write.table(cluster.sigOnly.genes,paste0(outname1,'.sigOnly.genes.xls'),sep='\t',row.names=F,quote=F) 



###################### Dotplot analysis for cluster6

# create bubble plot for the top 10 pathways

MM.cluster6.Reactome.results <- read.csv('MM.cluster6.Reactome.results.revised.csv',sep=',',header=T)


MM.cluster6.Reactome.results$negativ.log10 <- -log10(MM.cluster6.Reactome.results$upload_1..raw.P.value.)

colnames(MM.cluster6.Reactome.results)[6] <- 'Fold.Enrichment'


dim(MM.cluster6.Reactome.results)

head(MM.cluster6.Reactome.results)

tail(MM.cluster6.Reactome.results)

nrow(MM.cluster6.Reactome.results)

#write.table(MM.cluster6.Reactome.results,'MM.cluster6.Reactome.results.revised.xls',sep='\t',row.names=F,quote=F)


MM.cluster6.Reactome.results$Reactome.pathways <- factor(MM.cluster6.Reactome.results$Reactome.pathways,levels=MM.cluster6.Reactome.results$Reactome.pathways)

class(MM.cluster6.Reactome.results$Reactome.pathways)

MM.cluster6.Reactome.results$Reactome.pathways

range(MM.cluster6.Reactome.results$Fold.Enrichment)

range(MM.cluster6.Reactome.results$negativ.log10)


tiff('Dotplot.MM.cluster6.Reactome.tiff', width = 12000, height = 12000, units = "px", res = 600)

sp <- ggplot(MM.cluster6.Reactome.results, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(1, 10))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=seq(0,30,by=5), limits=c(0, 30)) +
  theme(axis.text.x = element_text(angle = 90))
  			#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank()) +
 				# + scale_x_discrete(name="Reactome pathways",labels=c(rep("",nrow(MM.cluster6.Reactome.results))))  
 				  
dev.off()


# cluster6 top10 only

MM.cluster6.Reactome.results.top10 <- MM.cluster6.Reactome.results[1:10,]

range(MM.cluster6.Reactome.results.top10$Fold.Enrichment)


tiff('Dotplot.MM.cluster6.Reactome.V2.tiff', width = 4000, height = 3000, units = "px", res = 600)

sp <- ggplot(MM.cluster6.Reactome.results.top10, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(6, 12))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=seq(0,30,by=5), limits=c(0, 30)) +
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) 
         				  
dev.off()


#############################################################################################################################################

# MM.cluster12.cluster12.revised DE analysis

MM.cluster12.cluster12.revised.de.markers <- FindMarkers(MM_3clusters_cluster12.revised, ident.1 = "12", ident.2 = NULL, only.pos = TRUE)

dim(MM.cluster12.cluster12.revised.de.markers)

MM.cluster12.cluster12.revised.DE <- FindMarkers(object = MM_3clusters_cluster12.revised, ident.1 = "12", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster12.cluster12.revised.DE <- cbind(rownames(MM.cluster12.cluster12.revised.DE),MM.cluster12.cluster12.revised.DE)

colnames(MM.cluster12.cluster12.revised.DE)[1] <- 'Gene'

dim(MM.cluster12.cluster12.revised.DE)

head(MM.cluster12.cluster12.revised.DE)

tail(MM.cluster12.cluster12.revised.DE)

write.table(MM.cluster12.cluster12.revised.DE,'MM.cluster12vsOthers.DEresults.xls',sep='\t',row.names=F,quote=F)


MM.cluster12.cluster12.revised.de.markers <- MM.cluster12.cluster12.revised.de.markers[order(MM.cluster12.cluster12.revised.de.markers$p_val),]

head(MM.cluster12.cluster12.revised.de.markers)


# Enriched cluster12 genes

range(MM.cluster12.cluster12.revised.DE$avg_logFC[MM.cluster12.cluster12.revised.DE$avg_logFC>0])

MM.cluster12.enriched <- MM.cluster12.cluster12.revised.DE[MM.cluster12.cluster12.revised.DE$avg_logFC > 0 & MM.cluster12.cluster12.revised.DE$p_val < 0.05,]

MM.cluster12.enriched <- MM.cluster12.enriched[order(-MM.cluster12.enriched$avg_logFC),]

dim(MM.cluster12.enriched)

head(MM.cluster12.enriched)

tail(MM.cluster12.enriched)

write.table(MM.cluster12.enriched,'MM.cluster12.enriched.xls',sep='\t',row.names=F,quote=F)


###################### Dotplot analysis for cluster12

# create bubble plot for the top 10 pathways

MM.cluster12.Reactome.results <- read.csv('analysis.MM.cluster12.Reactome.results.csv',sep=',',header=T)

MM.cluster12.Reactome.results$Reactome.pathways <- gsub(' \\(.*$','',MM.cluster12.Reactome.results$Reactome.pathways)

MM.cluster12.Reactome.results$negativ.log10 <- -log10(MM.cluster12.Reactome.results$upload_1..raw.P.value.)

colnames(MM.cluster12.Reactome.results)[6] <- 'Fold.Enrichment'


MM.cluster12.Reactome.results[,c(1,6)]


MM.cluster12.Reactome.results <- MM.cluster12.Reactome.results[-c(6,11,12,14,16,22,24),]

dim(MM.cluster12.Reactome.results)

head(MM.cluster12.Reactome.results)

tail(MM.cluster12.Reactome.results)

nrow(MM.cluster12.Reactome.results)

write.table(MM.cluster12.Reactome.results,'MM.cluster12.Reactome.results.revised.xls',sep='\t',row.names=F,quote=F)


MM.cluster12.Reactome.results$Reactome.pathways <- factor(MM.cluster12.Reactome.results$Reactome.pathways,levels=MM.cluster12.Reactome.results$Reactome.pathways)

class(MM.cluster12.Reactome.results$Reactome.pathways)

MM.cluster12.Reactome.results$Reactome.pathways

range(MM.cluster12.Reactome.results$Fold.Enrichment)

range(MM.cluster12.Reactome.results$negativ.log10)


tiff('Dotplot.MM.cluster12.Reactome.tiff', width = 12000, height = 12000, units = "px", res = 600)

sp <- ggplot(MM.cluster12.Reactome.results, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(1, 10))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=c(0:12), limits=c(0, 12)) +
  theme(axis.text.x = element_text(angle = 90))
  			#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank()) +
 				# + scale_x_discrete(name="Reactome pathways",labels=c(rep("",nrow(MM.cluster12.Reactome.results))))  
 				  
dev.off()


# cluster12 top10 only

MM.cluster12.Reactome.results.top10 <- MM.cluster12.Reactome.results[1:10,]

range(MM.cluster12.Reactome.results.top10$Fold.Enrichment)


tiff('Dotplot.MM.cluster12.Reactome.V2.tiff', width = 4000, height = 3000, units = "px", res = 600)

sp <- ggplot(MM.cluster12.Reactome.results.top10, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(4, 12))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=c(0:12), limits=c(0, 12)) +
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) 
         				  
dev.off()


#############################################################################################################################################

# MM.cluster23.cluster12.revised DE analysis

MM.cluster23.cluster12.revised.de.markers <- FindMarkers(MM_3clusters_cluster12.revised, ident.1 = "23", ident.2 = NULL, only.pos = TRUE)

MM.cluster23.cluster12.revised.de.markers <- MM.cluster23.cluster12.revised.de.markers[MM.cluster23.cluster12.revised.de.markers$p_val < 0.05,]

dim(MM.cluster23.cluster12.revised.de.markers)

range(MM.cluster23.cluster12.revised.de.markers$p_val)


MM.cluster23.cluster12.revised.de.markers.V2 <- MM.cluster23.cluster12.revised.de.markers

MM.cluster23.cluster12.revised.de.markers.V2$FoldChange <- MM.cluster23.cluster12.revised.de.markers.V2$pct.1/MM.cluster23.cluster12.revised.de.markers.V2$pct.2

MM.cluster23.cluster12.revised.de.markers.V2 <- MM.cluster23.cluster12.revised.de.markers.V2[order(MM.cluster23.cluster12.revised.de.markers.V2$FoldChange, decreasing=T),]

range(MM.cluster23.cluster12.revised.de.markers.V2$p_val)

MM.cluster23.cluster12.revised.de.markers.V2

MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes <- rownames(MM.cluster23.cluster12.revised.de.markers.V2[MM.cluster23.cluster12.revised.de.markers.V2$FoldChange > 2 & MM.cluster23.cluster12.revised.de.markers.V2$pct.1 > 0.7,])

MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes <- MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes[!MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes %in% c('ISG20','BIRC3','CLEC10A','FCN1','CD52','CXCL10')]

MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes




MM.cluster23.cluster12.revised.DE <- FindMarkers(object = MM_3clusters_cluster12.revised, ident.1 = "23", ident.2 = NULL, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster23.cluster12.revised.DE <- cbind(rownames(MM.cluster23.cluster12.revised.DE),MM.cluster23.cluster12.revised.DE)

colnames(MM.cluster23.cluster12.revised.DE)[1] <- 'Gene'

dim(MM.cluster23.cluster12.revised.DE)

head(MM.cluster23.cluster12.revised.DE)

tail(MM.cluster23.cluster12.revised.DE)

write.table(MM.cluster23.cluster12.revised.DE,'MM.cluster23vsOthers.DEresults.xls',sep='\t',row.names=F,quote=F)




# Enriched cluster23 genes

range(MM.cluster23.cluster12.revised.DE$avg_logFC[MM.cluster23.cluster12.revised.DE$avg_logFC>0])

MM.cluster23.enriched <- MM.cluster23.cluster12.revised.DE[MM.cluster23.cluster12.revised.DE$avg_logFC > 0 & MM.cluster23.cluster12.revised.DE$p_val < 0.05,]

MM.cluster23.enriched <- MM.cluster23.enriched[order(-MM.cluster23.enriched$avg_logFC),]

dim(MM.cluster23.enriched)

head(MM.cluster23.enriched)

tail(MM.cluster23.enriched)

write.table(MM.cluster23.enriched,'MM.cluster23.enriched.xls',sep='\t',row.names=F,quote=F)


###################### Dotplot analysis for cluster23

# create bubble plot for the top 10 pathways

MM.cluster23.Reactome.results <- read.csv('analysis.MM.cluster23.Reactome.results.csv',sep=',',header=T)

MM.cluster23.Reactome.results$Reactome.pathways <- gsub(' \\(.*$','',MM.cluster23.Reactome.results$Reactome.pathways)

MM.cluster23.Reactome.results$negativ.log10 <- -log10(MM.cluster23.Reactome.results$upload_1..raw.P.value.)

colnames(MM.cluster23.Reactome.results)[6] <- 'Fold.Enrichment'


MM.cluster23.Reactome.results[,c(1,6)]


MM.cluster23.Reactome.results <- MM.cluster23.Reactome.results[-c(1,2,5,6,12,14:20,25:30,32,34:36,38,40,42,44:46,48,49,54:60,63,65),]

dim(MM.cluster23.Reactome.results)

head(MM.cluster23.Reactome.results)

tail(MM.cluster23.Reactome.results)

nrow(MM.cluster23.Reactome.results)

write.table(MM.cluster23.Reactome.results,'MM.cluster23.Reactome.results.revised.xls',sep='\t',row.names=F,quote=F)


MM.cluster23.Reactome.results$Reactome.pathways <- factor(MM.cluster23.Reactome.results$Reactome.pathways,levels=MM.cluster23.Reactome.results$Reactome.pathways)

class(MM.cluster23.Reactome.results$Reactome.pathways)

MM.cluster23.Reactome.results$Reactome.pathways

range(MM.cluster23.Reactome.results$Fold.Enrichment)

range(MM.cluster23.Reactome.results$negativ.log10)


tiff('Dotplot.MM.cluster23.Reactome.tiff', width = 12000, height = 12000, units = "px", res = 600)

sp <- ggplot(MM.cluster23.Reactome.results, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(1, 10))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=seq(0,8,by=2), limits=c(0,8)) +
  theme(axis.text.x = element_text(angle = 90))
  			#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #axis.text.x=element_blank()) +
 				# + scale_x_discrete(name="Reactome pathways",labels=c(rep("",nrow(MM.cluster23.Reactome.results))))  
 				  
dev.off()


# cluster23 top10 only

MM.cluster23.Reactome.results.top10 <- MM.cluster23.Reactome.results[1:10,]

range(MM.cluster23.Reactome.results.top10$Fold.Enrichment)


tiff('Dotplot.MM.cluster23.Reactome.V2.tiff', width = 4000, height = 3000, units = "px", res = 600)

sp <- ggplot(MM.cluster23.Reactome.results.top10, aes(x = Reactome.pathways, y = negativ.log10)) + 
  geom_point(color='blue',aes(size = Fold.Enrichment)) +
    scale_size(range = c(5, 9))  # Adjust the range of points size

sp + scale_x_discrete(name="Reactome pathways") +
	scale_y_continuous(name="-Log10(p value)", breaks=seq(0,8,by=2), limits=c(0,8)) +
  theme(#axis.title.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) 
         				  
dev.off()


#############################################################################################################################################

###### MM_cluster12.cluster12.revised NR vs R analysis

MM_cluster12.cluster12.revised <- subset(MM_3clusters_cluster12.revised, idents='12')

table(Idents(MM_cluster12.cluster12.revised))

str(MM_cluster12.cluster12.revised@meta.data)

Idents(MM_cluster12.cluster12.revised) <- MM_cluster12.cluster12.revised@meta.data$Pheno
		
table(Idents(MM_cluster12.cluster12.revised))

tiff(file = "MM_cluster12.cluster12.revised_Heatmap_NRvR_deGenes.tiff", width = 3000, height = 2500, units = "px", res = 600)

DoHeatmap(object = MM_cluster12.cluster12.revised, features = rownames(MM.cluster12.cluster12.revised.de.markers), slot = "scale.data", hjust = 0, angle=0, size=4)

dev.off()


MM_cluster12.cluster12.revised_Responder_ID <- rownames(MM_cluster12.cluster12.revised@meta.data[MM_cluster12.cluster12.revised@meta.data$Pheno=='Responder',])

length(MM_cluster12.cluster12.revised_Responder_ID)

MM_cluster12.cluster12.revised_Responder_ID


MM_cluster12.cluster12.revised_NonResponder_ID <- rownames(MM_cluster12.cluster12.revised@meta.data[MM_cluster12.cluster12.revised@meta.data$Pheno=='NonResponder',])

length(MM_cluster12.cluster12.revised_NonResponder_ID)

MM_cluster12.cluster12.revised_NonResponder_ID



MM.cluster12.cluster12.revised.NRvsR.DE <- FindMarkers(object = MM_cluster12.cluster12.revised, ident.1 = "NonResponder", ident.2 = "Responder", logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)

MM.cluster12.cluster12.revised.NRvsR.DE <- MM.cluster12.cluster12.revised.NRvsR.DE[order(MM.cluster12.cluster12.revised.NRvsR.DE$avg_logFC,decreasing=T),]

dim(MM.cluster12.cluster12.revised.NRvsR.DE)

head(MM.cluster12.cluster12.revised.NRvsR.DE)

tail(MM.cluster12.cluster12.revised.NRvsR.DE)

MM.cluster12.cluster12.revised.NRvsR.DE[rownames(MM.cluster12.cluster12.revised.NRvsR.DE)=='TREM2',]


MM.cluster12.cluster12.revised.NRvsR.DE.sig <- MM.cluster12.cluster12.revised.NRvsR.DE[MM.cluster12.cluster12.revised.NRvsR.DE$p_val < 0.05,]

dim(MM.cluster12.cluster12.revised.NRvsR.DE.sig)

head(MM.cluster12.cluster12.revised.NRvsR.DE.sig)

tail(MM.cluster12.cluster12.revised.NRvsR.DE.sig)


tiff(file = "MM_cluster12.cluster12.revised_Heatmap_NRvR_deGenes_DE.sig.tiff", width = 3000, height = 2500, units = "px", res = 600)

DoHeatmap(object = MM_cluster12.cluster12.revised, features = rownames(MM.cluster12.cluster12.revised.NRvsR.DE.sig)[1:20], slot = "scale.data", hjust = 0, angle=90, size=4)

dev.off()


MM.cluster12.cluster12.revised.overlap.cluster.pheno.markers <- intersect(rownames(MM.cluster12.cluster12.revised.de.markers),rownames(MM.cluster12.cluster12.revised.NRvsR.DE.sig))

length(MM.cluster12.cluster12.revised.overlap.cluster.pheno.markers)


tiff(file = "MM_cluster12.cluster12.revised_Heatmap_cluster.pheno.markers.tiff", width = 3000, height = 2000, units = "px", res = 600)

DoHeatmap(object = MM_cluster12.cluster12.revised, features = MM.cluster12.cluster12.revised.overlap.cluster.pheno.markers, slot = "scale.data", hjust = 0, angle=90, size=3)

dev.off()


# Feature plot of MM_cluster12.cluster12.revised genes

tiff('FeaturePlot_MM_cluster12.cluster12.revised_Res.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.Res, features = rownames(MM.cluster12.cluster12.revised.de.markers), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Feature plot of MM_cluster12.cluster12.revised genes

tiff('FeaturePlot_MM_cluster12.cluster12.revised_NonRes.tiff', width = 10000, height = 10000, units = "px", res = 300)

FeaturePlot(immune.combined.NonRes, features = rownames(MM.cluster12.cluster12.revised.de.markers), min.cutoff = "q9", label.size = 8, cols = c("lightgrey", "red"))

dev.off()


# Feature plot of MM_cluster12.cluster12.revised_selected_markers genes

MM.cluster12.cluster12.revised.de.markers.selected <- c("TREM2","SPP1","RNASE1","MT1G","SEPP1","FOLR2","NUPR1","KLHDC8B","CCL18","MMP12","APOC2")

length(MM.cluster12.cluster12.revised.de.markers.selected)


tiff('FeaturePlot_MM_cluster12.cluster12.revised_selected_markers_Res.tiff', width = 6000, height = 6000, units = "px", res = 300)

FeaturePlot(immune.combined.more.clusters.V2.Res, features = MM.cluster12.cluster12.revised.de.markers.selected, min.cutoff = "q99", label=T, label.size = 4, cols = c("lightgrey", "red"))

dev.off()


tiff('FeaturePlot_MM_cluster12.cluster12.revised_selected_markers_NonRes.tiff', width = 6000, height = 6000, units = "px", res = 300)

FeaturePlot(immune.combined.more.clusters.V2.NonRes, features = MM.cluster12.cluster12.revised.de.markers.selected, min.cutoff = "q9", label=T, label.size = 4, cols = c("lightgrey", "red"))

dev.off()


# DoHeatmap for re-selected genes for cluster 6 and 12

topMMgenes.reselected <- unique(c(MM.cluster6.cluster12.revised.de.markers.V2.topFCgenes,MM.cluster12.cluster12.revised.de.markers.selected,MM.cluster23.cluster12.revised.de.markers.V2.topFCgenes))

topMMgenes.reselected <- topMMgenes.reselected[!topMMgenes.reselected %in% c('RPS4Y1','TRBC2')]

length(topMMgenes.reselected)

topMMgenes.reselected

#MM.3clusters.topMMgenes.reselected <- MM_3clusters_cluster12.revised


tiff("MM_3clusters.topMMgenes.reselected.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_3clusters_cluster12.revised, features = topMMgenes.reselected, hjust = 0.5, angle = 0) + NoLegend()

dev.off()


# DoHeatmap for re-selected genes V2 for 3 clusters


cluster6.new <- rownames(MM.cluster6.cluster12.revised.de.markers[MM.cluster6.cluster12.revised.de.markers$avg_logFC > log2(1.5),])

cluster6.new <- c(cluster6.new, "FCER1A", "IDO1", "S100A12", "CD1C")

cluster6.new


tiff("MM_3clusters.cluster6.posMarkers.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_3clusters_cluster12.revised, features = cluster6.new, hjust = 0.5, angle = 0) + NoLegend()

dev.off()


cluster12.new <- rownames(MM.cluster12.cluster12.revised.de.markers[MM.cluster12.cluster12.revised.de.markers$avg_logFC > log2(1.5),])

cluster12.new <- c(cluster12.new, "KLHDC8B", "CCL18", "MMP12", "APOC2")

cluster12.new


cluster23.new <- rownames(MM.cluster23.cluster12.revised.de.markers[MM.cluster23.cluster12.revised.de.markers$avg_logFC > log2(2.5),])

cluster23.new


tiff("MM_3clusters.cluster12.posMarkers.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_3clusters_cluster12.revised, features = cluster12.new, hjust = 0.5, angle = 0) + NoLegend()

dev.off()


MM.3clusters.new <- c(cluster6.new, cluster12.new, cluster23.new)

MM.3clusters.new <- MM.3clusters.new[!MM.3clusters.new %in% c('FCN1','IL1B','CLEC10A','SELL')]

MM.3clusters.new <- MM.3clusters.new[!MM.3clusters.new %in% c('C1QC','APOC1','C1QA','C1QB','GPNMB')]

length(MM.3clusters.new)

length(unique(MM.3clusters.new))

MM.3clusters.new <- c('SELL','CXCL10',MM.3clusters.new,'CD206','LGALS2')

MM.3clusters.new



tiff("MM_3clusters.newMarkers.forPaper.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_3clusters_cluster12.revised, features = MM.3clusters.new, hjust = 0.5, angle = 0, slot = "scale.data") + NoLegend()

dev.off()


write.table(MM.3clusters.new,'MM.3clusters.new.xls',sep='\t',row.names=F,col.names=F,quote=F)


# pie chart for three MM clsuters

table(Idents(MM_3clusters_cluster12.revised))

df <- data.frame(table(Idents(MM_3clusters_cluster12.revised)))

colnames(df) <- c('Cluster','Value')

df 

df$Label <- paste0(Cluster, paste0(" (",round(((df$Value/sum(df$Value))*100),1),"%)"))

library(ggplot2)

p <-  ggplot(df, aes(x = 1, y = Value, fill = Cluster)) + geom_bar(stat = "identity")
p <- p + coord_polar(theta = 'y') + theme_void()
p <- p + geom_text(aes(label = ''), position = position_stack(vjust = 0.5))


tiff("MM_3clusters.Piechart.tiff", width = 3000, height = 3000, units = "px", res = 600)
  
p

dev.off()


q <-  ggplot(df, aes(x = 1, y = Value, fill = Cluster)) + geom_bar(stat = "identity")
q <- q + coord_polar(theta = 'y') + theme_void()
q <- q + geom_text(aes(label = Label), position = position_stack(vjust = 0.5))

tiff("MM_3clusters.Piechart.withPercents.tiff", width = 3000, height = 3000, units = "px", res = 600)
  
q

dev.off()



# use unrevised cluster 12 data to draw heatmap

tiff("MM_3clusters.newMarkers.NotCluster12.revised.forPaper.heatmap.tiff", width = 7000, height = 7000, units = "px", res = 600)
  
DoHeatmap(MM_cluster6_12_23, features = MM.3clusters.new, hjust = 0.5, angle = 0, slot = "scale.data") + NoLegend()

dev.off()


### Caution: Actually did not use the above three lines codes of # use unrevised cluster 12 data to draw heatmap
### Instead used codes of previous section around Line 1252 to 1260to draw heatmap for the 3 MM clusters

##################################################################################################################################################################################
##################################################################################################################################################################################

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/R_NR_separate_analysis')

saveRDS(immune.combined.more.clusters, file='immune.combined.more.clusters.rds')

saveRDS(immune.combined.more.clusters.V2, file='immune.combined.more.clusters.V2.MM2clustersRev.rds')

##################################################################################################################################################################################
##################################################################################################################################################################################

# Some other analyses TBD

################################################################# Using T/NK divisions of previous work ##############################################################

#T_NK_cells.FineClustering <- readRDS('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/T_NK_cells.FineClustering.rds')
#class(T_NK_cells.FineClustering)
#T_NK_cells.FineClustering.Seurat3 = UpdateSeuratObject(object = T_NK_cells.FineClustering)
#table(Idents(T_NK_cells.FineClustering.Seurat3))
#str(T_NK_cells.FineClustering.Seurat3@meta.data)

# following line is from 'E:\Mito_Immunity_onHDD\DefineT_CellPaper_Analysis\SeuratBasedClustering_MarkGenes_Pathway_Analyses\Seurat_V3_Analysis_T_NK\No4_T_NK_clusters_percent_analysis.R'

T_NK_cells_Annotated <- readRDS('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/SingleR_11MainGroups_DefineT_CellPaper_out/No6_T_NK_Seurat_V3/T_NK_cells_Annotated.rds')

dim(T_NK_cells_Annotated@meta.data)

head(T_NK_cells_Annotated@meta.data)


table(Idents(T_NK_cells_Annotated))
 
str(T_NK_cells_Annotated@meta.data)
 
table(T_NK_cells_Annotated@meta.data$T_NK_original_clusterID)

table(T_NK_cells_Annotated@meta.data$T_NK_detailed_clusterID)

length(T_NK_cells_Annotated@meta.data$T_NK_detailed_clusterID)

table(T_NK_cells_Annotated@meta.data$cluster_treatment)


 
table(immune.combined.more.clusters@meta.data$RNA_snn_res.1.V2)
 
class(immune.combined.more.clusters@meta.data)

dim(immune.combined.more.clusters@meta.data)

head(immune.combined.more.clusters@meta.data)


immune.combined.more.clusters@meta.data$cellID <- rownames(immune.combined.more.clusters@meta.data)

T_NK_cells_Annotated_shortTable <- T_NK_cells_Annotated@meta.data[,c('T_NK_original_clusterID','T_NK_detailed_clusterID','T_NK_main_group','cell_treatment','cluster_treatment','cellID')]

dim(T_NK_cells_Annotated_shortTable)

head(T_NK_cells_Annotated_shortTable)

immune.combined.more.clusters@meta.data <- merge(immune.combined.more.clusters@meta.data,T_NK_cells_Annotated_shortTable,by='cellID',all.x=T)


immune.combined.more.clusters@meta.data$IncludeFormerTNKids <- immune.combined.more.clusters@meta.data$T_NK_detailed_clusterID

for (i in 1:nrow(immune.combined.more.clusters@meta.data)) {

	if (is.na(immune.combined.more.clusters@meta.data$IncludeFormerTNKids[i])) {
	
		immune.combined.more.clusters@meta.data$IncludeFormerTNKids[i] <- immune.combined.more.clusters@meta.data$RNA_snn_res.1.V2[i]
	
	}

}

dim(immune.combined.more.clusters@meta.data)

head(immune.combined.more.clusters@meta.data)

table(immune.combined.more.clusters@meta.data$IncludeFormerTNKids)


immune.combined.more.clusters@meta.data[immune.combined.more.clusters@meta.data$IncludeFormerTNKids == '1',]
 


######## Heatmap for T_NK_cells_Annotated_noTgd


T_NK_cells_Annotated_noTgd <- T_NK_cells_Annotated

Idents(T_NK_cells_Annotated_noTgd) <- 'T_NK_detailed_clusterID'

table(Idents(T_NK_cells_Annotated_noTgd))

str(T_NK_cells_Annotated_noTgd@meta.data)

wanted_T_NK <- c(paste0('CD4_c',1:3), paste0('CD8_c',1:9), 'NK_c')

wanted_T_NK


T_NK_cells_Annotated_noTgd <- subset(T_NK_cells_Annotated_noTgd, idents=wanted_T_NK)

str(T_NK_cells_Annotated_noTgd@meta.data)


Idents(T_NK_cells_Annotated_noTgd) <- factor(Idents(T_NK_cells_Annotated_noTgd), levels=wanted_T_NK)

table(Idents(T_NK_cells_Annotated_noTgd))


T_NK_cells_Annotated_noTgd.markers <- FindAllMarkers(object = T_NK_cells_Annotated_noTgd, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
    
library(dplyr)
 

top10 <- T_NK_cells_Annotated_noTgd.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff("T_NK_cells_Annotated_noTgd_Heatmap.tiff", width = 10000, height = 12000, units = "px", res = 600)
  
DoHeatmap(T_NK_cells_Annotated_noTgd, features = top10$gene, hjust = 0.3, angle = 45) + NoLegend()

dev.off()
 


################################################################# Find much more clusters than 23 clusters from 717 ##############################################################

# Want much more clusters

str(immune.combined.more.clusters@meta.data)


# Run the standard workflow for visualization and clustering

immune.combined.more.clusters.V3 <- immune.combined.more.clusters

immune.combined.more.clusters.V3 <- ScaleData(object = immune.combined.more.clusters.V3, verbose = FALSE)

immune.combined.more.clusters.V3 <- RunPCA(object = immune.combined.more.clusters.V3, npcs = 30, verbose = FALSE)


# t-SNE and Clustering

immune.combined.more.clusters.V3 <- RunUMAP(object = immune.combined.more.clusters.V3, reduction = "pca", dims = 1:20)

immune.combined.more.clusters.V3 <- FindNeighbors(object = immune.combined.more.clusters.V3, reduction = "pca", dims = 1:20)

immune.combined.more.clusters.V3 <- FindClusters(immune.combined.more.clusters.V3, resolution = 1.2)

str(immune.combined.more.clusters.V3@meta.data)


# Visualization

new.ids <- as.character(1:23)

names(new.ids) <- as.character(0:22)

new.ids

immune.combined.more.clusters.V3 <- RenameIdents(immune.combined.more.clusters.V3, new.ids)

table(Idents(immune.combined.more.clusters.V3))

immune.combined.more.clusters.V3[['RNA_snn_res.1.V2']] <- Idents(immune.combined.more.clusters.V3)

colnames(immune.combined.more.clusters.V3@meta.data)[7] <- 'integrated_snn_res.0.5.V2'

str(immune.combined.more.clusters.V3@meta.data)



tiff(file = "FigM11_ImmuneCells_Clustering_More_V2.tiff", width = 12000, height = 7000, units = "px", res = 600)

p1 <- DimPlot(object = immune.combined.more.clusters.V3, reduction = "umap", group.by = "Pheno")
p2 <- DimPlot(object = immune.combined.more.clusters.V3, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

dev.off()


##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################
##################################################################################################################################################################################


tiff(file = "Fig3_ImmuneCells_FeaturePlot_SeveralMarkers.tiff", width = 12000, height = 7000, units = "px", res = 600)

FeaturePlot(object = immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", 
    "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")

dev.off()


# assign to original cluster labels of Cell Paper

immune.combined.OriginalCluster <- immune.combined

for (i in 1:11) {

immune.combined.OriginalCluster <- SetIdent(object = immune.combined.OriginalCluster, cells = clustered.cells[[i]], value=names(clustered.cells)[i])

}

str(immune.combined.OriginalCluster)

str(immune.combined.OriginalCluster@active.ident)

levels(immune.combined.OriginalCluster)


############################## Follow "Seurat2_Using sctransform in Seurat.pdf" ################################ 

library(sctransform)

# store mitochondrial percentage in object meta data

CD45overall <- CreateSeuratObject(counts = CD45ImmuneCells_tpm)

CD45overall <- PercentageFeatureSet(object = CD45overall, pattern = "^MT-", col.name = "percent.mt")

str(CD45overall)


# run sctransform

CD45overall <- SCTransform(object = CD45overall, vars.to.regress = "percent.mt", verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering

CD45overall <- RunPCA(object = CD45overall, verbose = FALSE)

CD45overall <- RunUMAP(object = CD45overall, dims = 1:30, verbose = FALSE)

CD45overall <- FindNeighbors(object = CD45overall, dims = 1:30, verbose = FALSE)

CD45overall <- FindClusters(object = CD45overall, verbose = FALSE)

tiff("Fig4_CD45overall_Dimplot.tiff", width=7000, height=6000, units = "px", res = 800)

DimPlot(object = CD45overall, label = TRUE) + NoLegend()

dev.off()


# These are now standard steps in the Seurat workflow for visualization and clustering Visualize
# canonical marker genes as violin plots.

tiff("Fig5_VlnPlot_CD45overall_canonical_markers.tiff", width=12000, height=6000, units = "px", res = 800)

VlnPlot(object = CD45overall, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", 
    "CD3D"), pt.size = 0.2, ncol = 4)
    
dev.off()


# Visualize canonical marker genes on the sctransform embedding.

tiff("Fig6_FeaturePlot_CD45overall_canonical_markers.tiff", width=12000, height=6000, units = "px", res = 800)

FeaturePlot(object = CD45overall, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), pt.size = 0.2, 
    ncol = 3)

dev.off()


tiff("Fig7_FeaturePlot_CD45overall_canonical_markers_II.tiff", width=12000, height=6000, units = "px", res = 800)

FeaturePlot(object = CD45overall, features = c("CD3D", "ISG15", "TCL1A", "FCER2", "XCL1", "FCGR3A"), pt.size = 0.2, 
    ncol = 3)
    
dev.off()

# Examine and visualize PCA results a few different ways

print(x = CD45overall[["pca"]], dims = 1:5, nfeatures = 5)


tiff("Fig8_VizDimLoadings_CD45overall.tiff")

VizDimLoadings(object = CD45overall, dims = 1:2, reduction = "pca")

dev.off()


tiff("Fig9_CD45overall_Dimplot_pcaReduction.tiff")

DimPlot(object = CD45overall, reduction = "pca")

dev.off()


tiff("Fig10_CD45overall_PC_heatmap.tiff")

DimHeatmap(object = CD45overall, dims = 1, cells = 500, balanced = TRUE)

dev.off()


tiff("Fig11_CD45overall_PC_heatmap_15dim.tiff")

DimHeatmap(object = CD45overall, dims = 1:15, cells = 500, balanced = TRUE)

dev.off()


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

CD45overall <- JackStraw(object = CD45overall, num.replicate = 100)

CD45overall <- ScoreJackStraw(object = CD45overall, dims = 1:20)


tiff("Fig12_CD45overall_JackStrawPlot.tiff")

JackStrawPlot(object = CD45overall, dims = 1:20)

dev.off()


tiff("Fig13_CD45overall_ElbowPlot.tiff")

ElbowPlot(object = CD45overall)

dev.off()


CD45overall <- FindNeighbors(object = CD45overall, dims = 1:20)

CD45overall <- FindClusters(object = CD45overall, resolution = 0.5)


# Look at cluster IDs of the first 5 cells
head(x = Idents(object = CD45overall), 5)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')

CD45overall <- RunUMAP(object = CD45overall, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

tiff("Fig14_CD45overall_Dimplot_20dims.tiff", width=7000, height=6000, units = "px", res = 800)

DimPlot(object = CD45overall, label = TRUE, reduction = "umap")

dev.off()

saveRDS(CD45overall, file = "CD45overall_tutorial.rds")


# find all markers of cluster 1

cluster1.markers <- FindMarkers(object = CD45overall, ident.1 = 1, min.pct = 0.25)

head(x = cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3

cluster5.markers <- FindMarkers(object = CD45overall, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

head(x = cluster5.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive ones

CD45overall.markers <- FindAllMarkers(object = CD45overall, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

CD45overall.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(object = CD45overall, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
    
    
tiff("Fig15_VlnPlot_2genes.tiff",width=14000, height=8000, units = "px", res = 800)
    
VlnPlot(object = CD45overall, features = c("MS4A1", "CD79A"))

dev.off()


# you can plot raw counts as well

tiff("Fig16_VlnPlot_rawCountsData.tiff",width=14000, height=8000, units = "px", res = 800)
  
VlnPlot(object = CD45overall, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

dev.off()

tiff("Fig17_FeaturePlot_someMarkers.tiff",width=14000, height=8000, units = "px", res = 800)
 
FeaturePlot(object = CD45overall, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
    
dev.off()    
    

tiff("Fig18_heatmap_CD45overall.tiff",width=18000, height=12000, units = "px", res = 800)
 
top20 <- CD45overall.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

DoHeatmap(object = CD45overall, features = top20$gene)

dev.off()    
  

#new.cluster.ids <- c("Memory CD4 T", "Naive CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
#   "NK", "DC", "Mk")
#names(x = new.cluster.ids) <- levels(x = CD45overall)
#CD45overall <- RenameIdents(object = CD45overall, new.cluster.ids)
#DimPlot(object = CD45overall, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

str(CD45overall)

saveRDS(CD45overall, file = "CD45overall_final.rds")


# assign to original cluster labels of Cell Paper

CD45overall.OriginalCluster <- CD45overall

for (i in 1:11) {

CD45overall.OriginalCluster <- SetIdent(object = CD45overall.OriginalCluster, cells = clustered.cells[[i]], value=names(clustered.cells)[i])

}

str(CD45overall.OriginalCluster)

summary(CD45overall.OriginalCluster)

str(CD45overall.OriginalCluster@active.ident)

levels(CD45overall.OriginalCluster)

 
tiff("Fig19_CD45overall_OriginalCluster_Dimplot.tiff", width=12000, height=8000, units = "px", res = 800)
 
DimPlot(object = CD45overall.OriginalCluster, reduction = "umap", label = TRUE, pt.size = 0.5)  

dev.off()


dim(CD45overall.OriginalCluster@assays$RNA@data)

head(CD45overall.OriginalCluster@assays$RNA@data[,1:10])

tail(CD45overall.OriginalCluster@assays$RNA@data[,1:10])

#str(CD45overall.OriginalCluster@reductions$tsne$umap@cell.embeddings)


saveRDS(CD45overall.OriginalCluster, file = "CD45overall.OriginalCluster.rds")


# Assign new cell types based on SingleR results

SingleR.DefineT.cellID <- read.table('SingleR.DefineT.cellID.txt', sep='\t', header=T)

dim(SingleR.DefineT.cellID)

head(SingleR.DefineT.cellID)

tail(SingleR.DefineT.cellID)

sort(unique(SingleR.DefineT.cellID$type))


SingleR.cluster.names <- as.character(sort(unique(SingleR.DefineT.cellID$type)))

SingleR.cluster.cells <- list()

for (i in 1:length(SingleR.cluster.names)) {

		SingleR.cluster.cells[[i]] <- as.character(SingleR.DefineT.cellID$cell[SingleR.DefineT.cellID$type == SingleR.cluster.names[i]])
		
		names(SingleR.cluster.cells)[i] <- SingleR.cluster.names[i]
		
		print(names(SingleR.cluster.cells)[i])
		
		print(length(SingleR.cluster.cells[[i]]))

}

length(SingleR.cluster.cells)

lapply(SingleR.cluster.cells,length)

sub.immune.cells.enough.number <- lapply(SingleR.cluster.cells,length)[lapply(SingleR.cluster.cells,length)>=10]

length(sub.immune.cells.enough.number)

sub.immune.cells.enough.number

names(sub.immune.cells.enough.number)

sub.immune.cells.enough.number.V2 <- as.data.frame(sub.immune.cells.enough.number)

write.table(sub.immune.cells.enough.number.V2,"sub.immune.cells.enough.number.V2.xls",sep='\t',row.names=F,quote=F)



CD45overall.SingleRcluster <- CD45overall

for (i in 1:length(SingleR.cluster.cells)) {

CD45overall.SingleRcluster <- SetIdent(object = CD45overall.SingleRcluster, cells = SingleR.cluster.cells[[i]], value=names(SingleR.cluster.cells)[i])

}

str(CD45overall.SingleRcluster)

summary(CD45overall.SingleRcluster)

str(CD45overall.SingleRcluster@active.ident)

levels(CD45overall.SingleRcluster)

 
tiff("Fig20_CD45overall_SingleRcluster_Dimplot_All.tiff", width=16000, height=14000, units = "px", res = 600)
 
DimPlot(object = CD45overall.SingleRcluster, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

dev.off()


# plot those with significant # of cells

sigNumberClusters <- names(sub.immune.cells.enough.number)

sigNumberClusters <- sigNumberClusters[sigNumberClusters!='B_cell']

sigNumberClusters


sigNumberClusters.cells <- as.character(SingleR.DefineT.cellID$cell[SingleR.DefineT.cellID$type %in% sigNumberClusters])

length(sigNumberClusters.cells)

head(sigNumberClusters.cells)



CD45overall.SingleRcluster.35BigGroups <- SubsetData(object = CD45overall.SingleRcluster, cells = sigNumberClusters.cells)  

SingleR.cluster.cells.V2 <- SingleR.cluster.cells[names(SingleR.cluster.cells) %in% sigNumberClusters]

lapply(SingleR.cluster.cells.V2,length)

SingleR.cluster.cells.V2[[1]]

length(SingleR.cluster.cells.V2)

names(SingleR.cluster.cells.V2)

for (i in 1:length(SingleR.cluster.cells.V2)) {

		CD45overall.SingleRcluster.35BigGroups <- SetIdent(object = CD45overall.SingleRcluster.35BigGroups, 
		cells = SingleR.cluster.cells.V2[[i]], value=names(SingleR.cluster.cells.V2)[i])

}

str(CD45overall.SingleRcluster.35BigGroups)

summary(CD45overall.SingleRcluster.35BigGroups)

str(CD45overall.SingleRcluster.35BigGroups@active.ident)

levels(CD45overall.SingleRcluster.35BigGroups)

saveRDS(CD45overall.SingleRcluster.35BigGroups, file = "CD45overall.SingleRcluster.35BigGroups.rds")



tiff("Fig21_CD45overall_SingleRcluster_35BigGroups_Dimplot.tiff", width=13000, height=11000, units = "px", res = 600)
 
DimPlot(object = CD45overall.SingleRcluster.35BigGroups, reduction = "umap", label = TRUE, pt.size = 0.3)

dev.off()

 
################################################ DE analysis of each of the immune subpopulations

# first assign pheno information to each cell

unique(DefineT_PhenoInfo_V2$Response)

ResponderCells <- as.character(DefineT_PhenoInfo_V2$Cell.Name[DefineT_PhenoInfo_V2$Response=="Responder"])

length(ResponderCells)

head(ResponderCells)

 
NonResponderCells <- as.character(DefineT_PhenoInfo_V2$Cell.Name[DefineT_PhenoInfo_V2$Response=="Non-responder"])

length(NonResponderCells)

head(NonResponderCells)

 
modified.names.35population <- gsub(':','_',names(SingleR.cluster.cells.V2))

modified.names.35population <- gsub('-','_',modified.names.35population)

modified.names.35population <- gsub('/','_',modified.names.35population)

modified.names.35population


EachOf35ImmunePopulation <- list()


for (i in 1:length(SingleR.cluster.cells.V2)) {

currentResponderCells <- intersect(ResponderCells, SingleR.cluster.cells.V2[[i]])

currentNonResponderCells <- intersect(NonResponderCells, SingleR.cluster.cells.V2[[i]])

sum(length(currentResponderCells),length(currentNonResponderCells))==length(SingleR.cluster.cells.V2[[i]])


current.matrix <- CD45ImmuneCells_tpm[,c(currentResponderCells,currentNonResponderCells)]

EachOf35ImmunePopulation[[i]] <- apply(as.matrix.noquote(current.matrix),2,as.numeric) 

rownames(EachOf35ImmunePopulation[[i]]) <- rownames(current.matrix)


new.names <- c(rep("Responder", length(currentResponderCells)), rep("NonResponder", length(currentNonResponderCells)))

EachOf35ImmunePopulation[[i]] <- EachOf35ImmunePopulation[[i]][rowSums(EachOf35ImmunePopulation[[i]]) > 3,]

colnames(EachOf35ImmunePopulation[[i]]) <- new.names



current.output.file <- EachOf35ImmunePopulation[[i]]

current.output.file <- as.data.frame(current.output.file)

current.output.file$GeneID <- rownames(current.output.file)


ncolOfDF <- ncol(current.output.file)

nrowOfDF <- nrow(current.output.file)

current.output.file[nrowOfDF+1,] <- colnames(current.output.file)


current.output.file <- current.output.file[,c(ncolOfDF, (1:(ncolOfDF-1)))]

current.output.file <- current.output.file[c(nrowOfDF+1, (1:nrowOfDF)),]


names(EachOf35ImmunePopulation)[i] <- paste0('No',i,'_', modified.names.35population[i]) 
 
out.name <- paste0('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput/',names(EachOf35ImmunePopulation)[i]) 

write.table(current.output.file,out.name,sep='\t',row.names=F, col.names=F,quote=F)

}

str(EachOf35ImmunePopulation)

names(EachOf35ImmunePopulation)

D3Einput.names.35Populations <- names(EachOf35ImmunePopulation)

save(D3Einput.names.35Populations,file="D3Einput.names.35Populations.RData")


saveRDS(EachOf35ImmunePopulation,"EachOf35ImmunePopulation.matrix.rds")

EachOf35ImmunePopulation.ForHeatmap <- readRDS("EachOf35ImmunePopulation.matrix.rds")

str(EachOf35ImmunePopulation.ForHeatmap)

names(EachOf35ImmunePopulation.ForHeatmap)


############################################################## output non.log2 data for D3E ##########################################################################


setwd("/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput.non.log2")


EachOf35ImmunePopulation.non.log2.ForHeatmap <- list()


for (i in 1:length(EachOf35ImmunePopulation.ForHeatmap)) {

EachOf35ImmunePopulation.non.log2.ForHeatmap[[i]] <- 2^EachOf35ImmunePopulation.ForHeatmap[[i]]-1

names(EachOf35ImmunePopulation.non.log2.ForHeatmap)[i] <- names(EachOf35ImmunePopulation.ForHeatmap)[i]
		
current.output.non.log2.file <- EachOf35ImmunePopulation.non.log2.ForHeatmap[[i]]

current.output.non.log2.file <- as.data.frame(current.output.non.log2.file)

current.output.non.log2.file$GeneID <- rownames(current.output.non.log2.file)


ncolOfDF <- ncol(current.output.non.log2.file)

nrowOfDF <- nrow(current.output.non.log2.file)

current.output.non.log2.file[nrowOfDF+1,] <- colnames(current.output.non.log2.file)


current.output.non.log2.file <- current.output.non.log2.file[,c(ncolOfDF, (1:(ncolOfDF-1)))]

current.output.non.log2.file <- current.output.non.log2.file[c(nrowOfDF+1, (1:nrowOfDF)),]

		
out.name2 <- names(EachOf35ImmunePopulation.non.log2.ForHeatmap)[i] 

write.table(current.output.non.log2.file,out.name2,sep='\t',row.names=F, col.names=F,quote=F)

}

lapply(EachOf35ImmunePopulation.non.log2.ForHeatmap,nrow)

lapply(EachOf35ImmunePopulation.non.log2.ForHeatmap,dim)


###################################################################################### Heatmap on all DE genes (mode0) calculated on Cluster ######################################################################################

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput/Mode0_out')

(D3E.mode0.outfiles <- list.files())


(D3E.mode0.outfiles.names <- gsub('\\.out$','',D3E.mode0.outfiles))

D3E.mode0.29celltypes <- list()

header.names <- c("GeneID", "a1", "b1", "g1", "GOF1", "a2", "b2", "g2", "GOF2", "s1", 
"f1", "d1", "s2", "f2", "d2", "Rs", "Rf", "Rd", "p.value", "mu1", "cv1", "mu2", "cv2")


for (i in 1:length(D3E.mode0.outfiles)) {

		D3E.mode0.29celltypes[[i]] <- read.table(D3E.mode0.outfiles[[i]],sep='\t',header=F)
		
		colnames(D3E.mode0.29celltypes[[i]]) <- header.names
		
		D3E.mode0.29celltypes[[i]]$FoldChange <- D3E.mode0.29celltypes[[i]]$mu1/D3E.mode0.29celltypes[[i]]$mu2
		
		names(D3E.mode0.29celltypes)[i] <- D3E.mode0.outfiles.names[i]
		
		print(names(D3E.mode0.29celltypes)[i])
		
		print(dim(D3E.mode0.29celltypes[[i]]))
		
		print(head(D3E.mode0.29celltypes[[i]]))
		
}


str(D3E.mode0.29celltypes)

length(D3E.mode0.29celltypes)


D3E.mode0.29celltypes.sig.genes <- list()

for (i in 1:length(D3E.mode0.29celltypes)) {

		crit1 <- (D3E.mode0.29celltypes[[i]]$p.value < 0.05)
		
		crit2 <- (D3E.mode0.29celltypes[[i]]$cv1 < D3E.mode0.29celltypes[[i]]$mu1)
		
		crit3 <- (D3E.mode0.29celltypes[[i]]$cv2 < D3E.mode0.29celltypes[[i]]$mu2)
		
		crit4 <- (D3E.mode0.29celltypes[[i]]$FoldChange < 1/3) | (D3E.mode0.29celltypes[[i]]$FoldChange > 3)
		
		#D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit2 & crit3 & crit4,]
		
		#D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit4,]
		
		D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes[[i]][crit1,]
		
		D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes.sig.genes[[i]][!is.na(D3E.mode0.29celltypes.sig.genes[[i]]$GeneID),]
		
		
		names(D3E.mode0.29celltypes.sig.genes)[i] <- names(D3E.mode0.29celltypes)[i]
		
		
		crit5 <- (D3E.mode0.29celltypes.sig.genes[[i]]$mu1 < 1)
		
		crit6 <- (D3E.mode0.29celltypes.sig.genes[[i]]$mu2 < 1)
		
		
		#D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes.sig.genes[[i]][!(crit5 & crit6),]
		
		D3E.mode0.29celltypes.sig.genes[[i]] <- D3E.mode0.29celltypes.sig.genes[[i]][order(D3E.mode0.29celltypes.sig.genes[[i]]$FoldChange, D3E.mode0.29celltypes.sig.genes[[i]]$p.value),]
		
		print(names(D3E.mode0.29celltypes.sig.genes)[i])
		
		print(dim(D3E.mode0.29celltypes.sig.genes[[i]]))
		
		print(head(D3E.mode0.29celltypes.sig.genes[[i]]))

}

lapply(D3E.mode0.29celltypes.sig.genes,nrow) # mu filtering: mu1 < 1 and mu2 < 1 not used yet

lapply(D3E.mode0.29celltypes.sig.genes,nrow)[lapply(D3E.mode0.29celltypes.sig.genes,nrow)<=1]

names(lapply(D3E.mode0.29celltypes.sig.genes,nrow)[lapply(D3E.mode0.29celltypes.sig.genes,nrow)<=1])


D3E.mode0.29celltypes.sig.genes <- D3E.mode0.29celltypes.sig.genes[!names(D3E.mode0.29celltypes.sig.genes) %in% names(lapply(D3E.mode0.29celltypes.sig.genes,nrow)[lapply(D3E.mode0.29celltypes.sig.genes,nrow)<=1])]

length(D3E.mode0.29celltypes.sig.genes) # current 29


##################################################################### Heatmap step #####################################################################################

# Heatmaps for 29 celltypes based on sig genes identified from mode0

# frequency filtering sub-routine

row_sub = function(x) {

		selection <- apply(x, 1, function(row) sum(row > 3)<1)
		
		x2 <- rownames(x)[selection]
		
		x2

}


D3E.mode0.29celltypes.sig.genes <- D3E.mode0.29celltypes.sig.genes[!names(D3E.mode0.29celltypes.sig.genes) %in% c("No9_DC_monocyte_derived_mature","No16_Monocyte_anti_FcgRIIB")]

names(D3E.mode0.29celltypes.sig.genes)

EachOf29ImmunePopulation.ForHeatmap <- EachOf35ImmunePopulation.ForHeatmap[names(EachOf35ImmunePopulation.ForHeatmap) %in% names(D3E.mode0.29celltypes.sig.genes)]

length(EachOf29ImmunePopulation.ForHeatmap)

all(sort(names(D3E.mode0.29celltypes.sig.genes))==sort(names(EachOf29ImmunePopulation.ForHeatmap)))



EachOf29ImmunePopulation.heatmap.siggenes <- list()


# filter significant genes based on previous p < 0.05 and frequency filtering based on sub-routine: row_sub 

setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput/Mode0_outProcessed')

for (i in c(1:length(EachOf29ImmunePopulation.ForHeatmap))) {

labels.current <- colnames(EachOf29ImmunePopulation.ForHeatmap[i][[1]])

length(labels.current[labels.current=="Responder"])

length(labels.current[labels.current=="NonResponder"])

names(EachOf29ImmunePopulation.ForHeatmap)[i]

head(EachOf29ImmunePopulation.ForHeatmap[i][[1]])

dim(EachOf29ImmunePopulation.ForHeatmap[i][[1]])


groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

current.cell.type <- names(EachOf29ImmunePopulation.ForHeatmap)[i]


matrix.current <- EachOf29ImmunePopulation.ForHeatmap[i][[1]]

genes.use <- as.character(D3E.mode0.29celltypes.sig.genes[current.cell.type][[1]]$GeneID)

genes.use <- genes.use[!is.na(genes.use)]

print(length(genes.use))

print(head(genes.use))


matrix.current <- matrix.current[genes.use,]

dim(matrix.current)

head(matrix.current)


R.col.end <- length(labels.current[labels.current=="Responder"])

NR.col.start <- length(labels.current[labels.current=="Responder"])+1


matrix.current.Responder <- matrix.current[,1:R.col.end]

dim(matrix.current.Responder)


matrix.current.NonResponder <- matrix.current[,NR.col.start:ncol(matrix.current)]

dim(matrix.current.NonResponder)



NonResponderOnly.gene <- row_sub(matrix.current.Responder)

length(NonResponderOnly.gene)

head(NonResponderOnly.gene)


ResponderOnly.gene <- row_sub(matrix.current.NonResponder)

length(ResponderOnly.gene)

head(ResponderOnly.gene)


genes.use2 <- unique(c(NonResponderOnly.gene,ResponderOnly.gene))

EachOf29ImmunePopulation.heatmap.siggenes[[i]] <- matrix.current[genes.use2,]

print(dim(EachOf29ImmunePopulation.heatmap.siggenes[[i]]))

print(head(EachOf29ImmunePopulation.heatmap.siggenes[[i]]))

names(EachOf29ImmunePopulation.heatmap.siggenes)[i] <- names(EachOf29ImmunePopulation.ForHeatmap)[i]


column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(EachOf29ImmunePopulation.ForHeatmap)[i],".sig.genes.tiff")

tiff(file = Heatmap.file.name, width = 8000, height = 8000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(EachOf29ImmunePopulation.heatmap.siggenes[[i]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.6,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

}


lapply(EachOf29ImmunePopulation.heatmap.siggenes,nrow)

length(EachOf29ImmunePopulation.heatmap.siggenes)

names(EachOf29ImmunePopulation.heatmap.siggenes)


EachOf29ImmunePopulation.heatmap.siggenes[[18]] <- 
EachOf29ImmunePopulation.heatmap.siggenes[[18]][rownames(EachOf29ImmunePopulation.heatmap.siggenes[[18]]) != "PDXDC2P",]


saveRDS(EachOf29ImmunePopulation.heatmap.siggenes,"Individual27ImmunePopulation.heatmap.siggenes.1stRound.rds")

#EachOf29ImmunePopulation.heatmap.siggenes <- readRDS("Individual27ImmunePopulation.heatmap.siggenes.1stRound.rds")



## Among above, redraw No25

labels.current <- colnames(EachOf29ImmunePopulation.heatmap.siggenes[[18]])

groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(EachOf29ImmunePopulation.heatmap.siggenes)[18],".sig.genes.tiff")

Heatmap.file.name


tiff(file = Heatmap.file.name, width = 6000, height = 12000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(EachOf29ImmunePopulation.heatmap.siggenes[[18]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.2,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

dim(EachOf29ImmunePopulation.heatmap.siggenes[[18]])


################################################################################### further refine and draw heatmap: Round 2############################################################

FurtherMapping.siggenes <- EachOf29ImmunePopulation.heatmap.siggenes[!names(EachOf29ImmunePopulation.heatmap.siggenes) %in% c("No25_Pre_B_cell_CD34_")]

length(FurtherMapping.siggenes)

names(FurtherMapping.siggenes)


D3E.mode0.29celltypes.sig.genes.V2 <- list()

for (i in 1:length(D3E.mode0.29celltypes)) {

		crit1 <- (D3E.mode0.29celltypes[[i]]$p.value < 0.05)
		
		crit2 <- (D3E.mode0.29celltypes[[i]]$cv1 < D3E.mode0.29celltypes[[i]]$mu1)
		
		crit3 <- (D3E.mode0.29celltypes[[i]]$cv2 < D3E.mode0.29celltypes[[i]]$mu2)
		
		crit4 <- (D3E.mode0.29celltypes[[i]]$FoldChange < 1/3) | (D3E.mode0.29celltypes[[i]]$FoldChange > 3)
		
		#D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit2 & crit3 & crit4,]
		
		#D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit4,]
		
		D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes[[i]][crit1,]
		
		D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes.sig.genes.V2[[i]][!is.na(D3E.mode0.29celltypes.sig.genes.V2[[i]]$GeneID),]
		
		
		names(D3E.mode0.29celltypes.sig.genes.V2)[i] <- names(D3E.mode0.29celltypes)[i]
		
		
		crit5 <- (D3E.mode0.29celltypes.sig.genes.V2[[i]]$mu1 < 1)
		
		crit6 <- (D3E.mode0.29celltypes.sig.genes.V2[[i]]$mu2 < 1)
		
		
		D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes.sig.genes.V2[[i]][!(crit5 & crit6),]
		
		D3E.mode0.29celltypes.sig.genes.V2[[i]] <- D3E.mode0.29celltypes.sig.genes.V2[[i]][order(D3E.mode0.29celltypes.sig.genes.V2[[i]]$FoldChange, D3E.mode0.29celltypes.sig.genes.V2[[i]]$p.value),]
		
		print(names(D3E.mode0.29celltypes.sig.genes.V2)[i])
		
		print(dim(D3E.mode0.29celltypes.sig.genes.V2[[i]]))
		
		print(head(D3E.mode0.29celltypes.sig.genes.V2[[i]]))

}

lapply(D3E.mode0.29celltypes.sig.genes.V2,nrow)

lapply(D3E.mode0.29celltypes.sig.genes.V2,nrow)[lapply(D3E.mode0.29celltypes.sig.genes.V2,nrow)<=1]

names(lapply(D3E.mode0.29celltypes.sig.genes.V2,nrow)[lapply(D3E.mode0.29celltypes.sig.genes.V2,nrow)<=1])


D3E.mode0.29celltypes.sig.genes.V2 <- D3E.mode0.29celltypes.sig.genes.V2[names(D3E.mode0.29celltypes.sig.genes.V2) %in% names(FurtherMapping.siggenes)]

length(D3E.mode0.29celltypes.sig.genes.V2)

all(sort(names(D3E.mode0.29celltypes.sig.genes.V2))==sort(names(FurtherMapping.siggenes)))




FurtherMapping.siggenes.V2 <- list()

FurtherMapping.siggenes.SelectForLoop <- FurtherMapping.siggenes[!names(FurtherMapping.siggenes) %in% c("No3_B_cell_Memory","No4_B_cell_Naive",
"No5_B_cell_Plasma_cell","No22_NK_cell","No27_T_cell_CD4+_central_memory","No28_T_cell_CD4+_effector_memory","No30_T_cell_CD8+",
"No31_T_cell_CD8+_Central_memory","No33_T_cell_CD8+_effector_memory_RA","No35_T_cell_gamma_delta")] 

length(FurtherMapping.siggenes.SelectForLoop)


FurtherMapping.siggenes.NotSelectForLoop <- FurtherMapping.siggenes[names(FurtherMapping.siggenes) %in% c("No3_B_cell_Memory","No4_B_cell_Naive",
"No5_B_cell_Plasma_cell","No22_NK_cell","No27_T_cell_CD4+_central_memory","No28_T_cell_CD4+_effector_memory","No30_T_cell_CD8+",
"No31_T_cell_CD8+_Central_memory","No33_T_cell_CD8+_effector_memory_RA","No35_T_cell_gamma_delta")] 

length(FurtherMapping.siggenes.NotSelectForLoop)


setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput/Mode0_outProcessed.V2')

for (i in c(1:length(FurtherMapping.siggenes.SelectForLoop))) {

labels.current <- colnames(FurtherMapping.siggenes.SelectForLoop[i][[1]])

length(labels.current[labels.current=="Responder"])

length(labels.current[labels.current=="NonResponder"])

names(FurtherMapping.siggenes.SelectForLoop)[i]

head(FurtherMapping.siggenes.SelectForLoop[i][[1]])


groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

current.cell.type <- names(FurtherMapping.siggenes.SelectForLoop)[i]


matrix.current <- FurtherMapping.siggenes.SelectForLoop[i][[1]]

genes.use <- intersect(rownames(FurtherMapping.siggenes.SelectForLoop[i][[1]]),D3E.mode0.29celltypes.sig.genes.V2[current.cell.type][[1]]$GeneID)

genes.use <- genes.use[!is.na(genes.use)]

print(length(genes.use))

print(head(genes.use))


matrix.current <- matrix.current[genes.use,]

dim(matrix.current)

head(matrix.current)


FurtherMapping.siggenes.V2[[i]] <- matrix.current

print(dim(FurtherMapping.siggenes.V2[[i]]))

print(head(FurtherMapping.siggenes.V2[[i]]))

names(FurtherMapping.siggenes.V2)[i] <- names(FurtherMapping.siggenes.SelectForLoop)[i]


column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(FurtherMapping.siggenes.SelectForLoop)[i],".sig.genes.tiff")

tiff(file = Heatmap.file.name, width = 8000, height = 8000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(FurtherMapping.siggenes.V2[[i]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.6,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

}


lapply(FurtherMapping.siggenes.V2,nrow)

length(FurtherMapping.siggenes.V2)

names(FurtherMapping.siggenes.V2)


lapply(FurtherMapping.siggenes.V2,dim)


# Delete one inconsistent gene in No18

FurtherMapping.siggenes.V2[[8]] <- 
FurtherMapping.siggenes.V2[[8]][rownames(FurtherMapping.siggenes.V2[[8]]) != "RPL3P2",]


saveRDS(FurtherMapping.siggenes.V2,"Individual16ImmunePopulation.heatmap.siggenes.2ndRound.rds")

dim(FurtherMapping.siggenes.V2$No24_NK_cell_IL2)

head(FurtherMapping.siggenes.V2$No24_NK_cell_IL2)

dim(FurtherMapping.siggenes.V2[[14]])

head(FurtherMapping.siggenes.V2[[14]])


## Among above, redraw No2, No18

# No2

labels.current <- colnames(FurtherMapping.siggenes.V2[[2]])

groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(FurtherMapping.siggenes.V2)[2],".sig.genes.tiff")

Heatmap.file.name


tiff(file = Heatmap.file.name, width = 6000, height = 10000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(FurtherMapping.siggenes.V2[[2]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.2,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

dim(FurtherMapping.siggenes.V2[[2]])



# No18

labels.current <- colnames(FurtherMapping.siggenes.V2[[8]])

groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(FurtherMapping.siggenes.V2)[8],".sig.genes.tiff")

Heatmap.file.name


tiff(file = Heatmap.file.name, width = 5000, height = 4000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(FurtherMapping.siggenes.V2[[8]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.6,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

dim(FurtherMapping.siggenes.V2[[8]])



################################################################################### further refine and draw heatmap: Round 3 ############################################################

Round3Mapping.siggenes <- FurtherMapping.siggenes.V2[names(FurtherMapping.siggenes.V2) %in% c("No12_Macrophage_monocyte_derived_M_CSF","No13_Macrophage_monocyte_derived_M_CSF_IFNg",
"No14_Macrophage_monocyte_derived_M_CSF_IFNg_Pam3Cys","No15_Monocyte","No17_Monocyte_CD14+")]

length(Round3Mapping.siggenes)

names(Round3Mapping.siggenes)


D3E.mode0.29celltypes.sig.genes.V3 <- list()

for (i in 1:length(D3E.mode0.29celltypes)) {

		crit1 <- (D3E.mode0.29celltypes[[i]]$p.value < 0.05)
		
		crit2 <- (D3E.mode0.29celltypes[[i]]$cv1 < D3E.mode0.29celltypes[[i]]$mu1)
		
		crit3 <- (D3E.mode0.29celltypes[[i]]$cv2 < D3E.mode0.29celltypes[[i]]$mu2)
		
		crit4 <- (D3E.mode0.29celltypes[[i]]$FoldChange < 1/3) | (D3E.mode0.29celltypes[[i]]$FoldChange > 3)
		
		#D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit2 & crit3 & crit4,]
		
		#D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes[[i]][crit1 & crit4,]
		
		D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes[[i]][crit1,]
		
		D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes.sig.genes.V3[[i]][!is.na(D3E.mode0.29celltypes.sig.genes.V3[[i]]$GeneID),]
		
		
		names(D3E.mode0.29celltypes.sig.genes.V3)[i] <- names(D3E.mode0.29celltypes)[i]
		
		
		crit5 <- (D3E.mode0.29celltypes.sig.genes.V3[[i]]$mu1 < 4)
		
		crit6 <- (D3E.mode0.29celltypes.sig.genes.V3[[i]]$mu2 < 4)
		
		
		D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes.sig.genes.V3[[i]][!(crit5 & crit6),]
		
		D3E.mode0.29celltypes.sig.genes.V3[[i]] <- D3E.mode0.29celltypes.sig.genes.V3[[i]][order(D3E.mode0.29celltypes.sig.genes.V3[[i]]$FoldChange, D3E.mode0.29celltypes.sig.genes.V3[[i]]$p.value),]
		
		print(names(D3E.mode0.29celltypes.sig.genes.V3)[i])
		
		print(dim(D3E.mode0.29celltypes.sig.genes.V3[[i]]))
		
		print(head(D3E.mode0.29celltypes.sig.genes.V3[[i]]))

}

lapply(D3E.mode0.29celltypes.sig.genes.V3,nrow)

lapply(D3E.mode0.29celltypes.sig.genes.V3,nrow)[lapply(D3E.mode0.29celltypes.sig.genes.V3,nrow)<=1]

names(lapply(D3E.mode0.29celltypes.sig.genes.V3,nrow)[lapply(D3E.mode0.29celltypes.sig.genes.V3,nrow)<=1])


D3E.mode0.29celltypes.sig.genes.V3 <- D3E.mode0.29celltypes.sig.genes.V3[names(D3E.mode0.29celltypes.sig.genes.V3) %in% names(Round3Mapping.siggenes)]

length(D3E.mode0.29celltypes.sig.genes.V3)

all(sort(names(D3E.mode0.29celltypes.sig.genes.V3))==sort(names(Round3Mapping.siggenes)))



Round3Mapping.siggenes.V2 <- list()

Round3Mapping.siggenes.SelectForLoop <- Round3Mapping.siggenes



setwd('/rcc/temp_dxiong/DefineT_CellPaper_Analysis/D3Einput/Mode0_outProcessed.V3')

for (i in c(1:length(Round3Mapping.siggenes.SelectForLoop))) {

labels.current <- colnames(Round3Mapping.siggenes.SelectForLoop[i][[1]])

length(labels.current[labels.current=="Responder"])

length(labels.current[labels.current=="NonResponder"])

names(Round3Mapping.siggenes.SelectForLoop)[i]

head(Round3Mapping.siggenes.SelectForLoop[i][[1]])


groups.current <- factor(labels.current,levels=c("Responder","NonResponder"))

current.cell.type <- names(Round3Mapping.siggenes.SelectForLoop)[i]


matrix.current <- Round3Mapping.siggenes.SelectForLoop[i][[1]]

genes.use <- intersect(rownames(Round3Mapping.siggenes.SelectForLoop[i][[1]]),D3E.mode0.29celltypes.sig.genes.V3[current.cell.type][[1]]$GeneID)

genes.use <- genes.use[!is.na(genes.use)]

print(length(genes.use))

print(head(genes.use))


matrix.current <- matrix.current[genes.use,]

dim(matrix.current)

head(matrix.current)


Round3Mapping.siggenes.V2[[i]] <- matrix.current

print(dim(Round3Mapping.siggenes.V2[[i]]))

print(head(Round3Mapping.siggenes.V2[[i]]))

names(Round3Mapping.siggenes.V2)[i] <- names(Round3Mapping.siggenes.SelectForLoop)[i]


column_annotation <- c(rep("red",length(groups.current[groups.current=='Responder'])),
rep("blue",length(groups.current[groups.current=='NonResponder'])))

column_annotation <- as.matrix(column_annotation)

colfunc <- colorRampPalette(c("blue","white", "red"))


Heatmap.file.name <- paste0(names(Round3Mapping.siggenes.SelectForLoop)[i],".sig.genes.tiff")

tiff(file = Heatmap.file.name, width = 8000, height = 8000, units = "px", res = 800)

par(oma=c(1,0,0,6)) # down, left, up, right

heatmap.3(Round3Mapping.siggenes.V2[[i]],col=colfunc(120),scale="row", trace="none",Rowv=FALSE,Colv=FALSE,cexRow=1.6,cexCol=1.2,labCol = "",ColSideColors=column_annotation) # remove 'labCol = "",' to show sample names

dev.off()

}


lapply(Round3Mapping.siggenes.V2,nrow)

length(Round3Mapping.siggenes.V2)

names(Round3Mapping.siggenes.V2)


lapply(Round3Mapping.siggenes.V2,dim)


######################################################################################################## clustering CD4T.subset based on further subclustering OE scores or feature genes






