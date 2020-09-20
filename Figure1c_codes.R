
rm(list=ls())

setwd('C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes/PatientCompare')

cell.type.percent.summary.overall <- readRDS("cell.type.percent.patient.level.rds")

Ordered_23clusters_cellIDs <- readRDS("Ordered_23clusters_cellIDs.rds")

library(ggpubr)

library(coin)



# Plot significant ones with P values displayed

for (i in c(6,9,12,13,14,17,19,21,22)) {

		current.df <- cell.type.percent.summary.overall[cell.type.percent.summary.overall$cell.type == names(Ordered_23clusters_cellIDs)[i],]
				
		print(as.character(unique(current.df$cell.type)))
		
		wt <- wilcox_test(percent ~ pheno, data = current.df, distribution = "exact") 
		
		out.name <- paste0("No",i,"_",names(Ordered_23clusters_cellIDs)[i],".boxplot.tiff")
		
		current.df$percent <- 100*current.df$percent
		
		xlab_text <- paste0(as.character(unique(current.df$cell.type)),": P = ",round(pvalue(wt),4)) 
		

		tiff(file = out.name, width = 300, height = 400, units = "px", res = 120)
		
		# Add jitter points and errors (mean_se)

		print(ggboxplot(current.df, x = "pheno", y = "percent", ylab="(%) of CD45+", label.pos = "out", lab.nb.digits = 1, xlab=xlab_text, lab.size=3,
		add = c("jitter"), width=0.5, color="pheno", palette = c("#00AFBB", "#E7B800"), legend=NULL, shape = "pheno"))

		dev.off()

}

