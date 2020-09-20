rm(list=ls())

####### Must use R x64 3.2.3 or R version 3.4.1 (2017-06-30) -- "Single Candle"

####### o.rocc

library(rocc)

library(ggplot2)

library(Biobase)

library(edgeR)

library(limma)

library("RColorBrewer")

library(biomaRt)

#library(RUVSeq)

#library(EDASeq)

library(cancerclass)

rocdata <- function(grp, pred){
  # Produces x and y co-ordinates for ROC curve plot
  # Arguments: grp - labels classifying subject status
  #            pred - values of each observation
  # Output: List with 2 components:
  #         roc = data.frame with x and y co-ordinates of plot
  #         stats = data.frame containing: area under ROC curve, p value, upper and lower 95% confidence interval
 
  grp <- as.factor(grp)
  if (length(pred) != length(grp)) {
    stop("The number of classifiers must match the number of data points")
  } 
 
  if (length(levels(grp)) != 2) {
    stop("There must only be 2 values for the classifier")
  }
 
  cut <- unique(pred)
  tp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[2])))
  fn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[2])))
  fp <- sapply(cut, function(x) length(which(pred > x & grp == levels(grp)[1])))
  tn <- sapply(cut, function(x) length(which(pred < x & grp == levels(grp)[1])))
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  roc = data.frame(x = fpr, y = tpr)
  roc <- roc[order(roc$x, roc$y),]
 
  i <- 2:nrow(roc)
  auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
 
  pos <- pred[grp == levels(grp)[2]]
  neg <- pred[grp == levels(grp)[1]]
  q1 <- auc/(2-auc)
  q2 <- (2*auc^2)/(1+auc)
  se.auc <- sqrt(((auc * (1 - auc)) + ((length(pos) -1)*(q1 - auc^2)) + ((length(neg) -1)*(q2 - auc^2)))/(length(pos)*length(neg)))
  ci.upper <- auc + (se.auc * 0.96)
  ci.lower <- auc - (se.auc * 0.96)
 
  se.auc.null <- sqrt((1 + length(pos) + length(neg))/(12*length(pos)*length(neg)))
  z <- (auc - 0.5)/se.auc.null
  p <- 2*pnorm(-abs(z))
 
  stats <- data.frame (auc = auc,
                       p.value = p,
                       ci.upper = ci.upper,
                       ci.lower = ci.lower
                       )
 
  return (list(roc = roc, stats = stats))
}


# single ROC plot


rocplot.single.V2 <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC = ",signif(auc, 2), " (95% CI: ", signif(ci.lower, 2), "-", signif(ci.upper, 2), ")", sep=""))
  }
 
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
      geom_line(aes(colour = "")) +
      geom_abline (intercept = 0, slope = 1) +
      theme_bw() +
      scale_x_continuous("False Positive Rate (1-Specificity)") +
      scale_y_continuous("True Positive Rate (Sensitivity)") +
      scale_colour_manual(labels = annotation, values = "#000000") +
      theme(
           plot.title = element_text(face="bold", size=14), 
           axis.text.x = element_text(face="bold", size=18),
           axis.text.y = element_text(face="bold", size=18),
           axis.title.x = element_text(face="bold", size=18),
           axis.title.y = element_text(face="bold", size=18, angle=90),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.justification=c(1,0), 
           legend.position=c(0.95,0.15),
           legend.text = element_text(size = 14),
           legend.title=element_blank(),
           legend.key = element_blank()
           )+
      labs(title=title)
  return(p)
}


setwd("C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NatureCommunications_submission/NatCommu_revision/AppealSubmission_V2/Revision3/Codes")

ImSig.genes <- readRDS('ImSig.rds')


# the plot and analysis for NatMed 103 samples can be seen in the following file:
# C:/Donghai_Desktop2/TREM2validation/CancerCell_91samples/CC_73samples_prediction_results/Yr2020_NewRevision/NewValidation_NatMed/NatMed_Expression_SigValidation.R

###################################################################################################################################################################
###################################################################################################################################################################


### for GSE78220_TREM2.related.genes_ImSig.genes

GSE78220_AltAnalyze <- readRDS('GSE78220_expressionMatrix.rds')

GSE78220_PhenoInfo2 <- readRDS('GSE78220_PhenoInfo2.rds')

GSE78220_AltAnalyze_col_rearranged <- GSE78220_AltAnalyze[,c("Symbol",GSE78220_PhenoInfo2$sample)]

GSE78220_TREM2.related.genes_ImSig.genes <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% ImSig.genes,]


rownames(GSE78220_TREM2.related.genes_ImSig.genes) <- GSE78220_TREM2.related.genes_ImSig.genes$Symbol

GSE78220_TREM2.related.genes_ImSig.genes <- GSE78220_TREM2.related.genes_ImSig.genes[,-1]

GSE78220_TREM2.related.genes_ImSig.genes <- as.matrix(GSE78220_TREM2.related.genes_ImSig.genes)

GSE78220_TREM2.related.genes_ImSig.genes <- GSE78220_TREM2.related.genes_ImSig.genes[rowSums(GSE78220_TREM2.related.genes_ImSig.genes) > 1,]

dim(GSE78220_TREM2.related.genes_ImSig.genes)

head(GSE78220_TREM2.related.genes_ImSig.genes)

tail(GSE78220_TREM2.related.genes_ImSig.genes)

sort(rowSums(GSE78220_TREM2.related.genes_ImSig.genes))



GSE78220_pData = data.frame(class=GSE78220_PhenoInfo2$class, sample=GSE78220_PhenoInfo2$sample,
row.names=GSE78220_PhenoInfo2$sample)

GSE78220_pData

GSE78220_phenoData <- new("AnnotatedDataFrame",data=GSE78220_pData)



all(GSE78220_pData$sample[1:15]==colnames(GSE78220_TREM2.related.genes_ImSig.genes)[1:15])

all(GSE78220_pData$sample[16:28]==colnames(GSE78220_TREM2.related.genes_ImSig.genes)[16:28])


GSE78220_ExpSet_V5 <- ExpressionSet(assayData=as.matrix(GSE78220_TREM2.related.genes_ImSig.genes),phenoData=GSE78220_phenoData)

GSE78220_ExpSet_V5

pData(GSE78220_ExpSet_V5)

dim(pData(GSE78220_ExpSet_V5))

dim(exprs(GSE78220_ExpSet_V5))

nrow(GSE78220_TREM2.related.genes_ImSig.genes)

exprs(GSE78220_ExpSet_V5)


nvalidation <- nvalidate(GSE78220_ExpSet_V5, ngenes = 2:(nrow(GSE78220_TREM2.related.genes_ImSig.genes)+1), ntrain = "balanced", method = "welch.test", dist = "cor")


pData(GSE78220_ExpSet_V5)[["class"]]


predictor_GSE78220_V5 <- fit(GSE78220_ExpSet_V5, method = "welch.test") # must remove 0 or low expression genes < 10, otherwise, fit function does not work!!!

prediction_GSE78220_V5 <- predict(predictor_GSE78220_V5, GSE78220_ExpSet_V5,"PD", ngenes=nrow(GSE78220_TREM2.related.genes_ImSig.genes), dist = "cor")

str(prediction_GSE78220_V5)

prediction_GSE78220_V5@predictor



out_GSE78220_V5 <- as.factor(rep(c(1,2),c(15,13)))

out_GSE78220_V5 

z_GSE78220_V5 <- as.numeric(prediction_GSE78220_V5@prediction[,'z'])

z_GSE78220_V5


Test_GSE78220_V5 <- cbind(out_GSE78220_V5,z_GSE78220_V5)

colnames(Test_GSE78220_V5) <- c('grp','res')

Test_GSE78220_V5 <- as.data.frame(Test_GSE78220_V5)

dim(Test_GSE78220_V5)

head(Test_GSE78220_V5)

tail(Test_GSE78220_V5)

Test_GSE78220_V5


tiff("./output_test/p3_GSE78220_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_GSE78220_V5 <- rocplot.single.V2(Test_GSE78220_V5$grp, Test_GSE78220_V5$res, title = "") 
print(p3_GSE78220_V5)

dev.off()


p3_GSE78220_V5_data <- rocdata(Test_GSE78220_V5$grp, Test_GSE78220_V5$res)

str(p3_GSE78220_V5_data)

summary(prediction_GSE78220_V5)$tab

summary(prediction_GSE78220_V5)$sensitivity

summary(prediction_GSE78220_V5)$specificity

summary(prediction_GSE78220_V5)$AUC


###################################################################################################################################################################

# GSE91061 = BMS038

BMS038.Pre.CountTable.normalized.log <- readRDS('BMS038.Pre.CountTable.normalized.log.rds')

BMS038_phenoData <- readRDS('BMS038_phenoData.rds')


BMS038_TREM2.related.genes_ImSig.genes <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% ImSig.genes,]

BMS038_TREM2.related.genes_ImSig.genes <- BMS038_TREM2.related.genes_ImSig.genes[rowSums(BMS038_TREM2.related.genes_ImSig.genes) > 10,]


dim(BMS038_TREM2.related.genes_ImSig.genes)

head(BMS038_TREM2.related.genes_ImSig.genes)

tail(BMS038_TREM2.related.genes_ImSig.genes)

sort(rowSums(BMS038_TREM2.related.genes_ImSig.genes))


BMS038_ExpSet_V5 <- ExpressionSet(assayData=as.matrix(BMS038_TREM2.related.genes_ImSig.genes),phenoData=BMS038_phenoData)

BMS038_ExpSet_V5

pData(BMS038_ExpSet_V5)

dim(pData(BMS038_ExpSet_V5))

dim(exprs(BMS038_ExpSet_V5))

nrow(BMS038_TREM2.related.genes_ImSig.genes)



nvalidation <- nvalidate(BMS038_ExpSet_V5, ngenes = 2:nrow(BMS038_TREM2.related.genes_ImSig.genes), ntrain = "balanced", method = "welch.test", dist = "cor")


pData(BMS038_ExpSet_V5)[["class"]]


predictor_BMS038_V5 <- fit(BMS038_ExpSet_V5, method = "welch.test") # must remove 0 or low expression genes < 10, otherwise, fit function does not work!!!

prediction_BMS038_V5 <- predict(predictor_BMS038_V5, BMS038_ExpSet_V5,"NonResponder", ngenes=nrow(BMS038_TREM2.related.genes_ImSig.genes), dist = "cor")

str(prediction_BMS038_V5)

prediction_BMS038_V5@predictor



out_BMS038_V5 <- as.factor(rep(c(1,2),c(26,25)))

out_BMS038_V5 

z_BMS038_V5 <- as.numeric(prediction_BMS038_V5@prediction[,'z'])

z_BMS038_V5


Test_BMS038_V5 <- cbind(out_BMS038_V5,z_BMS038_V5)

colnames(Test_BMS038_V5) <- c('grp','res')

Test_BMS038_V5 <- as.data.frame(Test_BMS038_V5)

dim(Test_BMS038_V5)

head(Test_BMS038_V5)

tail(Test_BMS038_V5)

Test_BMS038_V5


tiff("./output_test/BMS038_GSE91061_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_BMS038_V5 <- rocplot.single.V2(Test_BMS038_V5$grp, Test_BMS038_V5$res, title = "")
print(p3_BMS038_V5)

dev.off()


p3_BMS038_V5_data <- rocdata(Test_BMS038_V5$grp, Test_BMS038_V5$res)

str(p3_BMS038_V5_data)

summary(prediction_BMS038_V5)$tab

summary(prediction_BMS038_V5)$sensitivity

summary(prediction_BMS038_V5)$specificity

summary(prediction_BMS038_V5)$AUC



###################################################################################################################################################################

# CC_73samples from PRJEB23709

CC_73samples_GE <- read.table('DATASET-PRJEB23709_Pre_73samples.txt',sep="\t",header=T)

CC_73samples_pData <- readRDS('PRJEB23709_Pre_73samples_phenoData.rds')


CC_73samples_GE <- CC_73samples_GE[CC_73samples_GE$Symbol != '',]

CC_73samples_GE <- CC_73samples_GE[,c(1:47,49:75,48,76:83)]

dim(CC_73samples_GE)

head(CC_73samples_GE,n=20)

tail(CC_73samples_GE)

any(duplicated(CC_73samples_GE$Symbol))

colnames(CC_73samples_GE)[2:74]


CC_73samples_GE_matrix <- CC_73samples_GE[,c('Symbol',CC_73samples_pData$sample)]

dim(CC_73samples_GE_matrix)

head(CC_73samples_GE_matrix)

tail(CC_73samples_GE_matrix)

all(CC_73samples_pData$sample==colnames(CC_73samples_GE_matrix)[-1])



CC_73samples_GE_matrix_ImSig.genes <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% ImSig.genes,]

rownames(CC_73samples_GE_matrix_ImSig.genes) <- CC_73samples_GE_matrix_ImSig.genes$Symbol

CC_73samples_GE_matrix_ImSig.genes <- CC_73samples_GE_matrix_ImSig.genes[,-1]

dim(CC_73samples_GE_matrix_ImSig.genes)

head(CC_73samples_GE_matrix_ImSig.genes)

tail(CC_73samples_GE_matrix_ImSig.genes)

all(CC_73samples_pData$sample==colnames(CC_73samples_GE_matrix_ImSig.genes)) # 108 signature genes to 101

sort(rowSums(CC_73samples_GE_matrix_ImSig.genes))


CC_73samples_phenoData <- new("AnnotatedDataFrame",data=CC_73samples_pData)


CC_73samples_ExpSet_V5 <- ExpressionSet(assayData=as.matrix(CC_73samples_GE_matrix_ImSig.genes),phenoData=CC_73samples_phenoData)


pData(CC_73samples_ExpSet_V5)


nvalidation <- nvalidate(CC_73samples_ExpSet_V5, ngenes = 2:(nrow(CC_73samples_GE_matrix_ImSig.genes)+1), ntrain = "prevalence", method = "welch.test", dist = "cor")
plot(nvalidation, type = "xy")


pData(CC_73samples_ExpSet_V5)[["class"]]



predictor_CC_73samples_V5 <- fit(CC_73samples_ExpSet_V5, method = "welch.test")

prediction_CC_73samples_V5 <- predict(predictor_CC_73samples_V5, CC_73samples_ExpSet_V5, "Nonresponder", ngenes=nrow(CC_73samples_GE_matrix_ImSig.genes), dist = "cor")

str(prediction_CC_73samples_V5)

prediction_CC_73samples_V5@predictor



out_CC_73samples_V5 <- as.factor(rep(c(1,2),c(46,27)))

out_CC_73samples_V5 

z_CC_73samples_V5 <- as.numeric(prediction_CC_73samples_V5@prediction[,'z'])

z_CC_73samples_V5

table(out_CC_73samples_V5)


Test_CC_73samples_V5 <- cbind(out_CC_73samples_V5,z_CC_73samples_V5)

colnames(Test_CC_73samples_V5) <- c('grp','res')

Test_CC_73samples_V5 <- as.data.frame(Test_CC_73samples_V5)

dim(Test_CC_73samples_V5)

head(Test_CC_73samples_V5)

tail(Test_CC_73samples_V5)

Test_CC_73samples_V5


tiff("./output_test/p3_CC_73samples_PRJEB23709_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_CC_73samples_V5_copy <- rocplot.single.V2(Test_CC_73samples_V5$grp, Test_CC_73samples_V5$res, title = "")
print(p3_CC_73samples_V5_copy)

dev.off()


p3_CC_73samples_V5_data <- rocdata(Test_CC_73samples_V5$grp, Test_CC_73samples_V5$res)

str(p3_CC_73samples_V5_data)

summary(prediction_CC_73samples_V5)$tab

summary(prediction_CC_73samples_V5)$sensitivity

summary(prediction_CC_73samples_V5)$specificity

summary(prediction_CC_73samples_V5)$AUC



###################################################################################################################################################################

#### MGSP project: 103 samples


NatMed_103samples_GE_matrix <- readRDS('NatMed_103samples_GE_matrix.rds')

NatMed_103samples_pData <- readRDS('NatMed_103samples_pData.rds')


NatMed_103samples_GE_ImSig.genes <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% ImSig.genes,]

rownames(NatMed_103samples_GE_ImSig.genes) <- NatMed_103samples_GE_ImSig.genes$Symbol

NatMed_103samples_GE_ImSig.genes <- NatMed_103samples_GE_ImSig.genes[,-1]

dim(NatMed_103samples_GE_ImSig.genes)

head(NatMed_103samples_GE_ImSig.genes)

tail(NatMed_103samples_GE_ImSig.genes)

all(NatMed_103samples_pData$sample==colnames(NatMed_103samples_GE_ImSig.genes))


NatMed_103samples_phenoData <- new("AnnotatedDataFrame",data=NatMed_103samples_pData)


NatMed_103samples_ExpSet <- ExpressionSet(assayData=as.matrix(NatMed_103samples_GE_ImSig.genes),phenoData=NatMed_103samples_phenoData)

NatMed_103samples_ExpSet

pData(NatMed_103samples_ExpSet)

dim(pData(NatMed_103samples_ExpSet))

dim(exprs(NatMed_103samples_ExpSet))

head(exprs(NatMed_103samples_ExpSet))

tail(exprs(NatMed_103samples_ExpSet))


nvalidation <- nvalidate(NatMed_103samples_ExpSet, ngenes = 2:99, ntrain = "prevalence", method = "welch.test", dist = "cor") # only 98 genes overlpping with 108 genes
plot(nvalidation, type = "xy")


pData(NatMed_103samples_ExpSet)[["class"]]


predictor_NatMed_103samples <- fit(NatMed_103samples_ExpSet, method = "welch.test")


prediction_NatMed_103samples <- predict(predictor_NatMed_103samples, NatMed_103samples_ExpSet, "Progressor", ngenes=98, dist = "cor")


str(prediction_NatMed_103samples)

prediction_NatMed_103samples@predictor



out_NatMed_103samples <- as.factor(rep(c(1,2),c(47,56)))

out_NatMed_103samples 

z_NatMed_103samples <- as.numeric(prediction_NatMed_103samples@prediction[,'z'])

z_NatMed_103samples

table(out_NatMed_103samples)


Test_NatMed_103samples <- cbind(out_NatMed_103samples,z_NatMed_103samples)

colnames(Test_NatMed_103samples) <- c('grp','res')

Test_NatMed_103samples <- as.data.frame(Test_NatMed_103samples)

dim(Test_NatMed_103samples)

head(Test_NatMed_103samples)

tail(Test_NatMed_103samples)

Test_NatMed_103samples


tiff("./output_test/NatMed_103samples_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_NatMed_V5 <- rocplot.single.V2(Test_NatMed_103samples$grp, Test_NatMed_103samples$res, title = "")
print(p3_NatMed_V5)

dev.off()


p3_NatMed_103samples_data <- rocdata(Test_NatMed_103samples$grp, Test_NatMed_103samples$res)

str(p3_NatMed_103samples_data)

summary(prediction_NatMed_103samples)$tab

summary(prediction_NatMed_103samples)$sensitivity

summary(prediction_NatMed_103samples)$specificity

summary(prediction_NatMed_103samples)$AUC


#########################################################################################################################################################################################
#########################################################################################################################################################################################
