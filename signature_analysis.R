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


human = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jul2015.archive.ensembl.org")
mouse = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host = "jul2015.archive.ensembl.org")


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

rocplot.single <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
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
           axis.title.x = element_text(face="bold", size=12),
           axis.title.y = element_text(face="bold", size=12, angle=90),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.justification=c(1,0), 
           legend.position=c(1,0),
           legend.title=element_blank(),
           legend.key = element_blank()
           )+
      labs(title=title)
  return(p)
}


rocplot.multiple <- function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
                   )
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
 
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
       geom_line(aes(colour = .id)) +
       geom_abline (intercept = 0, slope = 1) +
       theme_bw() +
       scale_x_continuous("False Positive Rate (1-Specificity)") +
       scale_y_continuous("True Positive Rate (Sensitivity)") +
       scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
       theme(
            plot.title = element_text(face="bold", size=14), 
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12, angle=90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.justification=c(1,0), 
            legend.position=c(1,0),
            legend.title=element_blank(),
            legend.key = element_blank()
            )+
      labs(title=title)
  return(p)
}


rocplot.single.V2 <- function(grp, pred, title = "ROC Plot", p.value = FALSE){
  require(ggplot2)
  plotdata <- rocdata(grp, pred)
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
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
           legend.position="none",
           legend.title=element_blank(),
           legend.key = element_blank()
           )+
      labs(title=title)
  return(p)
}



rocplot.multiple.V2 <- 
function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
                   )
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
 
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
       geom_line(aes(colour = .id)) +
       geom_abline (intercept = 0, slope = 1) +
       theme_bw() +
       scale_x_continuous("False Positive Rate (1-Specificity)") +
       scale_y_continuous("True Positive Rate (Sensitivity)") +
       scale_colour_brewer(palette="Set1", breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
       theme(
            plot.title = element_text(face="bold", size=14), 
            axis.text.x = element_text(face="bold", size=18),
          	axis.text.y = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=18),
            axis.title.y = element_text(face="bold", size=18, angle=90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.justification=c(1,0), 
            legend.position=c(0.99,0.02),
            legend.title=element_blank(),
            legend.text = element_text(face="bold",size = 12),
            legend.key = element_blank()
            )+
      labs(title=title)
  return(p)
}


rocplot.multiple.V3 <- 
function(test.data.list, groupName = "grp", predName = "res", title = "ROC Plot", p.value = TRUE){
  require(plyr)
  require(ggplot2)
  plotdata <- llply(test.data.list, function(x) with(x, rocdata(grp = eval(parse(text = groupName)), pred = eval(parse(text = predName)))))
  plotdata <- list(roc = ldply(plotdata, function(x) x$roc),
                   stats = ldply(plotdata, function(x) x$stats)
                   )
 
  if (p.value == TRUE){
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (P=", signif(p.value, 2), ")", sep=""))
  } else {
    annotation <- with(plotdata$stats, paste("AUC=",signif(auc, 2), " (95%CI ", signif(ci.upper, 2), " - ", signif(ci.lower, 2), ")", sep=""))
  }
  
  mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
 
  p <- ggplot(plotdata$roc, aes(x = x, y = y)) +
       geom_line(aes(colour = .id)) +
       geom_abline (intercept = 0, slope = 1) +
       theme_bw() +
       scale_x_continuous("False Positive Rate (1-Specificity)") +
       scale_y_continuous("True Positive Rate (Sensitivity)") +
       scale_color_manual(values = mycolors, breaks = names(test.data.list), labels = paste(names(test.data.list), ": ", annotation, sep = "")) +
       theme(
            plot.title = element_text(face="bold", size=14), 
            axis.text.x = element_text(face="bold", size=18),
          	axis.text.y = element_text(face="bold", size=18),
            axis.title.x = element_text(face="bold", size=18),
            axis.title.y = element_text(face="bold", size=18, angle=90),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.justification=c(1,0), 
            legend.position=c(0.995,0.01),
            legend.title=element_blank(),
            legend.text = element_text(face="bold",size = 12),
            legend.key = element_blank()
            )+
      labs(title=title)
  return(p)
}


#setwd("C:/Donghai_Desktop2/TREM2validation/NewDataSets_Validation")

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


p3_GSE78220_V5_plot1 <- rocplot.single(Test_GSE78220_V5$grp, Test_GSE78220_V5$res, title = "")
print(p3_GSE78220_V5_plot1)

tiff("./output/p3_GSE78220_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_GSE78220_V5 <- rocplot.single.V2(Test_GSE78220_V5$grp, Test_GSE78220_V5$res, title = "")
print(p3_GSE78220_V5)

dev.off()


saveRDS(p3_GSE78220_V5, file='./output/rocplot_GSE78220_ImSig.rds')


p3_GSE78220_V5_data <- rocdata(Test_GSE78220_V5$grp, Test_GSE78220_V5$res)

str(p3_GSE78220_V5_data)

summary(prediction_GSE78220_V5)$tab

summary(prediction_GSE78220_V5)$sensitivity

summary(prediction_GSE78220_V5)$specificity

summary(prediction_GSE78220_V5)$AUC

prediction_GSE78220_V5@predictor


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

p3_BMS038_V5_plot1 <- rocplot.single(Test_BMS038_V5$grp, Test_BMS038_V5$res, title = "")
print(p3_BMS038_V5_plot1)

tiff("./output/BMS038_GSE91061_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_BMS038_V5 <- rocplot.single.V2(Test_BMS038_V5$grp, Test_BMS038_V5$res, title = "")
print(p3_BMS038_V5)

dev.off()


saveRDS(p3_BMS038_V5, 
file='./output/rocplot_BMS038_ImSig.rds')


p3_BMS038_V5_data <- rocdata(Test_BMS038_V5$grp, Test_BMS038_V5$res)

str(p3_BMS038_V5_data)

summary(prediction_BMS038_V5)$tab

summary(prediction_BMS038_V5)$sensitivity

summary(prediction_BMS038_V5)$specificity

summary(prediction_BMS038_V5)$AUC

prediction_BMS038_V5@predictor

nrow(BMS038_TREM2.related.genes_ImSig.genes)


###################################################################################################################################################################

# CC_73samples from PRJEB23709

CC_73samples_GE <- read.table('DATASET-PRJEB23709_Pre_73samples.txt',sep="\t",header=T)

CC_73samples_GE <- CC_73samples_GE[CC_73samples_GE$Symbol != '',]

CC_73samples_GE <- CC_73samples_GE[,c(1:47,49:75,48,76:83)]

dim(CC_73samples_GE)

head(CC_73samples_GE,n=20)

tail(CC_73samples_GE)

any(duplicated(CC_73samples_GE$Symbol))

colnames(CC_73samples_GE)[2:74]


CC_73samples_pData <- readRDS('PRJEB23709_Pre_73samples_phenoData.rds')


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

saveRDS(Test_CC_73samples_V5, './output/rocplot_CC_73samples_PRJEB23709_ImSig.rds')



p3_CC_73samples_V5 <- rocplot.single(Test_CC_73samples_V5$grp, Test_CC_73samples_V5$res, title = "")
print(p3_CC_73samples_V5)


tiff("./output/p3_CC_73samples_PRJEB23709_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_CC_73samples_V5_copy <- rocplot.single.V2(Test_CC_73samples_V5$grp, Test_CC_73samples_V5$res, title = "")
print(p3_CC_73samples_V5_copy)

dev.off()


p3_CC_73samples_V5_data <- rocdata(Test_CC_73samples_V5$grp, Test_CC_73samples_V5$res)

str(p3_CC_73samples_V5_data)

summary(prediction_CC_73samples_V5)$tab

summary(prediction_CC_73samples_V5)$sensitivity

summary(prediction_CC_73samples_V5)$specificity

summary(prediction_CC_73samples_V5)$AUC

prediction_CC_73samples_V5@predictor

dim(prediction_CC_73samples_V5@predictor)



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


p3_NatMed_103samples <- rocplot.single(Test_NatMed_103samples$grp, Test_NatMed_103samples$res, title = "")
print(p3_NatMed_103samples)


tiff("./output/NatMed_103samples_ImSig.tiff", width = 5000, height = 5000, units = "px", res = 800)
 
p3_NatMed_V5 <- rocplot.single.V2(Test_NatMed_103samples$grp, Test_NatMed_103samples$res, title = "")
print(p3_NatMed_V5)

dev.off()


saveRDS(p3_NatMed_V5, file='./output/rocplot_NatMed_MGSP_ImSig.rds')


p3_NatMed_103samples_data <- rocdata(Test_NatMed_103samples$grp, Test_NatMed_103samples$res)

str(p3_NatMed_103samples_data)

summary(prediction_NatMed_103samples)$tab

summary(prediction_NatMed_103samples)$sensitivity

summary(prediction_NatMed_103samples)$specificity

summary(prediction_NatMed_103samples)$AUC

prediction_NatMed_103samples@predictor


#########################################################################################################################################################################################
#########################################################################################################################################################################################

### multiple ROC plots for different datasets 


# Other_1 IFNG.Sig

IFNG.Sig <- c('IFNG', 'STAT1', 'IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA')

IFNG.Sig

# Other_2 CD8.Sig

CD8.Sig <- c("CD8A", "CD8B", "CD3D", "CD3E", "CD3G")

CD8.Sig

# Other_3 PDL1.Sig

# PDL1 i.e., PDL1

PDL1.Sig <- c('PDL1','PDCD1LG2','PDCD1')

PDL1.Sig

# Other_4 CRMA.Sig 

CRMA.Sig <- c('CT1.2', 'MAGEA2', 'MAGEA2A', 'MAGEA2B', 'MAGEA3', 'MAGEA6', 'MAGEA12')

CRMA.Sig

# Other_5 IMPRES.Sig

IMPRES.Sig <- c("BTLA", "CD200", "CD200R1", "CD27", "CD276", "CD28", "CD40", "CD80", "CD86", "CEACAM1", "CTLA4", "IDO1",
"IL2RB", "LAG3", "PVR", "PVRL2", "TIGIT", "TNFRSF18", "TNFRSF4", "TNFRSF9", "PDL1", "HAVCR2", "PDCD1", "PDCD1LG2", "TNFRSF14", "TNFSF4", "TNFSF9", "C10orf54")

IMPRES.Sig


# Other_6 IRG.Sig

IRG.Sig <- c('LEPR','PRLHR','NR2F2','PRL','NRP1','ANGPTL5','IGF1','TNFRSF10B','TNFRSF10A','PLAU','IFI30') # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'

IRG.Sig


# Other_7 LRRC15.CAF.Sig

LRRC15.CAF.Sig <- c('MMP11','COL11A1','C1QTNF3','CTHRC1','COL12A1','COL10A1','COL5A2','GJB2','THBS2','AEBP1','MFAP2','LRRC15','PLAU','ITGA11') # Alias for 'PRLHR' are:'GR3','GPR10','PrRPR'

LRRC15.CAF.Sig


# Other_8_T.cell.inflamed.Sig

T.cell.inflamed.Sig <- c('CD3D','IDO1','CIITA','CD3E','CCL5','GZMK','CD2','HLA-DRA','CXCL13','IL2RG','NKG7','HLA-E','CXCR6','LAG3','TAGAP','CXCL10','STAT1','GZMB')

T.cell.inflamed.Sig


# Other_9 IPRES.Sig

IPRES.Sig <- c('ANGPT2','AXL','CCL13','CCL2','CCL7','CDH1','FAP','FLT1','1L10','LOXL2','RORS','TAGLN','TWIST2','VEGFA','VEGFC','WNT5A')

IPRES.Sig


# Other_10 Inflammatory.Sig

Inflammatory.Sig <- c('CCL5','CCR5','PDL1','CD3D','CD3E','CD8A','CIITA','CTLA4','CXCL10','CXCL11','CXCL13','CXCL9','GZMA','GZMB','HLA-DRA','HKA.DRB1','HLA-E',
'IDO1','IL2RG','ITGAL','LAG3','NKG7','PDCD1','PRF1','PTPRC','STAT1','TAGAP')

Inflammatory.Sig


# Other_11 EMT.Sig

EMT.Sig <- c('CDH1','CDH3','CLDN4','EPCAM','ST14','MAL2','VIM','SNAI2','ZEB2','FN1','MMP2','AGER')

EMT.Sig


# Other_12 Blood.Sig

Blood.Sig <- c('ADAM17', 'CDK2', 'CDKN2A', 'DPP4', 'ERBB2', 'HLA-DRA', 'ICOS', 'ITGA4', 'LARGE', 'MYC', 'NAB2', 'NRAS', 'RHOC', 'TGFB1', 'TIMP1')

Blood.Sig


### No1 - GSE78220

### ImmuneCells.Sig

ImmuneCells.Sig.GSE78220 <- Test_GSE78220_V4 


### Other signatures for comparison in GSE78220

# Other_1 IFNG.Sig

GSE78220.IFNG.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% IFNG.Sig,]

rownames(GSE78220.IFNG.Sig.Exp) <- GSE78220.IFNG.Sig.Exp$Symbol

GSE78220.IFNG.Sig.Exp <- GSE78220.IFNG.Sig.Exp[,-1]

GSE78220.IFNG.Sig.Exp <- as.matrix(GSE78220.IFNG.Sig.Exp)

GSE78220.IFNG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.IFNG.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.IFNG.Sig.Exp.Set <- fit(GSE78220.IFNG.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.IFNG.Sig.Exp.Set <- predict(predictor_GSE78220.IFNG.Sig.Exp.Set, GSE78220.IFNG.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.IFNG.Sig.Exp), dist = "cor")

z_GSE78220_IFNG.Sig <- as.numeric(prediction_GSE78220.IFNG.Sig.Exp.Set@prediction[,'z'])

IFNG.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_IFNG.Sig)

colnames(IFNG.Sig.GSE78220) <- c('grp','res')

IFNG.Sig.GSE78220 <- as.data.frame(IFNG.Sig.GSE78220)


# Other_2 CD8.Sig

GSE78220.CD8.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% CD8.Sig,]

rownames(GSE78220.CD8.Sig.Exp) <- GSE78220.CD8.Sig.Exp$Symbol

GSE78220.CD8.Sig.Exp <- GSE78220.CD8.Sig.Exp[,-1]

GSE78220.CD8.Sig.Exp <- as.matrix(GSE78220.CD8.Sig.Exp)

GSE78220.CD8.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.CD8.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.CD8.Sig.Exp.Set <- fit(GSE78220.CD8.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.CD8.Sig.Exp.Set <- predict(predictor_GSE78220.CD8.Sig.Exp.Set, GSE78220.CD8.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.CD8.Sig.Exp), dist = "cor")

z_GSE78220_CD8.Sig <- as.numeric(prediction_GSE78220.CD8.Sig.Exp.Set@prediction[,'z'])

CD8.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_CD8.Sig)

colnames(CD8.Sig.GSE78220) <- c('grp','res')

CD8.Sig.GSE78220 <- as.data.frame(CD8.Sig.GSE78220)


# Other_3 PDL1.Sig

GSE78220.PDL1.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% PDL1.Sig,]

rownames(GSE78220.PDL1.Sig.Exp) <- GSE78220.PDL1.Sig.Exp$Symbol

GSE78220.PDL1.Sig.Exp <- GSE78220.PDL1.Sig.Exp[,-1]

GSE78220.PDL1.Sig.Exp <- as.matrix(GSE78220.PDL1.Sig.Exp)

GSE78220.PDL1.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.PDL1.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.PDL1.Sig.Exp.Set <- fit(GSE78220.PDL1.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.PDL1.Sig.Exp.Set <- predict(predictor_GSE78220.PDL1.Sig.Exp.Set, GSE78220.PDL1.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.PDL1.Sig.Exp), dist = "cor")

z_GSE78220_PDL1.Sig <- as.numeric(prediction_GSE78220.PDL1.Sig.Exp.Set@prediction[,'z'])

PDL1.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_PDL1.Sig)

colnames(PDL1.Sig.GSE78220) <- c('grp','res')

PDL1.Sig.GSE78220 <- as.data.frame(PDL1.Sig.GSE78220)


# Other_4 CRMA.Sig 

GSE78220.CRMA.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% CRMA.Sig,]

rownames(GSE78220.CRMA.Sig.Exp) <- GSE78220.CRMA.Sig.Exp$Symbol

GSE78220.CRMA.Sig.Exp <- GSE78220.CRMA.Sig.Exp[,-1]

GSE78220.CRMA.Sig.Exp <- as.matrix(GSE78220.CRMA.Sig.Exp)

GSE78220.CRMA.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.CRMA.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.CRMA.Sig.Exp.Set <- fit(GSE78220.CRMA.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.CRMA.Sig.Exp.Set <- predict(predictor_GSE78220.CRMA.Sig.Exp.Set, GSE78220.CRMA.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.CRMA.Sig.Exp), dist = "cor")

z_GSE78220_CRMA.Sig <- as.numeric(prediction_GSE78220.CRMA.Sig.Exp.Set@prediction[,'z'])

CRMA.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_CRMA.Sig)

colnames(CRMA.Sig.GSE78220) <- c('grp','res')

CRMA.Sig.GSE78220 <- as.data.frame(CRMA.Sig.GSE78220)


# Other_5 IMPRES.Sig 

GSE78220.IMPRES.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% IMPRES.Sig,]

rownames(GSE78220.IMPRES.Sig.Exp) <- GSE78220.IMPRES.Sig.Exp$Symbol

GSE78220.IMPRES.Sig.Exp <- GSE78220.IMPRES.Sig.Exp[,-1]

GSE78220.IMPRES.Sig.Exp <- as.matrix(GSE78220.IMPRES.Sig.Exp)

GSE78220.IMPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.IMPRES.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.IMPRES.Sig.Exp.Set <- fit(GSE78220.IMPRES.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.IMPRES.Sig.Exp.Set <- predict(predictor_GSE78220.IMPRES.Sig.Exp.Set, GSE78220.IMPRES.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.IMPRES.Sig.Exp), dist = "cor")

z_GSE78220_IMPRES.Sig <- as.numeric(prediction_GSE78220.IMPRES.Sig.Exp.Set@prediction[,'z'])

IMPRES.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_IMPRES.Sig)

colnames(IMPRES.Sig.GSE78220) <- c('grp','res')

IMPRES.Sig.GSE78220 <- as.data.frame(IMPRES.Sig.GSE78220)


# Other_6 IRG.Sig

GSE78220.IRG.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% IRG.Sig,]

rownames(GSE78220.IRG.Sig.Exp) <- GSE78220.IRG.Sig.Exp$Symbol

GSE78220.IRG.Sig.Exp <- GSE78220.IRG.Sig.Exp[,-1]

GSE78220.IRG.Sig.Exp <- as.matrix(GSE78220.IRG.Sig.Exp)

GSE78220.IRG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.IRG.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.IRG.Sig.Exp.Set <- fit(GSE78220.IRG.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.IRG.Sig.Exp.Set <- predict(predictor_GSE78220.IRG.Sig.Exp.Set, GSE78220.IRG.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.IRG.Sig.Exp), dist = "cor")

z_GSE78220_IRG.Sig <- as.numeric(prediction_GSE78220.IRG.Sig.Exp.Set@prediction[,'z'])

IRG.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_IRG.Sig)

colnames(IRG.Sig.GSE78220) <- c('grp','res')

IRG.Sig.GSE78220 <- as.data.frame(IRG.Sig.GSE78220)


# Other_7 LRRC15.CAF.Sig

GSE78220.LRRC15.CAF.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% LRRC15.CAF.Sig,]

rownames(GSE78220.LRRC15.CAF.Sig.Exp) <- GSE78220.LRRC15.CAF.Sig.Exp$Symbol

GSE78220.LRRC15.CAF.Sig.Exp <- GSE78220.LRRC15.CAF.Sig.Exp[,-1]

GSE78220.LRRC15.CAF.Sig.Exp <- as.matrix(GSE78220.LRRC15.CAF.Sig.Exp)

GSE78220.LRRC15.CAF.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.LRRC15.CAF.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.LRRC15.CAF.Sig.Exp.Set <- fit(GSE78220.LRRC15.CAF.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.LRRC15.CAF.Sig.Exp.Set <- predict(predictor_GSE78220.LRRC15.CAF.Sig.Exp.Set, GSE78220.LRRC15.CAF.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.LRRC15.CAF.Sig.Exp), dist = "cor")

z_GSE78220_LRRC15.CAF.Sig <- as.numeric(prediction_GSE78220.LRRC15.CAF.Sig.Exp.Set@prediction[,'z'])

LRRC15.CAF.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_LRRC15.CAF.Sig)

colnames(LRRC15.CAF.Sig.GSE78220) <- c('grp','res')

LRRC15.CAF.Sig.GSE78220 <- as.data.frame(LRRC15.CAF.Sig.GSE78220)


# Other_8_T.cell.inflamed.Sig

GSE78220.T.cell.inflamed.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% T.cell.inflamed.Sig,]

rownames(GSE78220.T.cell.inflamed.Sig.Exp) <- GSE78220.T.cell.inflamed.Sig.Exp$Symbol

GSE78220.T.cell.inflamed.Sig.Exp <- GSE78220.T.cell.inflamed.Sig.Exp[,-1]

GSE78220.T.cell.inflamed.Sig.Exp <- as.matrix(GSE78220.T.cell.inflamed.Sig.Exp)

GSE78220.T.cell.inflamed.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.T.cell.inflamed.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.T.cell.inflamed.Sig.Exp.Set <- fit(GSE78220.T.cell.inflamed.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.T.cell.inflamed.Sig.Exp.Set <- predict(predictor_GSE78220.T.cell.inflamed.Sig.Exp.Set, GSE78220.T.cell.inflamed.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.T.cell.inflamed.Sig.Exp), dist = "cor")

z_GSE78220_T.cell.inflamed.Sig <- as.numeric(prediction_GSE78220.T.cell.inflamed.Sig.Exp.Set@prediction[,'z'])

T.cell.inflamed.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_T.cell.inflamed.Sig)

colnames(T.cell.inflamed.Sig.GSE78220) <- c('grp','res')

T.cell.inflamed.Sig.GSE78220 <- as.data.frame(T.cell.inflamed.Sig.GSE78220)


# Other_9 IPRES.Sig

GSE78220.IPRES.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% IPRES.Sig,]

rownames(GSE78220.IPRES.Sig.Exp) <- GSE78220.IPRES.Sig.Exp$Symbol

GSE78220.IPRES.Sig.Exp <- GSE78220.IPRES.Sig.Exp[,-1]

GSE78220.IPRES.Sig.Exp <- as.matrix(GSE78220.IPRES.Sig.Exp)

GSE78220.IPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.IPRES.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.IPRES.Sig.Exp.Set <- fit(GSE78220.IPRES.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.IPRES.Sig.Exp.Set <- predict(predictor_GSE78220.IPRES.Sig.Exp.Set, GSE78220.IPRES.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.IPRES.Sig.Exp), dist = "cor")

z_GSE78220_IPRES.Sig <- as.numeric(prediction_GSE78220.IPRES.Sig.Exp.Set@prediction[,'z'])

IPRES.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_IPRES.Sig)

colnames(IPRES.Sig.GSE78220) <- c('grp','res')

IPRES.Sig.GSE78220 <- as.data.frame(IPRES.Sig.GSE78220)


# Other_10 Inflammatory.Sig

GSE78220.Inflammatory.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% Inflammatory.Sig,]

rownames(GSE78220.Inflammatory.Sig.Exp) <- GSE78220.Inflammatory.Sig.Exp$Symbol

GSE78220.Inflammatory.Sig.Exp <- GSE78220.Inflammatory.Sig.Exp[,-1]

GSE78220.Inflammatory.Sig.Exp <- as.matrix(GSE78220.Inflammatory.Sig.Exp)

GSE78220.Inflammatory.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.Inflammatory.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.Inflammatory.Sig.Exp.Set <- fit(GSE78220.Inflammatory.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.Inflammatory.Sig.Exp.Set <- predict(predictor_GSE78220.Inflammatory.Sig.Exp.Set, GSE78220.Inflammatory.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.Inflammatory.Sig.Exp), dist = "cor")

z_GSE78220_Inflammatory.Sig <- as.numeric(prediction_GSE78220.Inflammatory.Sig.Exp.Set@prediction[,'z'])

Inflammatory.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_Inflammatory.Sig)

colnames(Inflammatory.Sig.GSE78220) <- c('grp','res')

Inflammatory.Sig.GSE78220 <- as.data.frame(Inflammatory.Sig.GSE78220)


# Other_11 EMT.Sig

GSE78220.EMT.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% EMT.Sig,]

rownames(GSE78220.EMT.Sig.Exp) <- GSE78220.EMT.Sig.Exp$Symbol

GSE78220.EMT.Sig.Exp <- GSE78220.EMT.Sig.Exp[,-1]

GSE78220.EMT.Sig.Exp <- as.matrix(GSE78220.EMT.Sig.Exp)

GSE78220.EMT.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.EMT.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.EMT.Sig.Exp.Set <- fit(GSE78220.EMT.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.EMT.Sig.Exp.Set <- predict(predictor_GSE78220.EMT.Sig.Exp.Set, GSE78220.EMT.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.EMT.Sig.Exp), dist = "cor")

z_GSE78220_EMT.Sig <- as.numeric(prediction_GSE78220.EMT.Sig.Exp.Set@prediction[,'z'])

EMT.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_EMT.Sig)

colnames(EMT.Sig.GSE78220) <- c('grp','res')

EMT.Sig.GSE78220 <- as.data.frame(EMT.Sig.GSE78220)


# Other_12 Blood.Sig

GSE78220.Blood.Sig.Exp <- GSE78220_AltAnalyze_col_rearranged[GSE78220_AltAnalyze_col_rearranged$Symbol %in% Blood.Sig,]

rownames(GSE78220.Blood.Sig.Exp) <- GSE78220.Blood.Sig.Exp$Symbol

GSE78220.Blood.Sig.Exp <- GSE78220.Blood.Sig.Exp[,-1]

GSE78220.Blood.Sig.Exp <- as.matrix(GSE78220.Blood.Sig.Exp)

GSE78220.Blood.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE78220.Blood.Sig.Exp),phenoData=GSE78220_phenoData)

predictor_GSE78220.Blood.Sig.Exp.Set <- fit(GSE78220.Blood.Sig.Exp.Set, method = "welch.test") 

prediction_GSE78220.Blood.Sig.Exp.Set <- predict(predictor_GSE78220.Blood.Sig.Exp.Set, GSE78220.Blood.Sig.Exp.Set,"PD", ngenes=nrow(GSE78220.Blood.Sig.Exp), dist = "cor")

z_GSE78220_Blood.Sig <- as.numeric(prediction_GSE78220.Blood.Sig.Exp.Set@prediction[,'z'])

Blood.Sig.GSE78220 <- cbind(out_GSE78220_V5,z_GSE78220_Blood.Sig)

colnames(Blood.Sig.GSE78220) <- c('grp','res')

Blood.Sig.GSE78220 <- as.data.frame(Blood.Sig.GSE78220)


# Multiplot - GSE78220

GSE78220_Data <- list(ImmuneCells.Sig=ImmuneCells.Sig.GSE78220, IMPRES.Sig=IMPRES.Sig.GSE78220, IRG.Sig=IRG.Sig.GSE78220, IPRES.Sig=IPRES.Sig.GSE78220, 
T.cell.inflamed.Sig=T.cell.inflamed.Sig.GSE78220, LRRC15.CAF.Sig=LRRC15.CAF.Sig.GSE78220, EMT.Sig=EMT.Sig.GSE78220, Inflammatory.Sig=Inflammatory.Sig.GSE78220, 
Blood.Sig=Blood.Sig.GSE78220, IFNG.Sig=IFNG.Sig.GSE78220, CD8.Sig=CD8.Sig.GSE78220, PDL1.Sig=PDL1.Sig.GSE78220, CRMA.Sig=CRMA.Sig.GSE78220)
                  
       
tiff("./output/No1_GSE78220_28samples.tiff", width = 7600, height = 7600, units = "px", res = 710)
 
p <- rocplot.multiple.V3(GSE78220_Data, title = "", p.value = FALSE)
print(p)

dev.off()



### No2 - GSE91061 = BMS038

### ImmuneCells.Sig

ImmuneCells.Sig.GSE91061 <- Test_BMS038_V5 


### Other signatures for comparison in GSE91061

# Other_1 IFNG.Sig

GSE91061.IFNG.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% IFNG.Sig,]

GSE91061.IFNG.Sig.Exp <- as.matrix(GSE91061.IFNG.Sig.Exp)

GSE91061.IFNG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.IFNG.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.IFNG.Sig.Exp.Set <- fit(GSE91061.IFNG.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.IFNG.Sig.Exp.Set <- predict(predictor_GSE91061.IFNG.Sig.Exp.Set, GSE91061.IFNG.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.IFNG.Sig.Exp), dist = "cor")

z_GSE91061_IFNG.Sig <- as.numeric(prediction_GSE91061.IFNG.Sig.Exp.Set@prediction[,'z'])

IFNG.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_IFNG.Sig)

colnames(IFNG.Sig.GSE91061) <- c('grp','res')

IFNG.Sig.GSE91061 <- as.data.frame(IFNG.Sig.GSE91061)


# Other_2 CD8.Sig

GSE91061.CD8.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% CD8.Sig,]

GSE91061.CD8.Sig.Exp <- as.matrix(GSE91061.CD8.Sig.Exp)

GSE91061.CD8.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.CD8.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.CD8.Sig.Exp.Set <- fit(GSE91061.CD8.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.CD8.Sig.Exp.Set <- predict(predictor_GSE91061.CD8.Sig.Exp.Set, GSE91061.CD8.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.CD8.Sig.Exp), dist = "cor")

z_GSE91061_CD8.Sig <- as.numeric(prediction_GSE91061.CD8.Sig.Exp.Set@prediction[,'z'])

CD8.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_CD8.Sig)

colnames(CD8.Sig.GSE91061) <- c('grp','res')

CD8.Sig.GSE91061 <- as.data.frame(CD8.Sig.GSE91061)


# Other_3 PDL1.Sig

GSE91061.PDL1.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% PDL1.Sig,]

GSE91061.PDL1.Sig.Exp <- as.matrix(GSE91061.PDL1.Sig.Exp)

GSE91061.PDL1.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.PDL1.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.PDL1.Sig.Exp.Set <- fit(GSE91061.PDL1.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.PDL1.Sig.Exp.Set <- predict(predictor_GSE91061.PDL1.Sig.Exp.Set, GSE91061.PDL1.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.PDL1.Sig.Exp), dist = "cor")

z_GSE91061_PDL1.Sig <- as.numeric(prediction_GSE91061.PDL1.Sig.Exp.Set@prediction[,'z'])

PDL1.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_PDL1.Sig)

colnames(PDL1.Sig.GSE91061) <- c('grp','res')

PDL1.Sig.GSE91061 <- as.data.frame(PDL1.Sig.GSE91061)


# Other_4 CRMA.Sig 

GSE91061.CRMA.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% CRMA.Sig,]

GSE91061.CRMA.Sig.Exp <- as.matrix(GSE91061.CRMA.Sig.Exp)

GSE91061.CRMA.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.CRMA.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.CRMA.Sig.Exp.Set <- fit(GSE91061.CRMA.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.CRMA.Sig.Exp.Set <- predict(predictor_GSE91061.CRMA.Sig.Exp.Set, GSE91061.CRMA.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.CRMA.Sig.Exp), dist = "cor")

z_GSE91061_CRMA.Sig <- as.numeric(prediction_GSE91061.CRMA.Sig.Exp.Set@prediction[,'z'])

CRMA.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_CRMA.Sig)

colnames(CRMA.Sig.GSE91061) <- c('grp','res')

CRMA.Sig.GSE91061 <- as.data.frame(CRMA.Sig.GSE91061)

CRMA.Sig.GSE91061 <- CRMA.Sig.GSE91061[!is.na(CRMA.Sig.GSE91061$res),]


# Other_5 IMPRES.Sig 

GSE91061.IMPRES.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% IMPRES.Sig,]

GSE91061.IMPRES.Sig.Exp <- as.matrix(GSE91061.IMPRES.Sig.Exp)

GSE91061.IMPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.IMPRES.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.IMPRES.Sig.Exp.Set <- fit(GSE91061.IMPRES.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.IMPRES.Sig.Exp.Set <- predict(predictor_GSE91061.IMPRES.Sig.Exp.Set, GSE91061.IMPRES.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.IMPRES.Sig.Exp), dist = "cor")

z_GSE91061_IMPRES.Sig <- as.numeric(prediction_GSE91061.IMPRES.Sig.Exp.Set@prediction[,'z'])

IMPRES.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_IMPRES.Sig)

colnames(IMPRES.Sig.GSE91061) <- c('grp','res')

IMPRES.Sig.GSE91061 <- as.data.frame(IMPRES.Sig.GSE91061)


# Other_6 IRG.Sig

GSE91061.IRG.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% IRG.Sig,]

GSE91061.IRG.Sig.Exp <- as.matrix(GSE91061.IRG.Sig.Exp)

GSE91061.IRG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.IRG.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.IRG.Sig.Exp.Set <- fit(GSE91061.IRG.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.IRG.Sig.Exp.Set <- predict(predictor_GSE91061.IRG.Sig.Exp.Set, GSE91061.IRG.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.IRG.Sig.Exp), dist = "cor")

z_GSE91061_IRG.Sig <- as.numeric(prediction_GSE91061.IRG.Sig.Exp.Set@prediction[,'z'])

IRG.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_IRG.Sig)

colnames(IRG.Sig.GSE91061) <- c('grp','res')

IRG.Sig.GSE91061 <- as.data.frame(IRG.Sig.GSE91061)


# Other_7 LRRC15.CAF.Sig

GSE91061.LRRC15.CAF.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% LRRC15.CAF.Sig,]

GSE91061.LRRC15.CAF.Sig.Exp <- as.matrix(GSE91061.LRRC15.CAF.Sig.Exp)

GSE91061.LRRC15.CAF.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.LRRC15.CAF.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.LRRC15.CAF.Sig.Exp.Set <- fit(GSE91061.LRRC15.CAF.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.LRRC15.CAF.Sig.Exp.Set <- predict(predictor_GSE91061.LRRC15.CAF.Sig.Exp.Set, GSE91061.LRRC15.CAF.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.LRRC15.CAF.Sig.Exp), dist = "cor")

z_GSE91061_LRRC15.CAF.Sig <- as.numeric(prediction_GSE91061.LRRC15.CAF.Sig.Exp.Set@prediction[,'z'])

LRRC15.CAF.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_LRRC15.CAF.Sig)

colnames(LRRC15.CAF.Sig.GSE91061) <- c('grp','res')

LRRC15.CAF.Sig.GSE91061 <- as.data.frame(LRRC15.CAF.Sig.GSE91061)


# Other_8_T.cell.inflamed.Sig

GSE91061.T.cell.inflamed.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% T.cell.inflamed.Sig,]

GSE91061.T.cell.inflamed.Sig.Exp <- as.matrix(GSE91061.T.cell.inflamed.Sig.Exp)

GSE91061.T.cell.inflamed.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.T.cell.inflamed.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.T.cell.inflamed.Sig.Exp.Set <- fit(GSE91061.T.cell.inflamed.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.T.cell.inflamed.Sig.Exp.Set <- predict(predictor_GSE91061.T.cell.inflamed.Sig.Exp.Set, GSE91061.T.cell.inflamed.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.T.cell.inflamed.Sig.Exp), dist = "cor")

z_GSE91061_T.cell.inflamed.Sig <- as.numeric(prediction_GSE91061.T.cell.inflamed.Sig.Exp.Set@prediction[,'z'])

T.cell.inflamed.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_T.cell.inflamed.Sig)

colnames(T.cell.inflamed.Sig.GSE91061) <- c('grp','res')

T.cell.inflamed.Sig.GSE91061 <- as.data.frame(T.cell.inflamed.Sig.GSE91061)


# Other_9 IPRES.Sig

GSE91061.IPRES.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% IPRES.Sig,]

GSE91061.IPRES.Sig.Exp <- as.matrix(GSE91061.IPRES.Sig.Exp)

GSE91061.IPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.IPRES.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.IPRES.Sig.Exp.Set <- fit(GSE91061.IPRES.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.IPRES.Sig.Exp.Set <- predict(predictor_GSE91061.IPRES.Sig.Exp.Set, GSE91061.IPRES.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.IPRES.Sig.Exp), dist = "cor")

z_GSE91061_IPRES.Sig <- as.numeric(prediction_GSE91061.IPRES.Sig.Exp.Set@prediction[,'z'])

IPRES.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_IPRES.Sig)

colnames(IPRES.Sig.GSE91061) <- c('grp','res')

IPRES.Sig.GSE91061 <- as.data.frame(IPRES.Sig.GSE91061)


# Other_10 Inflammatory.Sig

GSE91061.Inflammatory.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% Inflammatory.Sig,]

GSE91061.Inflammatory.Sig.Exp <- as.matrix(GSE91061.Inflammatory.Sig.Exp)

GSE91061.Inflammatory.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.Inflammatory.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.Inflammatory.Sig.Exp.Set <- fit(GSE91061.Inflammatory.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.Inflammatory.Sig.Exp.Set <- predict(predictor_GSE91061.Inflammatory.Sig.Exp.Set, GSE91061.Inflammatory.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.Inflammatory.Sig.Exp), dist = "cor")

z_GSE91061_Inflammatory.Sig <- as.numeric(prediction_GSE91061.Inflammatory.Sig.Exp.Set@prediction[,'z'])

Inflammatory.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_Inflammatory.Sig)

colnames(Inflammatory.Sig.GSE91061) <- c('grp','res')

Inflammatory.Sig.GSE91061 <- as.data.frame(Inflammatory.Sig.GSE91061)


# Other_11 EMT.Sig

GSE91061.EMT.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% EMT.Sig,]

GSE91061.EMT.Sig.Exp <- as.matrix(GSE91061.EMT.Sig.Exp)

GSE91061.EMT.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.EMT.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.EMT.Sig.Exp.Set <- fit(GSE91061.EMT.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.EMT.Sig.Exp.Set <- predict(predictor_GSE91061.EMT.Sig.Exp.Set, GSE91061.EMT.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.EMT.Sig.Exp), dist = "cor")

z_GSE91061_EMT.Sig <- as.numeric(prediction_GSE91061.EMT.Sig.Exp.Set@prediction[,'z'])

EMT.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_EMT.Sig)

colnames(EMT.Sig.GSE91061) <- c('grp','res')

EMT.Sig.GSE91061 <- as.data.frame(EMT.Sig.GSE91061)


# Other_12 Blood.Sig

GSE91061.Blood.Sig.Exp <- BMS038.Pre.CountTable.normalized.log[rownames(BMS038.Pre.CountTable.normalized.log) %in% Blood.Sig,]

GSE91061.Blood.Sig.Exp <- as.matrix(GSE91061.Blood.Sig.Exp)

GSE91061.Blood.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(GSE91061.Blood.Sig.Exp),phenoData=BMS038_phenoData)

predictor_GSE91061.Blood.Sig.Exp.Set <- fit(GSE91061.Blood.Sig.Exp.Set, method = "welch.test") 

prediction_GSE91061.Blood.Sig.Exp.Set <- predict(predictor_GSE91061.Blood.Sig.Exp.Set, GSE91061.Blood.Sig.Exp.Set,"NonResponder", ngenes=nrow(GSE91061.Blood.Sig.Exp), dist = "cor")

z_GSE91061_Blood.Sig <- as.numeric(prediction_GSE91061.Blood.Sig.Exp.Set@prediction[,'z'])

Blood.Sig.GSE91061 <- cbind(out_BMS038_V5,z_GSE91061_Blood.Sig)

colnames(Blood.Sig.GSE91061) <- c('grp','res')

Blood.Sig.GSE91061 <- as.data.frame(Blood.Sig.GSE91061)


# Multiplot - GSE91061

GSE91061_Data <- list(ImmuneCells.Sig=ImmuneCells.Sig.GSE91061, IMPRES.Sig=IMPRES.Sig.GSE91061, IRG.Sig=IRG.Sig.GSE91061, IPRES.Sig=IPRES.Sig.GSE91061, 
T.cell.inflamed.Sig=T.cell.inflamed.Sig.GSE91061, LRRC15.CAF.Sig=LRRC15.CAF.Sig.GSE91061, EMT.Sig=EMT.Sig.GSE91061, Inflammatory.Sig=Inflammatory.Sig.GSE91061, 
Blood.Sig=Blood.Sig.GSE91061, IFNG.Sig=IFNG.Sig.GSE91061, CD8.Sig=CD8.Sig.GSE91061, PDL1.Sig=PDL1.Sig.GSE91061, CRMA.Sig=CRMA.Sig.GSE91061)
                  
       
tiff("./output/No2_GSE91061_51samples.tiff", width = 7600, height = 7600, units = "px", res = 710)
 
p <- rocplot.multiple.V3(GSE91061_Data, title = "", p.value = FALSE)
print(p)

dev.off()


### No3 - PRJEB23709

### ImmuneCells.Sig

ImmuneCells.Sig.PRJEB23709 <- Test_CC_73samples_V5 


### Other signatures for comparison in PRJEB23709

# Other_1 IFNG.Sig

PRJEB23709.IFNG.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% IFNG.Sig,]

rownames(PRJEB23709.IFNG.Sig.Exp) <- PRJEB23709.IFNG.Sig.Exp$Symbol

PRJEB23709.IFNG.Sig.Exp <- PRJEB23709.IFNG.Sig.Exp[,-1]

PRJEB23709.IFNG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.IFNG.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.IFNG.Sig.Exp.Set <- fit(PRJEB23709.IFNG.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.IFNG.Sig.Exp.Set <- predict(predictor_PRJEB23709.IFNG.Sig.Exp.Set, PRJEB23709.IFNG.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.IFNG.Sig.Exp), dist = "cor")

z_PRJEB23709_IFNG.Sig <- as.numeric(prediction_PRJEB23709.IFNG.Sig.Exp.Set@prediction[,'z'])

IFNG.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_IFNG.Sig)

colnames(IFNG.Sig.PRJEB23709) <- c('grp','res')

IFNG.Sig.PRJEB23709 <- as.data.frame(IFNG.Sig.PRJEB23709)


# Other_2 CD8.Sig

PRJEB23709.CD8.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% CD8.Sig,]

rownames(PRJEB23709.CD8.Sig.Exp) <- PRJEB23709.CD8.Sig.Exp$Symbol

PRJEB23709.CD8.Sig.Exp <- PRJEB23709.CD8.Sig.Exp[,-1]

PRJEB23709.CD8.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.CD8.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.CD8.Sig.Exp.Set <- fit(PRJEB23709.CD8.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.CD8.Sig.Exp.Set <- predict(predictor_PRJEB23709.CD8.Sig.Exp.Set, PRJEB23709.CD8.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.CD8.Sig.Exp), dist = "cor")

z_PRJEB23709_CD8.Sig <- as.numeric(prediction_PRJEB23709.CD8.Sig.Exp.Set@prediction[,'z'])

CD8.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_CD8.Sig)

colnames(CD8.Sig.PRJEB23709) <- c('grp','res')

CD8.Sig.PRJEB23709 <- as.data.frame(CD8.Sig.PRJEB23709)


# Other_3 PDL1.Sig

PRJEB23709.PDL1.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% PDL1.Sig,]

rownames(PRJEB23709.PDL1.Sig.Exp) <- PRJEB23709.PDL1.Sig.Exp$Symbol

PRJEB23709.PDL1.Sig.Exp <- PRJEB23709.PDL1.Sig.Exp[,-1]

PRJEB23709.PDL1.Sig.Exp <- as.matrix(PRJEB23709.PDL1.Sig.Exp)

PRJEB23709.PDL1.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.PDL1.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.PDL1.Sig.Exp.Set <- fit(PRJEB23709.PDL1.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.PDL1.Sig.Exp.Set <- predict(predictor_PRJEB23709.PDL1.Sig.Exp.Set, PRJEB23709.PDL1.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.PDL1.Sig.Exp), dist = "cor")

z_PRJEB23709_PDL1.Sig <- as.numeric(prediction_PRJEB23709.PDL1.Sig.Exp.Set@prediction[,'z'])

PDL1.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_PDL1.Sig)

colnames(PDL1.Sig.PRJEB23709) <- c('grp','res')

PDL1.Sig.PRJEB23709 <- as.data.frame(PDL1.Sig.PRJEB23709)


# Other_4 CRMA.Sig 

PRJEB23709.CRMA.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% CRMA.Sig,]

rownames(PRJEB23709.CRMA.Sig.Exp) <- PRJEB23709.CRMA.Sig.Exp$Symbol

PRJEB23709.CRMA.Sig.Exp <- PRJEB23709.CRMA.Sig.Exp[,-1]

PRJEB23709.CRMA.Sig.Exp <- as.matrix(PRJEB23709.CRMA.Sig.Exp)

PRJEB23709.CRMA.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.CRMA.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.CRMA.Sig.Exp.Set <- fit(PRJEB23709.CRMA.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.CRMA.Sig.Exp.Set <- predict(predictor_PRJEB23709.CRMA.Sig.Exp.Set, PRJEB23709.CRMA.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.CRMA.Sig.Exp), dist = "cor")

z_PRJEB23709_CRMA.Sig <- as.numeric(prediction_PRJEB23709.CRMA.Sig.Exp.Set@prediction[,'z'])

CRMA.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_CRMA.Sig)

colnames(CRMA.Sig.PRJEB23709) <- c('grp','res')

CRMA.Sig.PRJEB23709 <- as.data.frame(CRMA.Sig.PRJEB23709)

CRMA.Sig.PRJEB23709 <- CRMA.Sig.PRJEB23709[!is.na(CRMA.Sig.PRJEB23709$res),]


# Other_5 IMPRES.Sig 

PRJEB23709.IMPRES.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% IMPRES.Sig,]

rownames(PRJEB23709.IMPRES.Sig.Exp) <- PRJEB23709.IMPRES.Sig.Exp$Symbol

PRJEB23709.IMPRES.Sig.Exp <- PRJEB23709.IMPRES.Sig.Exp[,-1]

PRJEB23709.IMPRES.Sig.Exp <- as.matrix(PRJEB23709.IMPRES.Sig.Exp)

PRJEB23709.IMPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.IMPRES.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.IMPRES.Sig.Exp.Set <- fit(PRJEB23709.IMPRES.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.IMPRES.Sig.Exp.Set <- predict(predictor_PRJEB23709.IMPRES.Sig.Exp.Set, PRJEB23709.IMPRES.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.IMPRES.Sig.Exp), dist = "cor")

z_PRJEB23709_IMPRES.Sig <- as.numeric(prediction_PRJEB23709.IMPRES.Sig.Exp.Set@prediction[,'z'])

IMPRES.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_IMPRES.Sig)

colnames(IMPRES.Sig.PRJEB23709) <- c('grp','res')

IMPRES.Sig.PRJEB23709 <- as.data.frame(IMPRES.Sig.PRJEB23709)


# Other_6 IRG.Sig

PRJEB23709.IRG.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% IRG.Sig,]

rownames(PRJEB23709.IRG.Sig.Exp) <- PRJEB23709.IRG.Sig.Exp$Symbol

PRJEB23709.IRG.Sig.Exp <- PRJEB23709.IRG.Sig.Exp[,-1]

PRJEB23709.IRG.Sig.Exp <- as.matrix(PRJEB23709.IRG.Sig.Exp)

PRJEB23709.IRG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.IRG.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.IRG.Sig.Exp.Set <- fit(PRJEB23709.IRG.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.IRG.Sig.Exp.Set <- predict(predictor_PRJEB23709.IRG.Sig.Exp.Set, PRJEB23709.IRG.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.IRG.Sig.Exp), dist = "cor")

z_PRJEB23709_IRG.Sig <- as.numeric(prediction_PRJEB23709.IRG.Sig.Exp.Set@prediction[,'z'])

IRG.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_IRG.Sig)

colnames(IRG.Sig.PRJEB23709) <- c('grp','res')

IRG.Sig.PRJEB23709 <- as.data.frame(IRG.Sig.PRJEB23709)


# Other_7 LRRC15.CAF.Sig

PRJEB23709.LRRC15.CAF.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% LRRC15.CAF.Sig,]

rownames(PRJEB23709.LRRC15.CAF.Sig.Exp) <- PRJEB23709.LRRC15.CAF.Sig.Exp$Symbol

PRJEB23709.LRRC15.CAF.Sig.Exp <- PRJEB23709.LRRC15.CAF.Sig.Exp[,-1]

PRJEB23709.LRRC15.CAF.Sig.Exp <- as.matrix(PRJEB23709.LRRC15.CAF.Sig.Exp)

PRJEB23709.LRRC15.CAF.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.LRRC15.CAF.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.LRRC15.CAF.Sig.Exp.Set <- fit(PRJEB23709.LRRC15.CAF.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.LRRC15.CAF.Sig.Exp.Set <- predict(predictor_PRJEB23709.LRRC15.CAF.Sig.Exp.Set, PRJEB23709.LRRC15.CAF.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.LRRC15.CAF.Sig.Exp), dist = "cor")

z_PRJEB23709_LRRC15.CAF.Sig <- as.numeric(prediction_PRJEB23709.LRRC15.CAF.Sig.Exp.Set@prediction[,'z'])

LRRC15.CAF.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_LRRC15.CAF.Sig)

colnames(LRRC15.CAF.Sig.PRJEB23709) <- c('grp','res')

LRRC15.CAF.Sig.PRJEB23709 <- as.data.frame(LRRC15.CAF.Sig.PRJEB23709)


# Other_8_T.cell.inflamed.Sig


PRJEB23709.T.cell.inflamed.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% T.cell.inflamed.Sig,]

rownames(PRJEB23709.T.cell.inflamed.Sig.Exp) <- PRJEB23709.T.cell.inflamed.Sig.Exp$Symbol

PRJEB23709.T.cell.inflamed.Sig.Exp <- PRJEB23709.T.cell.inflamed.Sig.Exp[,-1]

PRJEB23709.T.cell.inflamed.Sig.Exp <- as.matrix(PRJEB23709.T.cell.inflamed.Sig.Exp)

PRJEB23709.T.cell.inflamed.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.T.cell.inflamed.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.T.cell.inflamed.Sig.Exp.Set <- fit(PRJEB23709.T.cell.inflamed.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.T.cell.inflamed.Sig.Exp.Set <- predict(predictor_PRJEB23709.T.cell.inflamed.Sig.Exp.Set, PRJEB23709.T.cell.inflamed.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.T.cell.inflamed.Sig.Exp), dist = "cor")

z_PRJEB23709_T.cell.inflamed.Sig <- as.numeric(prediction_PRJEB23709.T.cell.inflamed.Sig.Exp.Set@prediction[,'z'])

T.cell.inflamed.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_T.cell.inflamed.Sig)

colnames(T.cell.inflamed.Sig.PRJEB23709) <- c('grp','res')

T.cell.inflamed.Sig.PRJEB23709 <- as.data.frame(T.cell.inflamed.Sig.PRJEB23709)


# Other_9 IPRES.Sig

PRJEB23709.IPRES.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% IPRES.Sig,]

rownames(PRJEB23709.IPRES.Sig.Exp) <- PRJEB23709.IPRES.Sig.Exp$Symbol

PRJEB23709.IPRES.Sig.Exp <- PRJEB23709.IPRES.Sig.Exp[,-1]

PRJEB23709.IPRES.Sig.Exp <- as.matrix(PRJEB23709.IPRES.Sig.Exp)

PRJEB23709.IPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.IPRES.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.IPRES.Sig.Exp.Set <- fit(PRJEB23709.IPRES.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.IPRES.Sig.Exp.Set <- predict(predictor_PRJEB23709.IPRES.Sig.Exp.Set, PRJEB23709.IPRES.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.IPRES.Sig.Exp), dist = "cor")

z_PRJEB23709_IPRES.Sig <- as.numeric(prediction_PRJEB23709.IPRES.Sig.Exp.Set@prediction[,'z'])

IPRES.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_IPRES.Sig)

colnames(IPRES.Sig.PRJEB23709) <- c('grp','res')

IPRES.Sig.PRJEB23709 <- as.data.frame(IPRES.Sig.PRJEB23709)


# Other_10 Inflammatory.Sig

PRJEB23709.Inflammatory.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% Inflammatory.Sig,]

rownames(PRJEB23709.Inflammatory.Sig.Exp) <- PRJEB23709.Inflammatory.Sig.Exp$Symbol

PRJEB23709.Inflammatory.Sig.Exp <- PRJEB23709.Inflammatory.Sig.Exp[,-1]

PRJEB23709.Inflammatory.Sig.Exp <- as.matrix(PRJEB23709.Inflammatory.Sig.Exp)

PRJEB23709.Inflammatory.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.Inflammatory.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.Inflammatory.Sig.Exp.Set <- fit(PRJEB23709.Inflammatory.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.Inflammatory.Sig.Exp.Set <- predict(predictor_PRJEB23709.Inflammatory.Sig.Exp.Set, PRJEB23709.Inflammatory.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.Inflammatory.Sig.Exp), dist = "cor")

z_PRJEB23709_Inflammatory.Sig <- as.numeric(prediction_PRJEB23709.Inflammatory.Sig.Exp.Set@prediction[,'z'])

Inflammatory.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_Inflammatory.Sig)

colnames(Inflammatory.Sig.PRJEB23709) <- c('grp','res')

Inflammatory.Sig.PRJEB23709 <- as.data.frame(Inflammatory.Sig.PRJEB23709)


# Other_11 EMT.Sig


PRJEB23709.EMT.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% EMT.Sig,]

rownames(PRJEB23709.EMT.Sig.Exp) <- PRJEB23709.EMT.Sig.Exp$Symbol

PRJEB23709.EMT.Sig.Exp <- PRJEB23709.EMT.Sig.Exp[,-1]

PRJEB23709.EMT.Sig.Exp <- as.matrix(PRJEB23709.EMT.Sig.Exp)

PRJEB23709.EMT.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.EMT.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.EMT.Sig.Exp.Set <- fit(PRJEB23709.EMT.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.EMT.Sig.Exp.Set <- predict(predictor_PRJEB23709.EMT.Sig.Exp.Set, PRJEB23709.EMT.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.EMT.Sig.Exp), dist = "cor")

z_PRJEB23709_EMT.Sig <- as.numeric(prediction_PRJEB23709.EMT.Sig.Exp.Set@prediction[,'z'])

EMT.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_EMT.Sig)

colnames(EMT.Sig.PRJEB23709) <- c('grp','res')

EMT.Sig.PRJEB23709 <- as.data.frame(EMT.Sig.PRJEB23709)


# Other_12 Blood.Sig

PRJEB23709.Blood.Sig.Exp <- CC_73samples_GE_matrix[CC_73samples_GE_matrix$Symbol %in% Blood.Sig,]

rownames(PRJEB23709.Blood.Sig.Exp) <- PRJEB23709.Blood.Sig.Exp$Symbol

PRJEB23709.Blood.Sig.Exp <- PRJEB23709.Blood.Sig.Exp[,-1]

PRJEB23709.Blood.Sig.Exp <- as.matrix(PRJEB23709.Blood.Sig.Exp)

PRJEB23709.Blood.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(PRJEB23709.Blood.Sig.Exp),phenoData=CC_73samples_phenoData)

predictor_PRJEB23709.Blood.Sig.Exp.Set <- fit(PRJEB23709.Blood.Sig.Exp.Set, method = "welch.test") 

prediction_PRJEB23709.Blood.Sig.Exp.Set <- predict(predictor_PRJEB23709.Blood.Sig.Exp.Set, PRJEB23709.Blood.Sig.Exp.Set,"Nonresponder", ngenes=nrow(PRJEB23709.Blood.Sig.Exp), dist = "cor")

z_PRJEB23709_Blood.Sig <- as.numeric(prediction_PRJEB23709.Blood.Sig.Exp.Set@prediction[,'z'])

Blood.Sig.PRJEB23709 <- cbind(out_CC_73samples_V5,z_PRJEB23709_Blood.Sig)

colnames(Blood.Sig.PRJEB23709) <- c('grp','res')

Blood.Sig.PRJEB23709 <- as.data.frame(Blood.Sig.PRJEB23709)


# Multiplot - PRJEB23709

PRJEB23709_Data <- list(ImmuneCells.Sig=ImmuneCells.Sig.PRJEB23709, IMPRES.Sig=IMPRES.Sig.PRJEB23709, IRG.Sig=IRG.Sig.PRJEB23709, IPRES.Sig=IPRES.Sig.PRJEB23709, 
T.cell.inflamed.Sig=T.cell.inflamed.Sig.PRJEB23709, LRRC15.CAF.Sig=LRRC15.CAF.Sig.PRJEB23709, EMT.Sig=EMT.Sig.PRJEB23709, Inflammatory.Sig=Inflammatory.Sig.PRJEB23709, 
Blood.Sig=Blood.Sig.PRJEB23709, IFNG.Sig=IFNG.Sig.PRJEB23709, CD8.Sig=CD8.Sig.PRJEB23709, PDL1.Sig=PDL1.Sig.PRJEB23709, CRMA.Sig=CRMA.Sig.PRJEB23709)
                  
       
tiff("./output/No3_PRJEB23709_73samples.tiff", width = 7600, height = 7600, units = "px", res = 710)
 
p <- rocplot.multiple.V3(PRJEB23709_Data, title = "", p.value = FALSE)
print(p)

dev.off()



### No4 - MGSP

### ImmuneCells.Sig

ImmuneCells.Sig.MGSP <- Test_NatMed_103samples 


### Other signatures for comparison in MGSP

# Other_1 IFNG.Sig

MGSP.IFNG.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% IFNG.Sig,]

rownames(MGSP.IFNG.Sig.Exp) <- MGSP.IFNG.Sig.Exp$Symbol

MGSP.IFNG.Sig.Exp <- MGSP.IFNG.Sig.Exp[,-1]

MGSP.IFNG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.IFNG.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.IFNG.Sig.Exp.Set <- fit(MGSP.IFNG.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.IFNG.Sig.Exp.Set <- predict(predictor_MGSP.IFNG.Sig.Exp.Set, MGSP.IFNG.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.IFNG.Sig.Exp), dist = "cor")

z_MGSP_IFNG.Sig <- as.numeric(prediction_MGSP.IFNG.Sig.Exp.Set@prediction[,'z'])

IFNG.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_IFNG.Sig)

colnames(IFNG.Sig.MGSP) <- c('grp','res')

IFNG.Sig.MGSP <- as.data.frame(IFNG.Sig.MGSP)


# Other_2 CD8.Sig

MGSP.CD8.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% CD8.Sig,]

rownames(MGSP.CD8.Sig.Exp) <- MGSP.CD8.Sig.Exp$Symbol

MGSP.CD8.Sig.Exp <- MGSP.CD8.Sig.Exp[,-1]

MGSP.CD8.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.CD8.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.CD8.Sig.Exp.Set <- fit(MGSP.CD8.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.CD8.Sig.Exp.Set <- predict(predictor_MGSP.CD8.Sig.Exp.Set, MGSP.CD8.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.CD8.Sig.Exp), dist = "cor")

z_MGSP_CD8.Sig <- as.numeric(prediction_MGSP.CD8.Sig.Exp.Set@prediction[,'z'])

CD8.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_CD8.Sig)

colnames(CD8.Sig.MGSP) <- c('grp','res')

CD8.Sig.MGSP <- as.data.frame(CD8.Sig.MGSP)

CD8.Sig.MGSP <- CD8.Sig.MGSP[!is.na(CD8.Sig.MGSP$res),]


# Other_3 PDL1.Sig

MGSP.PDL1.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% PDL1.Sig,]

rownames(MGSP.PDL1.Sig.Exp) <- MGSP.PDL1.Sig.Exp$Symbol

MGSP.PDL1.Sig.Exp <- MGSP.PDL1.Sig.Exp[,-1]

MGSP.PDL1.Sig.Exp <- as.matrix(MGSP.PDL1.Sig.Exp)

MGSP.PDL1.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.PDL1.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.PDL1.Sig.Exp.Set <- fit(MGSP.PDL1.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.PDL1.Sig.Exp.Set <- predict(predictor_MGSP.PDL1.Sig.Exp.Set, MGSP.PDL1.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.PDL1.Sig.Exp), dist = "cor")

z_MGSP_PDL1.Sig <- as.numeric(prediction_MGSP.PDL1.Sig.Exp.Set@prediction[,'z'])

PDL1.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_PDL1.Sig)

colnames(PDL1.Sig.MGSP) <- c('grp','res')

PDL1.Sig.MGSP <- as.data.frame(PDL1.Sig.MGSP)

PDL1.Sig.MGSP <- PDL1.Sig.MGSP[!is.na(PDL1.Sig.MGSP$res),]


# Other_4 CRMA.Sig 

MGSP.CRMA.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% CRMA.Sig,]

rownames(MGSP.CRMA.Sig.Exp) <- MGSP.CRMA.Sig.Exp$Symbol

MGSP.CRMA.Sig.Exp <- MGSP.CRMA.Sig.Exp[,-1]

MGSP.CRMA.Sig.Exp <- as.matrix(MGSP.CRMA.Sig.Exp)

MGSP.CRMA.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.CRMA.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.CRMA.Sig.Exp.Set <- fit(MGSP.CRMA.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.CRMA.Sig.Exp.Set <- predict(predictor_MGSP.CRMA.Sig.Exp.Set, MGSP.CRMA.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.CRMA.Sig.Exp), dist = "cor")

z_MGSP_CRMA.Sig <- as.numeric(prediction_MGSP.CRMA.Sig.Exp.Set@prediction[,'z'])

CRMA.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_CRMA.Sig)

colnames(CRMA.Sig.MGSP) <- c('grp','res')

CRMA.Sig.MGSP <- as.data.frame(CRMA.Sig.MGSP)

CRMA.Sig.MGSP <- CRMA.Sig.MGSP[!is.na(CRMA.Sig.MGSP$res),]


# Other_5 IMPRES.Sig 

MGSP.IMPRES.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% IMPRES.Sig,]

rownames(MGSP.IMPRES.Sig.Exp) <- MGSP.IMPRES.Sig.Exp$Symbol

MGSP.IMPRES.Sig.Exp <- MGSP.IMPRES.Sig.Exp[,-1]

MGSP.IMPRES.Sig.Exp <- as.matrix(MGSP.IMPRES.Sig.Exp)

MGSP.IMPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.IMPRES.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.IMPRES.Sig.Exp.Set <- fit(MGSP.IMPRES.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.IMPRES.Sig.Exp.Set <- predict(predictor_MGSP.IMPRES.Sig.Exp.Set, MGSP.IMPRES.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.IMPRES.Sig.Exp), dist = "cor")

z_MGSP_IMPRES.Sig <- as.numeric(prediction_MGSP.IMPRES.Sig.Exp.Set@prediction[,'z'])

IMPRES.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_IMPRES.Sig)

colnames(IMPRES.Sig.MGSP) <- c('grp','res')

IMPRES.Sig.MGSP <- as.data.frame(IMPRES.Sig.MGSP)


# Other_6 IRG.Sig

MGSP.IRG.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% IRG.Sig,]

rownames(MGSP.IRG.Sig.Exp) <- MGSP.IRG.Sig.Exp$Symbol

MGSP.IRG.Sig.Exp <- MGSP.IRG.Sig.Exp[,-1]

MGSP.IRG.Sig.Exp <- as.matrix(MGSP.IRG.Sig.Exp)

MGSP.IRG.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.IRG.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.IRG.Sig.Exp.Set <- fit(MGSP.IRG.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.IRG.Sig.Exp.Set <- predict(predictor_MGSP.IRG.Sig.Exp.Set, MGSP.IRG.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.IRG.Sig.Exp), dist = "cor")

z_MGSP_IRG.Sig <- as.numeric(prediction_MGSP.IRG.Sig.Exp.Set@prediction[,'z'])

IRG.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_IRG.Sig)

colnames(IRG.Sig.MGSP) <- c('grp','res')

IRG.Sig.MGSP <- as.data.frame(IRG.Sig.MGSP)


# Other_7 LRRC15.CAF.Sig

MGSP.LRRC15.CAF.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% LRRC15.CAF.Sig,]

rownames(MGSP.LRRC15.CAF.Sig.Exp) <- MGSP.LRRC15.CAF.Sig.Exp$Symbol

MGSP.LRRC15.CAF.Sig.Exp <- MGSP.LRRC15.CAF.Sig.Exp[,-1]

MGSP.LRRC15.CAF.Sig.Exp <- as.matrix(MGSP.LRRC15.CAF.Sig.Exp)

MGSP.LRRC15.CAF.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.LRRC15.CAF.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.LRRC15.CAF.Sig.Exp.Set <- fit(MGSP.LRRC15.CAF.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.LRRC15.CAF.Sig.Exp.Set <- predict(predictor_MGSP.LRRC15.CAF.Sig.Exp.Set, MGSP.LRRC15.CAF.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.LRRC15.CAF.Sig.Exp), dist = "cor")

z_MGSP_LRRC15.CAF.Sig <- as.numeric(prediction_MGSP.LRRC15.CAF.Sig.Exp.Set@prediction[,'z'])

LRRC15.CAF.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_LRRC15.CAF.Sig)

colnames(LRRC15.CAF.Sig.MGSP) <- c('grp','res')

LRRC15.CAF.Sig.MGSP <- as.data.frame(LRRC15.CAF.Sig.MGSP)


# Other_8_T.cell.inflamed.Sig


MGSP.T.cell.inflamed.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% T.cell.inflamed.Sig,]

rownames(MGSP.T.cell.inflamed.Sig.Exp) <- MGSP.T.cell.inflamed.Sig.Exp$Symbol

MGSP.T.cell.inflamed.Sig.Exp <- MGSP.T.cell.inflamed.Sig.Exp[,-1]

MGSP.T.cell.inflamed.Sig.Exp <- as.matrix(MGSP.T.cell.inflamed.Sig.Exp)

MGSP.T.cell.inflamed.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.T.cell.inflamed.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.T.cell.inflamed.Sig.Exp.Set <- fit(MGSP.T.cell.inflamed.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.T.cell.inflamed.Sig.Exp.Set <- predict(predictor_MGSP.T.cell.inflamed.Sig.Exp.Set, MGSP.T.cell.inflamed.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.T.cell.inflamed.Sig.Exp), dist = "cor")

z_MGSP_T.cell.inflamed.Sig <- as.numeric(prediction_MGSP.T.cell.inflamed.Sig.Exp.Set@prediction[,'z'])

T.cell.inflamed.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_T.cell.inflamed.Sig)

colnames(T.cell.inflamed.Sig.MGSP) <- c('grp','res')

T.cell.inflamed.Sig.MGSP <- as.data.frame(T.cell.inflamed.Sig.MGSP)


# Other_9 IPRES.Sig

MGSP.IPRES.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% IPRES.Sig,]

rownames(MGSP.IPRES.Sig.Exp) <- MGSP.IPRES.Sig.Exp$Symbol

MGSP.IPRES.Sig.Exp <- MGSP.IPRES.Sig.Exp[,-1]

MGSP.IPRES.Sig.Exp <- as.matrix(MGSP.IPRES.Sig.Exp)

MGSP.IPRES.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.IPRES.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.IPRES.Sig.Exp.Set <- fit(MGSP.IPRES.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.IPRES.Sig.Exp.Set <- predict(predictor_MGSP.IPRES.Sig.Exp.Set, MGSP.IPRES.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.IPRES.Sig.Exp), dist = "cor")

z_MGSP_IPRES.Sig <- as.numeric(prediction_MGSP.IPRES.Sig.Exp.Set@prediction[,'z'])

IPRES.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_IPRES.Sig)

colnames(IPRES.Sig.MGSP) <- c('grp','res')

IPRES.Sig.MGSP <- as.data.frame(IPRES.Sig.MGSP)


# Other_10 Inflammatory.Sig

MGSP.Inflammatory.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% Inflammatory.Sig,]

rownames(MGSP.Inflammatory.Sig.Exp) <- MGSP.Inflammatory.Sig.Exp$Symbol

MGSP.Inflammatory.Sig.Exp <- MGSP.Inflammatory.Sig.Exp[,-1]

MGSP.Inflammatory.Sig.Exp <- as.matrix(MGSP.Inflammatory.Sig.Exp)

MGSP.Inflammatory.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.Inflammatory.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.Inflammatory.Sig.Exp.Set <- fit(MGSP.Inflammatory.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.Inflammatory.Sig.Exp.Set <- predict(predictor_MGSP.Inflammatory.Sig.Exp.Set, MGSP.Inflammatory.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.Inflammatory.Sig.Exp), dist = "cor")

z_MGSP_Inflammatory.Sig <- as.numeric(prediction_MGSP.Inflammatory.Sig.Exp.Set@prediction[,'z'])

Inflammatory.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_Inflammatory.Sig)

colnames(Inflammatory.Sig.MGSP) <- c('grp','res')

Inflammatory.Sig.MGSP <- as.data.frame(Inflammatory.Sig.MGSP)


# Other_11 EMT.Sig


MGSP.EMT.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% EMT.Sig,]

rownames(MGSP.EMT.Sig.Exp) <- MGSP.EMT.Sig.Exp$Symbol

MGSP.EMT.Sig.Exp <- MGSP.EMT.Sig.Exp[,-1]

MGSP.EMT.Sig.Exp <- as.matrix(MGSP.EMT.Sig.Exp)

MGSP.EMT.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.EMT.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.EMT.Sig.Exp.Set <- fit(MGSP.EMT.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.EMT.Sig.Exp.Set <- predict(predictor_MGSP.EMT.Sig.Exp.Set, MGSP.EMT.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.EMT.Sig.Exp), dist = "cor")

z_MGSP_EMT.Sig <- as.numeric(prediction_MGSP.EMT.Sig.Exp.Set@prediction[,'z'])

EMT.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_EMT.Sig)

colnames(EMT.Sig.MGSP) <- c('grp','res')

EMT.Sig.MGSP <- as.data.frame(EMT.Sig.MGSP)


# Other_12 Blood.Sig

MGSP.Blood.Sig.Exp <- NatMed_103samples_GE_matrix[NatMed_103samples_GE_matrix$Symbol %in% Blood.Sig,]

rownames(MGSP.Blood.Sig.Exp) <- MGSP.Blood.Sig.Exp$Symbol

MGSP.Blood.Sig.Exp <- MGSP.Blood.Sig.Exp[,-1]

MGSP.Blood.Sig.Exp <- as.matrix(MGSP.Blood.Sig.Exp)

MGSP.Blood.Sig.Exp.Set <- ExpressionSet(assayData=as.matrix(MGSP.Blood.Sig.Exp),phenoData=NatMed_103samples_phenoData)

predictor_MGSP.Blood.Sig.Exp.Set <- fit(MGSP.Blood.Sig.Exp.Set, method = "welch.test") 

prediction_MGSP.Blood.Sig.Exp.Set <- predict(predictor_MGSP.Blood.Sig.Exp.Set, MGSP.Blood.Sig.Exp.Set,"Progressor", ngenes=nrow(MGSP.Blood.Sig.Exp), dist = "cor")

z_MGSP_Blood.Sig <- as.numeric(prediction_MGSP.Blood.Sig.Exp.Set@prediction[,'z'])

Blood.Sig.MGSP <- cbind(out_NatMed_103samples,z_MGSP_Blood.Sig)

colnames(Blood.Sig.MGSP) <- c('grp','res')

Blood.Sig.MGSP <- as.data.frame(Blood.Sig.MGSP)


# Multiplot - MGSP

MGSP_Data <- list(ImmuneCells.Sig=ImmuneCells.Sig.MGSP, IMPRES.Sig=IMPRES.Sig.MGSP, IRG.Sig=IRG.Sig.MGSP, IPRES.Sig=IPRES.Sig.MGSP, 
T.cell.inflamed.Sig=T.cell.inflamed.Sig.MGSP, LRRC15.CAF.Sig=LRRC15.CAF.Sig.MGSP, EMT.Sig=EMT.Sig.MGSP, Inflammatory.Sig=Inflammatory.Sig.MGSP, 
Blood.Sig=Blood.Sig.MGSP, IFNG.Sig=IFNG.Sig.MGSP, CD8.Sig=CD8.Sig.MGSP, PDL1.Sig=PDL1.Sig.MGSP, CRMA.Sig=CRMA.Sig.MGSP)
                  
       
tiff("./output/No4_MGSP_103samples.tiff", width = 7600, height = 7600, units = "px", res = 710)
 
p <- rocplot.multiple.V3(MGSP_Data, title = "", p.value = FALSE)
print(p)

dev.off()



#########################################################################################################################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################
#########################################################################################################################################################################################
