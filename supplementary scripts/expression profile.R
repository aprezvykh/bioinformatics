library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(org.Mm.eg.db)
library(gplots)
library(plyr)
library(dplyr)
library(pheatmap)
library(xlsx)
library(gage)
library(gageData)
library(topGO)
library(ggplot2)
### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO
pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5

### Statistical analysis
directory <- '~/GitHub/counts/ALS Mice/microglia/'
setwd(directory)
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('1', '1', '1', '1', '1', '1')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
cpm <- as.data.frame(cpm(y))
cpm$sum <- rowSums(cpm)
### ANNOTATE

cpm$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(cpm), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
cpm$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(cpm), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

write.csv(cpm, file = "microglia.csv")
