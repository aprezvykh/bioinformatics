library(AnnotationDbi)
library(org.Mm.eg.db)
library(Rcpp)
library(edgeR)
library(xlsx)

### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO
pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5

### Statistical analysis
directory <- '~/GitHub/counts/ALS Mice/motoneurons laser dissected/'
setwd(directory)
sampleFiles <- grep('moto',list.files(directory),value=TRUE)
sampleCondition <- c('1', '1')
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

write.csv(cpm, file = "moto_expression_profile.csv")
