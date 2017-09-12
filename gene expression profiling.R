library("dplyr")
library("DESeq2")
library("AnnotationDbi")
library("org.Mm.eg.db")

pval_cutoff <- 0.05
base_mean_cutoff_value <- 1000

### Statistical analysis
directory <- '~/motoneurons compare/Microglia/'
setwd('~/diffexp_reports/')
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'control',
                     'late', 'late', 'late')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)

## Simple tests ##
res <- results(dds, tidy = FALSE )
allpadj <- sum(res$padj < pval_cutoff, na.rm=TRUE)
res <- res[order(res$padj),]
resadj <- head(res, print(allpadj))
mcols(res,use.names=TRUE)
ddsh <- estimateSizeFactors(dds)
allgenes <- nrow(dds)
obm <- sum(res$baseMean)/allgenes
columns(org.Mm.eg.db)
res$symbol <- mapIds(org.Mm.eg.db, 
                        keys=row.names(res), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")

resj$entrez <- mapIds(org.Mm.eg.db, 
                        keys=row.names(res), 
                        column="ENTREZID", 
                        keytype="ENSEMBL",
                        multiVals="first")

res$name =   mapIds(org.Mm.eg.db,
                       keys=row.names(res), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")

res <- res[,colSums(is.na(res))<nrow(res)]
res <- res[complete.cases(res), ]
upper_bm_cutoff <- obm*base_mean_cutoff_value
resHBM <- as.data.frame(subset(res, baseMean > upper_bm_cutoff))
df <- as.data.frame(resHBM)
