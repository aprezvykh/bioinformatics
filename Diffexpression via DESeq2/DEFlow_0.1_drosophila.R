### USE THIS ON YOUR OWN RISK!! ###
###### WITH BEST REGARDS #########
### ALEXANDER P. REZVYKH ##########
######## EIMB RAS #################

### Installing required packages ###
source('http://bioconductor.org/biocLite.R')
biocLite('DESeq2')
biocLite('vsn')
biocLite('topGO')
biocLite("gage")
biocLite("gageData")
biocLite("Biostrings")
biocLite("tibble")
biocLite("Matrix")
biocLite("OrgDB")
biocLite("pathview")
biocLite("ReactomePA")
biocLite("reactome.db")
biocLite("EGSEA")
biocLite("org.Dm.eg.db")
biocLite("IRanges")
biocLite("S4Vectors")
biocLite("GenomicRanges")
install.packages("xlsx")
install.packages("checkmate")
install.packages("ggplot2")
install.packages("gplots")
install.packages("plyr")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("WriteXLS")
install.packages("magrittr")
install.packages("gapminder")
install.packages("xlsx")
install.packages("functional")
install.packages("rJava")
install.packages("psych")
install.packages("IRanges")
install.packages("ggplot2")
install.packages("reticulate")
install.packages("calibrate")
install.packages("ggthemes")
install.packages("calibrate")
install.packages("ggthemes")
install.packages("tibble")
library("GenomicRanges")
library("tibble")
library("rJava")
library("pathview")
library("gage")
library("gageData")
library("gapminder")
library("magrittr")
library("plyr")
library("IRanges")
library("checkmate")
library("dplyr")
library("S4Vectors")
library("DESeq2")
library('RColorBrewer')
library('gplots')
library("pheatmap")
library('vsn')
library( "genefilter" )
library("AnnotationDbi")
library("org.Dm.eg.db")
library("xlsx")
library("ReactomePA")
library("functional")
library("DOSE")
library("EGSEA")
library("ggplot2")
library("reshape2")
library("clusterProfiler")
library("reticulate")
library("topGO")
library("xlsx")
library("calibrate")
library("foreach")
library("ggthemes")

### Parameters manual###
# pval_cutoff - cutoffs your pvalue. default is 0.05
# base_mean_cutoff - cutoffs low-expressed genes. default is 100
# logfcup/down - cutoffs genes without diffexpression default is 1/-1
# gs_size - minimal size of according to one reactome cluster genes. default is 15
# base_mean_cutoff_value - coefficient by which the average value multiplies average BaseMean
# it is required to the genes expression profile vizualization (GEP). default is 5
# gene_list - put in it genes RefSeq name, and output will be a diffexpression plot
# like PlotCounts DESeq2 functions.
# hm_genes_count - genes, that will be on a diffexpression heatmap

### Parameters
pval_cutoff <- 0.05
base_mean_cutoff <- 100 
logfcup_cutoff <- 1
logfcdown_cutoff <- -1
gs_size <- 15
base_mean_cutoff_value <- 5
hm_genes_count <- 100

### Statistical analysis
directory <- '~/Fly memory project/experimental/K_vs_F24/'
setwd('~/diffexp_reports/')
sampleFiles <- grep('fly',list.files(directory), value=TRUE)
sampleCondition <- c('24', '24', 'k', 'k')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)

## Simple tests ##
res <- results(dds, tidy = FALSE, pAdjustMethod = "bonferroni" )
allpadj <- sum(res$padj < pval_cutoff, na.rm=TRUE)
res <- res[order(res$padj),]
resadj <- head(res, print(allpadj))
mcols(res,use.names=TRUE)
ddsh <- estimateSizeFactors(dds)
allgenes <- nrow(dds)  
obm <- sum(res$baseMean)/allgenes
logfcdown <- sum(resadj$log2FoldChange < logfcdown_cutoff, na.rm=TRUE)
logfcup <- sum(resadj$log2FoldChange > logfcup_cutoff, na.rm=TRUE)
obm <- sum(res$baseMean)/allgenes
header <- c('all genes', 'mean of overall baseMean', 'padj<0,05', 'genes with logfcup', 'genes with logfcdown')
meaning <- c(print(allgenes), print(obm), print(allpadj), print(logfcup), print(logfcdown))
df <- data.frame(header, meaning)
write.xlsx(df, file = "Results all.xlsx", sheetName = "Common info", append = TRUE)
# write.xlsx(res, file = "Results all.xlsx", sheetName = "All DEG no filter")
### Annotation ###
aaa <- NULL
resadj <- as.data.frame(resadj)
columns(org.Dm.eg.db)

resadj$symbol <- mapIds(org.Dm.eg.db, 
                        keys=row.names(resadj), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")

resadj$entrez <- mapIds(org.Dm.eg.db, 
                        keys=row.names(resadj), 
                        column="ENTREZID", 
                        keytype="ENSEMBL",
                        multiVals="first")

resadj$name =   mapIds(org.Dm.eg.db,
                       keys=row.names(resadj), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")

res <- as.data.frame(res)


### BaseMean Filtering
aaa <- as.data.frame(resadj)
resOrderedBM <- resadj[order(resadj$padj),]
resOrderedBM <- resOrderedBM[,colSums(is.na(resOrderedBM))<nrow(resOrderedBM)]
resOrderedBM <- resOrderedBM[complete.cases(resOrderedBM), ]
resHBM <- resOrderedBM
upper_bm_cutoff <- obm*base_mean_cutoff_value
resHBM <- as.data.frame(subset(resHBM, baseMean > upper_bm_cutoff))

### Fold Change & BaseMean filtering
resOrderedBM <- resOrderedBM[order(resOrderedBM$log2FoldChange),]
write.xlsx(resOrderedBM, file = "Results all.xlsx", sheetName = "DEG padj <0.05", append = TRUE)
dfx <- as.data.frame(resOrderedBM)
dfx <- as.data.frame(subset(resOrderedBM, baseMean > base_mean_cutoff))
dfx <- as.data.frame(subset(dfx, log2FoldChange > logfcup_cutoff | log2FoldChange < logfcdown_cutoff))
write.xlsx(resHBM, file = "Results all.xlsx", sheetName = "genes with Bm>base_mean_cutoff ", append = TRUE)
resOrderedBM <- dfx


## TRANSFORM ### 

rld<- rlogTransformation(dds, blind=TRUE)
pdf(file = "pcaplot.pdf", width = 12, height = 17, family = "Helvetica")
print(plotPCA(rld, intgroup=c('condition')))
dev.off()
pdf(file = "MAplot.pdf", width = 12, height = 17, family = "Helvetica")
plotMA(dds,ylim=c(-10,10),main='DESeq2')
dev.off()

### Genes heatmap
hm <- as.data.frame(resOrderedBM)
select <- order((hm$baseMean), decreasing=TRUE)[1:hm_genes_count]
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(hm_genes_count)
pdf(file = "topvargenes.pdf", width = 12, height = 17, family = "Helvetica")
ass <-as.data.frame(assay(rld, normalized = TRUE)[select,])

hdf <- hm[1:hm_genes_count,]
rownames(ass) <- NULL
rownames(hdf) <- hdf$symbol
rownames(ass) <- rownames(hdf)

pheatmap(ass,
         cellwidth = 40, scale="row", fontsize = 10,
         clustering_distance_rows = 'euclidean',
         cluster_rows=F, cluster_cols=F,
         legend = TRUE, main = "Transcripts differential expression, p <0,05",
         clustering_method = "complete", show_rownames = T)

dev.off()
## Estimate & sparsity

pdf(file = "dispestimate.pdf", width = 12, height = 17, family = "Helvetica")
plotDispEsts(dds)
dev.off()
pdf(file = "sparsity.pdf", width = 12, height = 17, family = "Helvetica")
plotSparsity(dds)
dev.off()

## MATRIX ### 
par(mar=c(1,1,1,1))
pdf(file = "matrix.pdf", width = 12, height = 17, family = "Helvetica")
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(50)
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=' : '))

hc <- hclust(distsRL)
par(mar=c(1,1,1,1))
pheatmap(mat)

dev.off()

### Volcano Plot
res <- as.data.frame(res)
res$threshold = as.factor(abs(res$log2FoldChange) > 2 & res$padj < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 12, height = 17, family = "Helvetica")
g = ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()


### PLOTS
cfds <- function(gene_name){
  ens <- rownames(resc[grep(gene_name, resc$symbol, ignore.case=TRUE),])
  m <- match(ens, rownames(resc))
  m <- as.data.frame(m)
  return(plotCounts(dds, gene = ens, main = gene_name, transform = FALSE, returnData = TRUE, normalized = TRUE))
}

resc <- as.data.frame(resadj)

gene_list = list("TotM",
                 "onecut")

for (i in gene_list){
  u <- cfds(i)
  u <- as.data.frame(u)
  png(file = paste(i, ".png", sep=""))
  plot(x = u$condition, y = u$count, type = "a", main = paste(i), xlab = "Experimental group", ylab = "Normalized counts",
       col = "blue", lwd = 2, bg = "white",  lwd = 5, yline( 0, lwd=2, col=4))
  grid( col = "lightgray", lty = "dotted",
        lwd = par("lwd"), equilogs = TRUE)
  dev.off()
}
dev.off()