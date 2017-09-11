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
biocLite("org.Mm.eg.db")
biocLite("IRanges")
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
library("DESeq2")
library('RColorBrewer')
library('gplots')
library("pheatmap")
library('vsn')
library( "genefilter" )
library("AnnotationDbi")
library("org.Mm.eg.db")
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
library(xlsx)
library("calibrate")
library("foreach")

### Parameters.
pval_cutoff <- 0.05
base_mean_cutoff <- 100 
logfcup_cutoff <- 1
logfcdown_cutoff <- -1
gs_size <- 15
base_mean_cutoff_value <- 5

### Statistical analysis
directory <- '~/counts_ens/5_ctrl_late_vs_ctrl_early//'
setwd('~/diffexp_reports/')
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'control', 'control', 'control',
                     'tg', 'tg', 'tg', 'tg', 'tg')
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
logfcdown <- sum(resadj$log2FoldChange < logfcdown_cutoff, na.rm=TRUE)
logfcup <- sum(resadj$log2FoldChange > logfcup_cutoff, na.rm=TRUE)
obm <- sum(res$baseMean)/allgenes

header <- c('all genes', 'mean of overall baseMean', 'padj<0,05', 'genes with logfcup', 'genes with logfcdown')
meaning <- c(print(allgenes), print(obm), print(allpadj), print(logfcup), print(logfcdown))
df <- data.frame(header, meaning)
write.xlsx(df, file = "report.xlsx", sheetName = "Common info")
sumres <- summary(res)
### Annotation ###

columns(org.Mm.eg.db)
resadj$symbol <- mapIds(org.Mm.eg.db, 
                        keys=row.names(resadj), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")

resadj$entrez <- mapIds(org.Mm.eg.db, 
                        keys=row.names(resadj), 
                        column="ENTREZID", 
                        keytype="ENSEMBL",
                        multiVals="first")

resadj$name =   mapIds(org.Mm.eg.db,
                       keys=row.names(resadj), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")


### BaseMean Filtering

resOrderedBM <- resadj[order(resadj$padj),]
resOrderedBM <- resOrderedBM[,colSums(is.na(resOrderedBM))<nrow(resOrderedBM)]
resOrderedBM <- resOrderedBM[complete.cases(resOrderedBM), ]
resHBM <- resOrderedBM
upper_bm_cutoff <- obm*base_mean_cutoff_value
resHBM <- as.data.frame(subset(resHBM, baseMean > upper_bm_cutoff))
write.xlsx(resHBM, file = "Top_basemean.xlsx", sheetName = "Top basemean padj <0.05")

### Fold Change & BaseMean filtering
resOrderedBM <- resOrderedBM[order(resOrderedBM$log2FoldChange),]
write.xlsx(resOrderedBM, file = "DEG_all.xlsx", sheetName = "DEG padj <0.05")
dfx <- as.data.frame(resOrderedBM)
dfx <- as.data.frame(subset(resOrderedBM, baseMean > base_mean_cutoff))
dfx <- as.data.frame(subset(dfx, log2FoldChange > logfcup_cutoff | log2FoldChange < logfcdown_cutoff))
write.xlsx(dfx, file = "DEG_filtered.xlsx", sheetName = "DEG padj <0.05 filtered")
resOrderedBM <- dfx

## GO ###

data(go.sets.mm)
data(go.subs.mm)
foldchanges = resOrderedBM$log2FoldChange
names(foldchanges) = resOrderedBM$entrez

gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE,set.size = c(10, 100), rank.test = TRUE)
lapply(gomfres, head)
gomfres <- as.data.frame(gomfres)
gomfres <- gomfres[complete.cases(gomfres), ]
write.xlsx(gomfres, file = "GO_MF.xlsx", sheetName = "GO_MF")

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- as.data.frame(gobpres)
gobpres <- gobpres[complete.cases(gobpres), ]
write.xlsx(gobpres, file = "GO_BP.xlsx", sheetName = "GO_BP")

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- as.data.frame(goccres)
goccres <- goccres[complete.cases(goccres), ]

write.xlsx(goccres, file = "GO_CC.xlsx", sheetName = "GO_CC")

## TRANSFORM ### 

rld<- rlogTransformation(dds, blind=TRUE)
pdf(file = "pcaplot.pdf", width = 12, height = 17, family = "Helvetica")
print(plotPCA(rld, intgroup=c('condition')))
dev.off()
pdf(file = "MAplot.pdf", width = 12, height = 17, family = "Helvetica")
plotMA(dds,ylim=c(-10,10),main='DESeq2')
dev.off()

### Genes heatmap

select <- order((resOrderedBM$log2FoldChange), decreasing=FALSE)[1:50]
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(50)
pdf(file = "topvargenes.pdf", width = 12, height = 17, family = "Helvetica")

pheatmap(assay(dds, normalized = TRUE)[select,],
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
# heatmap.2(mat, Rowv=as.dendrogram(hc),
#          symm=TRUE, trace='none',
#          col = rev(hmcol), margin=c(13, 13))
pheatmap(mat)

dev.off()


## REACTOME ## 

require(clusterProfiler)
require(reactome.db)
dfa <- as.character(resOrderedBM$entrez)
x <- enrichPathway(gene=dfa, organism = "mouse", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges, categorySize="pvalue", showCategory = 10)
dev.off()

## high expressed genes expression profile - reactome profiling
## of genes with baseMean = avg(baseMean)*base_mean_cutoff_value

foldchanges2 = resHBM$log2FoldChange
dfb <- as.character(resHBM$entrez)
y <- enrichPathway(gene=dfb, organism = "mouse", minGSSize=gs_size, readable = TRUE )
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_GEP.pdf", width = 12, height = 17, family = "Helvetica")
barplot(y, showCategory=30,  font.size = 9)
dev.off()

pdf(file = "enrichmap_GEP.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(y, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot_GEP.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges2, categorySize="pvalue", showCategory = 10)
dev.off()

### KEGG ###  
data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggres <- as.data.frame(keggres)
keggres <- keggres[complete.cases(keggres), ]
write.xlsx(keggres, file = "KEGG.xlsx", sheetName = "KEGG")

###KEGG expression profile (without lfc, but with generatio)

kk <- enrichKEGG(gene = dfa, organism = "mmu", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG_DEP.xlsx", sheetName = "KEGG")
pdf(file = "KEGG_enrichment_dotplot.pdf", width = 12, height = 17, family = "Helvetica")
dotplot(kk)
dev.off()

### Volcano Plot
resadj <- as.data.frame(resadj)
resadj$threshold = as.factor(abs(resadj$log2FoldChange) > 2 & resadj$padj < 0.05/allgenes)
pdf(file = "Volcano plot_late_vs_early.pdf", width = 12, height = 17, family = "Helvetica")
g = ggplot(data=resadj, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()

### Counts plot. Enter Refseq ID in "gene_name" variable!
resc <- as.data.frame(resadj)


gene_name <- "Sprr1a"

ens <- rownames(resc[grep(gene_name, resc$symbol, ignore.case=TRUE),])
m <- match(ens, rownames(resc))
m <- as.data.frame(m)

if(nrow(m)=0){

  print("Gene not found, change your query!")
  stop()
}

if(nrow(m)>1){
  
  print("Too much overlap, try another gene ID!")
  stop()

} else {
  pdf(file = paste(gene_name, ".pdf", sep=""), width = 12, height = 17, family = "Helvetica")
  plotCounts(dds, gene = ens, main = gene_name, transform = FALSE)
  dev.off()
  
}

