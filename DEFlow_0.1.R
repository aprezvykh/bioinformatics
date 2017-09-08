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

### Parameters
pval_cutoff <- 0.1
base_mean_cutoff <- 100 
logfcup_cutoff <- 1
logfcdown_cutoff <- -1
gs_size <- 15

### Statistical analysis
directory <- 'C://Users//alexander/Dropbox//ALS-mice project//counts_trimmed_geo//counts_ens//5_ctrl_late_vs_ctrl_early///'
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'control', 'control', 'control', 
                     'late', 'late', 'late', 'late', 'late')

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)

## Simple tests ##
res<-results(dds)
allpadj <- sum(res$padj < pval_cutoff, na.rm=TRUE)
padjcount <- sum(res$padj < pval_cutoff, na.rm=TRUE)
res<-res[order(res$padj),]
resadj <- head(res, print(padjcount))
print(resadj)
head(res, print(padjcount))
mcols(res,use.names=TRUE)
resadj <- resadj[order(resadj$log2FoldChange),]
ddsh <- estimateSizeFactors(dds)
allgenes <- nrow(dds)
logfcdown <- sum(resadj$log2FoldChange < logfcdown_cutoff, na.rm=TRUE)
logfcup <- sum(resadj$log2FoldChange > logfcup_cutoff, na.rm=TRUE)
header <- c('all genes','padj<0,05', 'genes with logfcup', 'genes with logfcdown')
meaning <- c(print(allgenes), print(allpadj), print(logfcup), print(logfcdown))
df <- data.frame(header, meaning)
write.xlsx(df, file = "report.xlsx", sheetName = "Common info")

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

### Genes clustering BY rld (not commited to logfc, basemean and pvalue filtration!
### Use to your own risk of unsignificant and low-expressed genes

pdf(file = "topvargenes.pdf", width = 12, height = 17, family = "Helvetica")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ),20)
pheatmap( assay(rld)[topVarGenes,], scale="row",
          trace="column", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), fontsize = 9
)

dev.off()
## Estimate & sparsity

pdf(file = "dispestimate.pdf", width = 12, height = 17, family = "Helvetica")
plotDispEsts(dds)
dev.off()
pdf(file = "sparsity.pdf", width = 12, height = 17, family = "Helvetica")
plotSparsity(dds)
dev.off()

## MATRIX ### 
pdf(file = "matrix.pdf", width = 12, height = 17, family = "Helvetica")
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(50)
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=' : '))

hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))

dev.off()


## REACTOME ## 

require(clusterProfiler)
require(reactome.db)
print(resOrderedBM)
dfa <- as.character(resOrderedBM$entrez)
x <- enrichPathway(gene=dfa, organism = "mouse", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()


par(mar=c(1,1,1,1))
pdf(file = "dotplot.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges, categorySize="pvalue", showCategory = 10)
dev.off()


### KEGG ### 

data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggres <- as.data.frame(keggres)
keggres <- keggres[complete.cases(keggres), ]
write.xlsx(keggres, file = "KEGG.xlsx", sheetName = "KEGG")


### Heatmap by significant genes
