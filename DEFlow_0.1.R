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
biocLite("org.Mm.eg.db")
biocLite("IRanges")
install.packages("IRanges")
install.packages("ggplot2")
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
library(reshape2)
library(clusterProfiler)
### Statistical analysis (add 1 tg_mid)
directory <- 'C://Users/alexander/Dropbox/ALS-mice project/counts_trimmed_geo/counts_ens/5_ctrl_late_vs_ctrl_early//'
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('early', 'early', 'early', 'early', 'early',  
                     'late', 'late', 'late', 'late', 'late')

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)

## Simple tests ##
res<-results(dds)
allpadj <- sum(res$padj < 0.05, na.rm=TRUE)
padjcount <- sum(res$padj < 0.05, na.rm=TRUE)
res<-res[order(res$padj),]
resadj <- head(res, print(padjcount))
print(resadj)
head(res, print(padjcount))
mcols(res,use.names=TRUE)
# write.csv(as.data.frame(res),file='results.csv')
resadj <- resadj[order(resadj$log2FoldChange),]
logfcdown <- sum(resadj$log2FoldChange < -1, na.rm=TRUE)
logfcup <- sum(resadj$log2FoldChange > 1, na.rm=TRUE)
header <- c('padj<0,05', 'genes with logfcup', 'genes with logfcdown')
meaning <- c(print(allpadj), print(logfcup), print(logfcdown))
df <- data.frame(header, meaning)
write.csv(df, file = "result_stattests.csv")

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

resOrdered <- resadj[order(resadj$padj),]
resOrderedBM <- resOrdered[order(resOrdered$baseMean),]
resOrderedBM <- resOrderedBM[complete.cases(resOrderedBM), ]
write.csv(as.data.frame(resOrderedBM), file = "result_DEG_annotated.csv")

### GO ###
data(go.sets.mm)
data(go.subs.mm)
gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE)
lapply(gomfres, head)
gomfres <- gomfres[complete.cases(gomfres), ]
write.csv(gomfres, file = "MF.csv")


gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- gobpres[complete.cases(gobpres), ]
write.csv(gobpres, file = "BP.csv")

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- goccres[complete.cases(goccres), ]
write.csv(goccres, file = "CC.csv")

## TRANSFORM ### 
rld<- rlogTransformation(dds, blind=TRUE)
print(plotPCA(rld, intgroup=c('condition')))
dev.copy(png, 'pcaplot.png')
dev.off()
plotMA(dds,ylim=c(-10,10),main='DESeq2')
dev.copy(png, 'maplot.png')
dev.off()
### Genes clustering

topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ),50)
pheatmap( assay(rld)[topVarGenes,], scale="row",
          trace="column", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), fontsize = 9
)
dev.copy(png, 'topvargenes.png')
dev.off()
plotDispEsts(dds)
dev.copy(png, 'dispestimate.png')
dev.off()
plotSparsity(dds)
dev.copy(png, 'sparcity.png')
dev.off()

## MATRIX ### 
dev.off()
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(50)
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=' : '))

hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))

dev.copy(png, 'distmatrix.png')
dev.off()


## REACTOME
require(clusterProfiler)
require(reactome.db)
print(resOrderedBM)
typeof(resOrderedBM)
dfa <- as.character(resOrderedBM$entrez)
print(dfa)
summary(dfa)
class(dfa)
x = enrichPathway(gene=dfa, organism = "mouse", minGSSize=15, readable = TRUE )
head(as.data.frame(x))
dev.off()
par(mar=c(1,1,1,1))
dotplot(x, showCategory=5)
dev.copy(png, 'dotplot.png')
dev.off()
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7 )
dev.copy(png, 'enrichmap.png')
dev.off()
cnetplot(x, showCategory = "10", foldChange = foldchanges, categorySize="pvalue")
