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
detach("package:tidyselect", unload=TRUE)
detach("package:pathview", unload=TRUE)
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
library(psych)
library(xlsx)
### Statistical analysis (add 1 tg_mid)
directory <- 'C://Users/rezvykh/Dropbox/ALS-mice project/counts_trimmed_geo/counts_ens/3_mid_tg_vs_ctrl_tg///'
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'control', 'control', 'control', 
                     'tg', 'tg', 'tg', 'tg')
                  
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
write.csv(df, file = "1_result_stattests.csv")

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
write.csv(as.data.frame(resOrderedBM), file = "result_DEG_annotated.csv")
lapply(resOrderedBM, head)
resOrderedBM <- resOrderedBM[complete.cases(resOrderedBM),]

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





rld<- rlogTransformation(dds, blind=TRUE)
vsd<-varianceStabilizingTransformation(dds, blind=TRUE)
ntd <- normTransform(dds)
print(plotPCA(rld, intgroup=c('condition')))
plotMA(dds,ylim=c(-10,10),main='DESeq2')
meanSdPlot(assay(ntd))
### variance graph ###
print(pca)
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))


### Dist Matrix
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition,sampleFiles , sep=' : '))

hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col = rev(hmcol), margin=c(13, 13))

### Dispersion plot
plotDispEsts(dds)
plotSparsity(dds)

### Heatmaps ###

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
select <- order((resadj$log2FoldChange), decreasing=TRUE)[1:50]
hmcol<- colorRampPalette(brewer.pal(11, 'RdYlBu'))(50)

pheatmap(assay(rld, normalized = TRUE)[select,],
         col = hmcol, cellwidth = 20, scale="row", fontsize = 8,
         clustering_distance_rows = 'euclidean',
         cluster_rows=F, cluster_cols=F,
         legend = TRUE, main = "Lens transcripts differential expression, p <0,05")

heatmap.2(assay(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = TRUE, Colv = FALSE, scale='none',
          dendrogram='row', trace='none', margin=c(4,4))

heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = TRUE, Colv = FALSE, scale='row',
          dendrogram='row', trace='none', margin=c(4,4))

heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = TRUE, Colv = FALSE, scale='row',
          dendrogram='row', trace='none', margin=c(4,4))


## Log variance
par(mai=ifelse(1:4 <= 2, par('mai'), 0))
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord]<150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c('blue', 'black')
matplot(px[ord], cbind(assay(vsd)[, 1], log2(px))[ord, ], type="l", lty=1, col=vstcol, xlab='n', ylab='f(n)')
legend('bottomright', legend = c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,'DESeq2_VST_and_log2.png')


## heatmap
resOrderedBM <- as.data.frame(resOrderedBM$log2FoldChange)


### Genes clustering

topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE ),50)
pheatmap( assay(rld)[topVarGenes,], scale="row",
          trace="column", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), fontsize = 9
)

a <- as.data.frame(assay(rld)[topVarGenes,])
write.csv(a, file = "tvg.csv")

plotCounts(dds, gene = "ENSMUSG00000021986", normalized = TRUE, transform = TRUE)
print(a)
a <- as.data.frame(assay(rld[topVarGenes,]))
row.names(a)

write.csv(res, file = "res.csv")

### KEGG ###  
data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 20)
foldchanges = resOrderedBM$log2FoldChange
names(foldchanges) = resOrderedBM$entrez
head(foldchanges)
# Results
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
names

lapply(keggres, head)
write.csv(as.data.frame(lap), file = "6_mid_tg_vs_late_tg.csv")

write.csv(head(keggres), file = "KEGGRESALL.csv")

keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)

tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))



### Reactome ### (foldchanges)
require(clusterProfiler)
require(reactome.db)
print(resOrderedBM)
typeof(resOrderedBM)
dfa <- as.character(resOrderedBM$entrez)
print(dfa)
summary(dfa)
class(dfa)
x = enrichPathway(gene=dfa, organism = "mouse", minGSSize=50, readable = TRUE )
write.csv(x, file = "GSEA.csv")
head(as.data.frame(x))
dotplot(x, showCategory=12)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7 )
