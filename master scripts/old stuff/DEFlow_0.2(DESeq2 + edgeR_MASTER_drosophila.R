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
install.packages("sqldf")
install.packages("RH2")
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
library("topGO")
library("edgeR")
library("tcltk")
library("sqldf")
library("RH2")
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
base_mean_cutoff <- 20
logfcup_cutoff <- 0.3
logfcdown_cutoff <- -0.3
gs_size <- 15
base_mean_cutoff_value <- 5
hm_genes_count <- 100
cpm_cutoff <- 1

### Statistical analysis
directory <- '~/Fly memory project/experimental_multimap/K_vs_F24/'
setwd(directory)
sampleFiles <- grep('fly',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'stress_24', 'stress_24')
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
write.xlsx(df, file = "Results diffexpression.xlsx", sheetName = "Common info", append = TRUE)
# write.xlsx(res, file = "Results all.xlsx", sheetName = "All DEG no filter")
sumres <- summary(res)
### Annotation ###

columns(org.Dm.eg.db)
resadj$symbol <- mapIds(org.Dm.eg.db, 
                        keys=row.names(resadj), 
                        column="SYMBOL", 
                        keytype="FLYBASE",
                        multiVals="first")

resadj$entrez <- mapIds(org.Dm.eg.db, 
                        keys=row.names(resadj), 
                        column="ENTREZID", 
                        keytype="FLYBASE",
                        multiVals="first")

resadj$name =   mapIds(org.Dm.eg.db,
                       keys=row.names(resadj), 
                       column="GENENAME",
                       keytype="FLYBASE",
                       multiVals="first")

write.csv(resadj, file = "K_vs_F24.csv")

### BaseMean Filtering

resOrderedBM <- resadj[order(resadj$padj),]
resOrderedBM <- resOrderedBM[,colSums(is.na(resOrderedBM))<nrow(resOrderedBM)]
resOrderedBM <- resOrderedBM[complete.cases(resOrderedBM), ]
resHBM <- resOrderedBM
upper_bm_cutoff <- obm*base_mean_cutoff_value
resHBM <- as.data.frame(subset(resHBM, baseMean > upper_bm_cutoff))
write.xlsx(resOrderedBM, file = "Results diffexpression.xlsx", sheetName = "high expressed genes, logfc not considered", append = TRUE)
### Fold Change & BaseMean filtering
resOrderedBM <- resOrderedBM[order(resOrderedBM$log2FoldChange),]
write.xlsx(resOrderedBM, file = "Results diffexpression.xlsx", sheetName = "DEG padj <0.05", append = TRUE)
dfx <- as.data.frame(resOrderedBM)
dfx <- as.data.frame(subset(resOrderedBM, baseMean > base_mean_cutoff))
dfx <- as.data.frame(subset(dfx, log2FoldChange > logfcup_cutoff | log2FoldChange < logfcdown_cutoff))
write.xlsx(dfx, file = "Results diffexpression.xlsx", sheetName = "genes with Bm > 100 & -1 > log2FC > 1 ", append = TRUE)
resOrderedBM <- dfx


## TRANSFORM ###
rld<- rlogTransformation(dds, blind=TRUE)
pdf(file = "pcaplot.pdf", width = 12, height = 17, family = "Helvetica")
print(plotPCA(rld, intgroup=c('condition')))
dev.off()
pdf(file = "MAplot.pdf", width = 12, height = 17, family = "Helvetica")
plotMA(dds,ylim=c(-10,10),main='DESeq2')
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
resadj$threshold = as.factor(abs(resadj$log2FoldChange) > 2 & resadj$padj < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 12, height = 17, family = "Helvetica")
g = ggplot(data=resadj, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()

# function - gets raw counts by gene name
cfds <- function(gene_name){
  ens <- rownames(resc[grep(gene_name, resc$symbol, ignore.case=TRUE),])
  m <- match(ens, rownames(resc))
  m <- as.data.frame(m)
  return(plotCounts(dds, gene = ens, main = gene_name, transform = FALSE, returnData = TRUE, normalized = TRUE))
}

resc <- as.data.frame(resadj)

gene_list = NULL

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

### EDGER PART!

pval_cutoff <- 0.05
base_mean_cutoff <- 20
logfcup_cutoff <- 0.3
logfcdown_cutoff <- -0.3
gs_size <- 15
base_mean_cutoff_value <- 5
hm_genes_count <- 100
cpm_cutoff <- 1

### Statistical analysis
directory <- '~/Fly memory project/experimental_multimap/K_vs_F24/'
setwd(directory)
sampleFiles <- grep('fly',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 'stress_24', 'stress_24')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

setwd(directory)
y <- readDGE(files = sampleFiles, group = sampleCondition, labels = sampleFiles)
y <- DGEList(y, group=class)
normalized_lib_sizes <- calcNormFactors(y)
keep <- rowSums(cpm(y) > 0.5) >= ncol(sampleTable)/2
y <- y[keep, , keep.lib.sizes=FALSE]
log_cpm <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
CountsTable <- as.data.frame(y$counts)
raw_counts <- as.data.frame(y$counts)
y <- calcNormFactors(y, method = "TMM")
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
et_annot <- as.data.frame(et$table)
et_annot_non_filtered <- as.data.frame(et$table)


### ANNOTATE LogFC

columns(org.Dm.eg.db)
et_annot$symbol <- mapIds(org.Dm.eg.db, 
                          keys=row.names(et_annot), 
                          column="SYMBOL", 
                          keytype="FLYBASE",
                          multiVals="first")

et_annot$name <- mapIds(org.Dm.eg.db, 
                        keys=row.names(et_annot), 
                        column="GENENAME", 
                        keytype="FLYBASE",
                        multiVals="first")


et_annot$entrez <- mapIds(org.Dm.eg.db, 
                          keys=row.names(et_annot), 
                          column="ENTREZID", 
                          keytype="FLYBASE",
                          multiVals="first")


### ANNOTATE COUNTS
CountsTable$symbol <- mapIds(org.Dm.eg.db, 
                             keys=row.names(CountsTable), 
                             column="SYMBOL", 
                             keytype="FLYBASE",
                             multiVals="first")

CountsTable$name <- mapIds(org.Dm.eg.db, 
                           keys=row.names(CountsTable), 
                           column="GENENAME", 
                           keytype="FLYBASE",
                           multiVals="first")
CountsTable <- CountsTable[rownames(et_annot),]

### TESTING CONTROL/CASE STATEMENT
fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
c <- grep(paste(lfgene), rownames(CountsTable))
dfc <- as.data.frame(CountsTable[c,])
dfc$meancontrol <- (dfc[1,1] + dfc[1,2])/2
dfc$meancase <- (dfc[1,3] + dfc[1,4])/2
if(dfc$meancase > dfc$meancontrol){
  et_annot$logFC <- et_annot$logFC*(-1)
}

### GO ### 
et_annot <- et_annot[complete.cases(et_annot), ]
et_annot <- as.data.frame(subset(et_annot, PValue < 0.05))
et_annot <- as.data.frame(subset(et_annot, logFC > logfcup_cutoff | logFC < logfcdown_cutoff))
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

GOFisherBP <- function(df, nodes, nrows, p){
  all_genes <- c(df$logFC)
  names(all_genes) <- rownames(df)
  go_data <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
                   p, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                 ID = "ENSEMBL", nodeSize = nodes)
  go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  go_table <- GenTable(go_data, weightFisher = go_test,
                       orderBy = "weightFisher", ranksOf = "weightFisher",
                       topNodes = nrows)
  return(go_table)
}
GOFisherMF <- function(df, nodes, nrows, p){
  all_genes <- c(df$logFC)
  names(all_genes) <- rownames(df)
  go_data <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(s) s < 
                   p, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                 ID = "ENSEMBL", nodeSize = nodes)
  go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  go_table <- GenTable(go_data, weightFisher = go_test,
                       orderBy = "weightFisher", ranksOf = "weightFisher",
                       topNodes = nrows)
  return(go_table)
}
GOFisherCC <- function(df, nodes, nrows, p){
  all_genes <- c(df$logFC)
  names(all_genes) <- rownames(df)
  go_data <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(s) s < 
                   p, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                 ID = "ENSEMBL", nodeSize = nodes)
  go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
  go_table <- GenTable(go_data, weightFisher = go_test,
                       orderBy = "weightFisher", ranksOf = "weightFisher",
                       topNodes = nrows)
  return(go_table)
}

gobp_h_50 <- GOFisherBP(et_annot_high, 5, 50, 0.05)
gomf_h_50 <- GOFisherMF(et_annot_high, 5, 50, 0.05)
gocc_h_50 <- GOFisherCC(et_annot_high, 5, 50, 0.05)

gobp_l_50 <- GOFisherBP(et_annot_low, 5, 50, 0.05)
gomf_l_50 <- GOFisherMF(et_annot_low, 5, 50, 0.05)
gocc_l_50 <- GOFisherCC(et_annot_low, 5, 50, 0.05)

gobp_h_100 <- GOFisherBP(et_annot_high, 5, 100, 0.05)
gomf_h_100 <- GOFisherMF(et_annot_high, 5, 100, 0.05)
gocc_h_100 <- GOFisherCC(et_annot_high, 5, 100, 0.05)

gobp_l_100 <- GOFisherBP(et_annot_low, 5, 100, 0.05)
gomf_l_100 <- GOFisherMF(et_annot_low, 5, 100, 0.05)
gocc_l_100 <- GOFisherCC(et_annot_low, 5, 100, 0.05)

gobp_h_1000 <- GOFisherBP(et_annot_high, 5, 1000, 0.05)
gomf_h_1000 <- GOFisherMF(et_annot_high, 5, 1000, 0.05)
gocc_h_1000 <- GOFisherCC(et_annot_high, 5, 1000, 0.05)

gobp_l_1000 <- GOFisherBP(et_annot_low, 5, 1000, 0.05)
gomf_l_1000 <- GOFisherMF(et_annot_low, 5, 1000, 0.05)
gocc_l_1000 <- GOFisherCC(et_annot_low, 5, 1000, 0.05)


write.xlsx(gobp_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 50", append = TRUE)
write.xlsx(gomf_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 50", append = TRUE)
write.xlsx(gocc_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 50", append = TRUE)
write.xlsx(gobp_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(gomf_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(gocc_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 100", append = TRUE)
write.xlsx(gobp_h_1000, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 1000", append = TRUE)
write.xlsx(gomf_h_1000, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 1000", append = TRUE)
write.xlsx(gocc_h_1000, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 1000", append = TRUE)

write.xlsx(gobp_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 50", append = TRUE)
write.xlsx(gomf_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 50", append = TRUE)
write.xlsx(gocc_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 50", append = TRUE)
write.xlsx(gobp_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(gomf_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(gocc_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 100", append = TRUE)
write.xlsx(gobp_l_1000, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 1000", append = TRUE)
write.xlsx(gomf_l_1000, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 1000", append = TRUE)
write.xlsx(gocc_l_1000, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 1000", append = TRUE)


### KEGG
require(clusterProfiler)
require(reactome.db)
dfa <- as.character(et_annot$entrez)
foldchanges <- et_annot$logFC
x <- enrichPathway(gene=dfa, organism = "fly", minGSSize=5, readable = TRUE )
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


### Simple summary
all <- nrow(et_annot_non_filtered)
allpadj <- sum(et_annot$PValue < pval_cutoff, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > logfcup_cutoff, na.rm=TRUE)
down <- sum(et_annot$logFC < logfcdown_cutoff, na.rm=TRUE)
header <- c('all genes', 'mean of logCPM', 'padj<0,05', 'genes with > high', 'genes with < low')
meaning <- c(print(all), print(avg_cpm), print(allpadj), print(up), print(down))
df <- data.frame(header, meaning)

dev.off()

# NEW HEATMAP
y$genes$Symbol <- mapIds(org.Dm.eg.db, 
                        keys=row.names(y), 
                        column="SYMBOL", 
                        keytype="FLYBASE",
                        multiVals="first")


logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol 
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$PValue)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()
for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}


et_annot <- as.data.frame(subset(et_annot, PValue < 0.1))
et_annot <- as.data.frame(subset(et_annot, logFC > logfcup_cutoff | logFC < logfcdown_cutoff))

write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)
write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table, logCPM>1", append = TRUE)







