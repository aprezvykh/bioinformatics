library(heatmap.2)
library(edgeR)
library(org.Dm.eg.db)
library(gplots)
library(plyr)
library(ggplot2)
### Statistical analysis
directory <- '~/Fly memory project/experimental_multimap/all//'
setwd(directory)
sampleFiles <- grep('fly',list.files(directory),value=TRUE)
sampleCondition <- c('control', 'control', 
                     'cross', 'cross', 
                     'stress', 'stress', 
                     'stress_24', 'stress_24')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)


y <- readDGE(files = sampleFiles, group = sampleCondition, labels = sampleFiles)
readqual <- as.data.frame(tail(y$counts, 3))

y <- calcNormFactors(y, method = "TMM", logratioTrim = TRUE)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

#MDS PLOT

pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition])
legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
dev.off()
keep <- rowSums(cpm(y) > 1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
et <- exactTest(y)
et_annot <- as.data.frame(et$table)

y$genes$Symbol <- mapIds(org.Dm.eg.db, 
                          keys=row.names(y), 
                          column="SYMBOL", 
                          keytype="FLYBASE",
                          multiVals="first")
y$genes$Name <- mapIds(org.Dm.eg.db, 
                          keys=row.names(y), 
                          column="GENENAME", 
                          keytype="FLYBASE",
                          multiVals="first")
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


## TOP 100 PVALUE
logCPM <- cpm(y, prior.count=2, log=TRUE, normalized.lib.sizes=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$PValue)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 Pvalue Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Transcripts differential
          expression, p < 0.05")
dev.off()

### TOP 100 LOGFC GENES 
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$logFC)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 LogFC Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Transcripts differential
          expression, p < 0.05")
dev.off()

### SEARCH AND PLOT!
let <- c("fruitless")
logCPM <- NULL
logCPM <- cpm(y, prior.count=2, log=TRUE)
nColCount <- ncol(logCPM)
logCPM <- as.data.frame(logCPM)
logCPM$Name <- mapIds(org.Dm.eg.db, 
                      keys=row.names(logCPM), 
                      column="GENENAME", 
                      keytype="FLYBASE",
                      multiVals="first")

rownames(logCPM) <- make.names(y$genes$Symbol, unique=TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM)[nColCount+1] <- c("Name")
sub <- logCPM[grepl(paste(let), logCPM$Name, ignore.case = TRUE),]
sub$Name <- NULL
sub <- t(scale(t(sub)))
pdf(file = paste(let,"_query_heatmap.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
heatmap.2(sub, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()


### MULTIPLE BOXPLOTS OF INTERESTING GENES
logdf <- as.data.frame(sub)
for (i in seq(1:nrow(logdf))){
  png(file = paste(rownames(logdf[i,]), "_CPM.png", sep=""))
  boxplot.default(logdf[i,],outline = TRUE,  main = paste(rownames(logdf[i,])))
  dev.off()
}

### Multiple MDPlots

for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}

pdf(file = "MDplot_common.pdf", width = 12, height = 17, family = "Helvetica")
plotMD(tr, values=c(1,-1), col=c("red","blue"),
       legend="topright")
dev.off()


### High Expressed Genes
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature")
cpm <-cpm(y)
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ]
cpm <- as.data.frame(cpm(y))
cpm$rowsum <- rowSums(cpm)
topcpm <- cpm[order(cpm$rowsum, decreasing = TRUE),]
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[1:100,]
topcpm$rowsum <- NULL
topcpm <- as.data.frame(topcpm)
colnames(topcpm) <- paste(y$samples$group, 1:2, sep="-")
topcpm$Symbol <- mapIds(org.Dm.eg.db, 
                      keys=row.names(topcpm), 
                      column="SYMBOL", 
                      keytype="FLYBASE",
                      multiVals="first")
rownames(topcpm) <- make.names(topcpm$Symbol, unique=TRUE)
topcpm$Symbol <- NULL
topcpm <- t(scale(t(topcpm)))
pdf(file = "Top 100 high expressed genes.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(topcpm, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Highest expressed genes")
dev.off()


