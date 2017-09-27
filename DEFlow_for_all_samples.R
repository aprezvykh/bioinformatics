library(edgeR)
library(org.Mm.eg.db)

pval_cutoff <- 0.05
base_mean_cutoff <- 20
logfcup_cutoff <- 1
logfcdown_cutoff <- -1
gs_size <- 15
base_mean_cutoff_value <- 5
hm_genes_count <- 100
cpm_cutoff <- 1

### Statistical analysis
directory <- '~/counts_ens/'
setwd('~/counts_ens/')
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control_early', 'control_early', 'control_early', 'control_early', 'control_early',
                     'control_late', 'control_late', 'control_late', 'control_late', 'control_late', 
                     'tg_early', 'tg_early', 'tg_early', 'tg_early', 'tg_early', 
                     'tg_mid', 'tg_mid', 'tg_mid', 'tg_mid',   
                     'tg_late', 'tg_late', 'tg_late', 'tg_late', 'tg_late')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleFiles, group = sampleCondition, labels = sampleFiles)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
et_annot <- as.data.frame(et$table)
y$genes$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(y), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
y$genes$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(y), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")
fit <- glmQLFit(y, robust=TRUE)
tr <- glmTreat(fit, lfc=log2(1.5))
logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$Symbol 
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:50],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

let <- c(" actin ")
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
sub <- logCPM[grep(paste(let), l,]
sub <- t(scale(t(sub)))
pdf(file = "paste(let).pdf", width = 12, height = 17, family = "Helvetica")
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

logdf <- as.data.frame(logCPM)
for (i in seq(1:nrow(logdf))){
  png(file = paste(rownames(logdf[i,]), "_CPM.png", sep=""))
  boxplot.default(logdf[i,],outline = TRUE,  main = paste(rownames(logdf[i,])))
  dev.off()
  