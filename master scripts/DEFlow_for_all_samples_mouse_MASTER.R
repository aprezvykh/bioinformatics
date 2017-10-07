source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
install.packages("Rcpp")
library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(org.Mm.eg.db)
library(gplots)
library(plyr)
library(dplyr)
library(pheatmap)
library(xlsx)
library(gage)
library(gageData)
library(topGO)
library(ggplot2)
### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO
analyze_all_samples <- FALSE
pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5
gr_control <- c("tg_1")
gr_case <- c("tg_2")

### Statistical analysis
directory <- '~/bioinformatics/counts/ALS Mice/experimental/'
setwd(directory)

if (analyze_all_samples == TRUE){
  sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
  sampleCondition <- c('control_early', 'control_early', 'control_early', 'control_early', 'control_early', 
                       'control_late', 'control_late', 'control_late', 'control_late', 'control_late', 
                       'tg_early', 'tg_early', 'tg_early', 'tg_early', 'tg_early', 
                       'tg_mid', 'tg_mid', 'tg_mid', 'tg_mid', 
                       'tg_late', 'tg_late', 'tg_late', 'tg_late', 'tg_late')
  
  sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
  y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
} else if (analyze_all_samples == FALSE){
        files_control <- grep(paste(gr_control),list.files(directory),value=TRUE)
        files_case <- grep(paste(gr_case),list.files(directory),value=TRUE)
        sampleFiles <- c(files_control, files_case)
        cond_control <- rep(paste(gr_control), length(files_control))
        cond_case <- rep(paste(gr_case), length(files_case))
        sampleCondition <- c(cond_control, cond_case)
        sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
        y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
}

readqual <- as.data.frame(tail(y$counts, 5))
libsize <- as.data.frame(t(y$samples$lib.size))
names(libsize) <- names(readqual)
readqual <- rbind(libsize, readqual)
rownames(readqual[1,]) <- c("lib_size")
normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
CountsTable <- as.data.frame(y$counts)
raw_counts <- as.data.frame(y$counts)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
top <- as.data.frame(topTags(et))
et_annot <- as.data.frame(et$table)
et_annot_non_filtered <- as.data.frame(et$table)

### ANNOTATE

y$genes$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
y$genes$Name <- mapIds(org.Mm.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="GENENAME", 
                         keytype="ENSEMBL",
                         multiVals="first")

et_annot$symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

et_annot$name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(et_annot), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")


et_annot$entrez <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot), 
                          column="ENTREZID", 
                          keytype="ENSEMBL",
                          multiVals="first")


et_annot_non_filtered$Symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")
et_annot_non_filtered$Symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="ENTREZID", 
                          keytype="ENSEMBL",
                          multiVals="first")

### TESTING A HYPOTESIS

fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
c <- grep(paste(lfgene), rownames(CountsTable))
dfc <- as.data.frame(CountsTable[c,])
dfc$meancontrol <- (dfc[1,1] + dfc[1,2])/2
dfc$meancase <- (dfc[1,3] + dfc[1,4])/2
if(dfc$meancase > dfc$meancontrol){
  et_annot$logFC <- et_annot$logFC*(-1)
}

### Simple summary
all <- nrow(raw_counts)
allpadj <- sum(et_annot$PValue < pvalue_cutoff, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > logfchigh_cutoff, na.rm=TRUE)
down <- sum(et_annot$logFC < logfclow_cutoff, na.rm=TRUE)
header <- c('all genes', 'mean of logCPM', 'padj<0,05', 'genes with > high', 'genes with < low')
meaning <- c(print(all), print(avg_cpm), print(allpadj), print(up), print(down))
df <- data.frame(header, meaning)


### ANNOTATE COUNTS
CountsTable$symbol <- mapIds(org.Mm.eg.db, 
                             keys=row.names(CountsTable), 
                             column="SYMBOL", 
                             keytype="ENSEMBL",
                             multiVals="first")

CountsTable$name <- mapIds(org.Mm.eg.db, 
                           keys=row.names(CountsTable), 
                           column="GENENAME", 
                           keytype="ENSEMBL",
                           multiVals="first")

CountsTable <- CountsTable[rownames(et_annot),]

### GO TESTS
data(go.sets.mm)
data(go.subs.mm)

foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE,set.size = c(10, 100), rank.test = TRUE)
lapply(gomfres, head)
gomfres <- as.data.frame(gomfres)
gomfres <- gomfres[complete.cases(gomfres), ]
write.xlsx(gomfres, file = "GO.xlsx", sheetName = "GO_MF", append = TRUE)

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- as.data.frame(gobpres)
gobpres <- gobpres[complete.cases(gobpres), ]
write.xlsx(gobpres, file = "GO.xlsx", sheetName = "GO_BP", append = TRUE)

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- as.data.frame(goccres)
goccres <- goccres[complete.cases(goccres), ]
write.xlsx(goccres, file = "GO.xlsx", sheetName = "GO_CC", append = TRUE)


### FISHER GO TESTS
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

GOFisherBP <- function(df, nodes, nrows, p){
  all_genes <- c(df$logFC)
  names(all_genes) <- rownames(df)
  go_data <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
                   p, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
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
                   p, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
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
                   p, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
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


### MDS PLOT
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition])
legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
dev.off()

### VOLCANO PLOT
allgenes <- nrow(et_annot_non_filtered)
et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > 1 & et_annot_non_filtered$PValue < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 12, height = 17, family = "Helvetica")
g = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()


et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))

write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)
write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table, logCPM>1", append = TRUE)

# TOP 100 PVALUE GENES
y$counts$Symbol <- NULL
y$counts$Symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(y$counts), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
rownames(logCPM) <- y$counts$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$PValue)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 Pvalue Heatmap_2.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top FDR genes, p < 0.05")
dev.off()

# TOP 100 LOGFC GENES
logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
rownames(logCPM) <- y$counts$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$logFC)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 logFC Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top Log2FoldChange genes, p < 0.05")
dev.off()


### SEARCH AND PLOT!
let <- c("sas")
logCPM <- NULL
normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
nColCount <- ncol(logCPM)
logCPM <- as.data.frame(logCPM)
logCPM$Name <- mapIds(org.Mm.eg.db, 
                          keys=row.names(logCPM), 
                          column="GENENAME", 
                          keytype="ENSEMBL",
                          multiVals="first")

rownames(logCPM) <- make.names(y$counts$Symbol, unique=TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM)[nColCount+1] <- c("Name")

sub <- logCPM[grepl(paste(let), logCPM$Name),]
sub$Name <- NULL
sub <- t(scale(t(sub)))
pdf(file = paste(let,"_query_heatmap.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
heatmap.2(sub, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()


### High Expressed Genes
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )
cpm <- cpm(y)
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ]
cpm <- as.data.frame(cpm(y))
cpm$rowsum <- rowSums(cpm)
topcpm <- cpm[order(cpm$rowsum, decreasing = TRUE),]
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[1:100,]
topcpm$rowsum <- NULL
topcpm <- as.data.frame(topcpm)
colnames(topcpm) <- paste(y$samples$group, 1:2, sep="-")
topcpm$Symbol <- mapIds(org.Mm.eg.db, 
                        keys=row.names(topcpm), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")
rownames(topcpm) <- make.names(topcpm$Symbol, unique=TRUE)
topcpm$Symbol <- NULL
topcpm <- t(scale(t(topcpm)))
pdf(file = "Top 100 high expressed genes.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(topcpm, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Highest expressed genes")
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
  dev.off()}

### DISEASE

