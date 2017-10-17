zz <- file("error.log", open="wt")
sink(zz, type="message")
library(beepr)
library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(org.Mm.eg.db)
library(gplots)
library(pheatmap)
library(xlsx)
library(gage)
library(gageData)
library(topGO)
library(ggplot2)
library(ReactomePA)
library(reshape)
require(clusterProfiler)
require(reactome.db)
library(data.table)
library(pathview)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(Matrix)
library(PANTHER.db)
library(GO.db)
library(Hmisc)
library(DESeq2)
library(checkmate)
### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO

heatmaps <- TRUE
custom_heatmap <- FALSE
custom_genes_plots <- FALSE
fisherGO <- FALSE
analyze_all_samples <- TRUE
disease_association <- TRUE
kegg_plots <- TRUE
panther_analysis <- TRUE
deseq2_part <- TRUE
qlm_test <- FALSE

pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5
gr_control <- c("tg_1")
gr_case <- c("tg_3")
gs_size <- 10
diseases_set <- 50
number_of_kegg_plots <- 50
go_terms_set <- 50
pathways_set <- 30
genes_in_term <- 3


### Statistical analysis
col.pan <- colorpanel(100, "blue", "white", "red")
directory <- '~/GitHub/counts/ALS Mice/experimental/'
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

stattest <- paste(gr_control, gr_case, sep = "-")
directory <- '~/GitHub/counts/ALS Mice/experimental/results/'
setwd(directory)

if (analyze_all_samples == FALSE){
dir.create(stattest)
setwd(stattest)
} else if (analyze_all_samples == TRUE){
  dir.create("all")
  setwd("all")
}

readqual <- as.data.frame(tail(y$counts, 5))
libsize <- as.data.frame(t(y$samples$lib.size))
names(libsize) <- names(readqual)
readqual <- rbind(libsize, readqual)
rownames(readqual[1,]) <- c("lib_size")


readqual[nrow(readqual)+1,] <- (readqual[2,]/readqual[1,])*100
readqual[nrow(readqual)+1,] <- (readqual[3,]/readqual[1,])*100

qual <- readqual[7:8,]
rownames(qual) <- c("% of no feature", "% of ambigous")

m1 <- melt(qual[1,])
m2 <- melt(qual[2,])
pdf(file = "No feature.pdf", width = 12, height = 17, family = "Helvetica")
ggplot(m1) + aes(x = variable, y = value) + geom_bar(stat = "identity")
dev.off()

pdf(file = "Ambigous.pdf", width = 12, height = 17, family = "Helvetica")
ggplot(m1) + aes(x = variable, y = value) + geom_bar(stat = "identity")
dev.off()

row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )




### QLM TESTS!
if (qlm_test == TRUE){
      a <- DGEList(counts=y, group = sampleTable$condition)
      cpm <- cpm(y)
      cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ]
      cpm <- as.data.frame(cpm(y))
      cpm$rowsum <- rowSums(cpm)
      keep <- rowSums(cpm > cpm_cutoff) >= ncol(sampleTable)
      a <- a[keep, , keep.lib.sizes=FALSE]
      a <- calcNormFactors(a)
      design <- model.matrix(~sampleTable$condition)
      a <- estimateDisp(a,design)
      fit <- glmQLFit(a,design)
      qlf <- glmQLFTest(fit,coef=2)
      et_annot <- as.data.frame(qlf$table)
      et_annot_non_filtered <- as.data.frame(qlf$table)

} else if (qlm_test == FALSE){

      normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
      CountsTable <- as.data.frame(y$counts)
      raw_counts <- as.data.frame(y$counts)
      y <- estimateCommonDisp(y)
      y <- estimateTagwiseDisp(y)
      keep <- rowSums(cpm(y) > cpm_cutoff) >= ncol(sampleTable)
      cpm <- cpm(y)
      cpm <- as.data.frame(cpm(y))
      cpm$rowsum <- rowSums(cpm)
      y <- y[keep, , keep.lib.sizes=FALSE]
      normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
      logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
      et <- exactTest(y)
      top <- as.data.frame(topTags(et))
      et_annot <- as.data.frame(et$table)
      et_annot_non_filtered <- as.data.frame(et$table)
}
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

et_annot$GOID <-     mapIds(org.Mm.eg.db, 
                         keys=row.names(et_annot), 
                         column="GO", 
                         keytype="ENSEMBL",
                         multiVals="first")

et_annot$term <- mapIds(GO.db, 
                       keys=et_annot$GOID, 
                       column="TERM", 
                       keytype="GOID",
                       multiVals="first")

et_annot$term <- as.character(et_annot$term)

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
### FILTRATION

et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))

### TESTING A HYPOTESIS
counts_control <- CountsTable[,grep(gr_control, names(CountsTable))]
counts_case <- CountsTable[,grep(gr_case, names(CountsTable))]
counts_control$rowsum.control <- rowSums(counts_control)
counts_case$rowsum.case <- rowSums(counts_case)
diff <- data.frame(counts_control$rowsum.control, counts_case$rowsum.case)
rownames(diff) <- rownames(CountsTable) 

fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
c <- grep(paste(lfgene), rownames(diff))
dfc <- as.data.frame(diff[c,])
lfgenefc <- fc[4,1]
stat <- dfc$counts_control.rowsum.control > dfc$counts_case.rowsum.case

if (stat == TRUE & lfgenefc < 0){
  print("No Correction Needed!")
} else {
  et_annot$logFC <- et_annot$logFC*(-1)
}



### Distribution
dir.create("distributions")
setwd("distributions")

png(file = "logCPM distribution, filtered.png")
hist(et_annot$logCPM, main = "logCPM distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM")
dev.off()
png(file = "p-value distribution, filtered.png")
hist(et_annot$PValue, main = "p-value distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value")
dev.off()
png(file = "LogFC distribution, filtered.png")
hist(et_annot$logFC, main = "logFC distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC")
dev.off()

png(file = "logCPM distribution, nonfiltered.png")
hist(et_annot_non_filtered$logCPM, main = "logCPM distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM")
dev.off()
png(file = "p-value distribution, nonfiltered.png")
hist(et_annot_non_filtered$PValue, main = "p-value distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value")
dev.off()
png(file = "LogFC distribution, nonfiltered.png")
hist(et_annot_non_filtered$logFC, main = "logFC distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC")
dev.off()

setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
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


write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)
# write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table, logCPM>1", append = TRUE)



### MY GO TESTS
tab <- as.data.frame(table(unlist(et_annot$GOID)))
rownames(tab) <- tab$Var1
my_go <- data.frame()
for (f in tab$Var1){
  sub <- grep(paste(f), et_annot$GOID, ignore.case = TRUE)
  st <- et_annot[sub,]
  lfc <- mean(st$logFC)
  df <- data.frame(f, lfc, nrow(st))
  my_go <- rbind(df, my_go)
}

rownames(my_go) <- my_go$f
my_go$term <- mapIds(GO.db, 
                         keys=rownames(my_go), 
                         column="TERM", 
                         keytype="GOID",
                         multiVals="first")
my_go$term <- as.character(my_go$term)
names(my_go) <- c("term", "logFC", "genes in term", "name of term")

for_hm <- my_go
my_go <- my_go[order(my_go$`genes in term`, decreasing = TRUE),]
write.xlsx(my_go, file = "My GO.xlsx", sheetName = "GO classes", append = TRUE)


### GO HEATMAP
cpmfgh <- as.data.frame(logCPM)

cpmfgh$GOID <-            mapIds(org.Mm.eg.db, 
                                 keys=row.names(cpmfgh), 
                                 column="GO", 
                                 keytype="ENSEMBL",
                                 multiVals="first")


tab <- as.data.frame(table(unlist(cpmfgh$GOID)))
rownames(tab) <- tab$Var1
my_go <- data.frame()

for (f in rownames(tab)){
  sub <- grep(paste(f), cpmfgh$GOID, ignore.case = TRUE)
  st <- cpmfgh[sub,]
  st$GOID <- NULL
  a <- as.data.frame(colMeans(st))
  a <- t(a)
  a <- as.data.frame(a)
  rownames(a) <- paste(f)
  my_go <- rbind(a, my_go)
}


for_hm <- for_hm[order(for_hm$logFC),]
for_hm <- as.data.frame(subset(for_hm, for_hm$`genes in term` > genes_in_term))
for_hm <- as.data.frame(subset(for_hm, for_hm$logFC > logfchigh_cutoff | for_hm$logFC < logfclow_cutoff))

go_for_heatmap <- data.frame()

for (f in for_hm$term){
  a <- grep(paste(f), rownames(my_go))
  d <- my_go[a,]
  go_for_heatmap <- rbind(d, go_for_heatmap)
  
}

go_for_heatmap$term <- mapIds(GO.db, 
                        keys=row.names(go_for_heatmap), 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")




rownames(go_for_heatmap) <- go_for_heatmap$term
go_for_heatmap$term <- NULL

go_for_heatmap <- t(scale(t(go_for_heatmap)))
pdf(file = "GO heatmap", width = 12, height = 17, family = "Helvetica")
heatmap.2(go_for_heatmap, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=0.7, cexCol=0.7, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,15))

dev.off()

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

if (fisherGO == TRUE){

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
}



## CORELLATION MARIX WITH DESEQ2
correl <- cpm(y)
x <- cor(correl)
png(file = "Corellation matrix.png")
heatmap.2(x, margins = c(10,10))
dev.off()
        
### MDS PLOT
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
# pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
png(file = "MDSPlot.png")
plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition])
legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
dev.off()



### VOLCANO PLOT
allgenes <- nrow(et_annot_non_filtered)
et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > logfchigh_cutoff & et_annot_non_filtered$PValue < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 12, height = 17, family = "Helvetica")
g = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()


### REACTOME PART
dfa <- as.character(et_annot$entrez)
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

#HIGH
df_high <- et_annot_high$entrez
x <- enrichPathway(gene=df_high, organism = "mouse", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "mouse", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()



###KEGG expression profile (without lfc, but with generatio)

kk <- enrichKEGG(gene = df_high, organism = "mmu", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

kk <- enrichKEGG(gene = df_low, organism = "mmu", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()


# TOP 100 PVALUE GENES
if (heatmaps == TRUE){

logCPM <- as.data.frame(logCPM)
rownames(logCPM) <- make.names(y$genes$Symbol, unique = TRUE)
# colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$PValue)
logCPMpval <- logCPM[o[1:100],]
logCPMpval <- t(scale(t(logCPMpval)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 Pvalue Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPMpval, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top FDR genes, p < 0.05")
dev.off()

# TOP 100 LOGFC GENES
o <- order(et$table$logFC)
logCPMfc <- logCPM[o[1:100],]
logCPMfc <- t(scale(t(logCPMfc)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 logFC Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPMfc, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top Log2FoldChange genes, p < 0.05")
dev.off()


### High Expressed Genes
topcpm <- cpm[order(cpm$rowsum, decreasing = TRUE),]
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[1:100,]
topcpm$rowsum <- NULL
topcpm <- as.data.frame(topcpm)
# colnames(topcpm) <- paste(y$samples$group, 1:2, sep="-")
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

}



### Multiple MDPlots
dir.create("mdplots")
setwd("mdplots")
for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}
pdf(file = "MDplot_common.pdf", width = 12, height = 17, family = "Helvetica")
plotMD(y, values=c(1,-1), col=c("red","blue"),
       legend="topright")
dev.off()
setwd(stattest)
### SEARCH AND PLOT!
if (custom_heatmap == TRUE) {
logCPM$Name <- y$genes$Name
let <- c("caspase", 
         "apoptosis",
         "neural",
         "neuron",
         "death",
         "mitochondrial",
         "ATP")
for (f in let){
  sub <- logCPM[grepl(paste(f), logCPM$Name),]
  sub$Name <- NULL
  sub <- t(scale(t(sub)))
  pdf(file = paste(f,"_query_heatmap.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
  heatmap.2(sub, col=col.pan, Rowv=TRUE, scale="column",
            trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
  dev.off()
  if (custom_genes_plots == TRUE){
  logdf <- as.data.frame(sub)
  for (i in seq(1:nrow(logdf))){
    png(file = paste(rownames(logdf[i,]), "_CPM.png", sep=""))
    boxplot.default(logdf[i,],outline = TRUE,  main = paste(rownames(logdf[i,])))
    dev.off()}
  }
  
}
}

# KEGG PLOTS
setwd(directory)
setwd(stattest)
data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggreswr <- as.data.frame(keggres)
write.xlsx(keggreswr, file = "KEGG pathview.xlsx", sheetName = "KEGG", append = TRUE)

if (kegg_plots == TRUE){
keggrespathways = data.frame(id=rownames(keggres$greater),
  keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=number_of_kegg_plots) %>% 
  .$id %>% 
  as.character()
keggresids = substr(keggrespathways, start=1, stop=8)

plot_pathway = function(pid){
         pathview(gene.data=foldchanges, 
         pathway.id=pid, 
         species="mmu", 
         new.signature=FALSE)
}
detach("package:dplyr", unload=TRUE)
dir.create("kegg")
setwd("kegg")
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
}

pathview(gene.data=foldchanges, 
         pathway.id="mmu03050", 
         species="mmu", 
         new.signature=FALSE)

setwd(directory)
### DISEASE ASSOCIATION
setwd(stattest)
if (disease_association == TRUE){
et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))

table <- read.delim(file = '~/GitHub/data/curated_gene_disease_associations.tsv')
table <- as.data.frame(table)

diseases_up <- data.frame()
diseases_down <- data.frame()
for (f in et_annot_high$symbol){
  sub <- NULL
  sub <- table[grepl(paste(f), table$geneSymbol, ignore.case = TRUE),]
  sub <- sub[order(sub$score, decreasing = TRUE),]
  sub <- sub[seq(1:5),]
  sub <- as.data.frame(sub$diseaseName)
  sub <- transpose(sub)
  sub$gene <- paste(f)
  diseases_up <- rbind(sub, diseases_up)
}
# write.xlsx(diseases_up, file = "Diseases association by Disgenet.xlsx", sheetName = "upreg", append = TRUE)

for (f in et_annot_low$symbol){
  sub <- NULL
  sub <- table[grepl(paste(f), table$geneSymbol, ignore.case = TRUE),]
  sub <- sub[order(sub$score, decreasing = TRUE),]
  sub <- sub[seq(1:5),]
  sub <- as.data.frame(sub$diseaseName)
  sub <- transpose(sub)
  sub$gene <- paste(f)
  diseases_down<- rbind(sub, diseases_down)
}
# write.xlsx(diseases_down, file = "Diseases association by Disgenet.xlsx", sheetName = "downreg", append = TRUE)


up <- as.data.frame(table(unlist(diseases_up)))
up <- up[order(up$Freq, decreasing = TRUE),]
up <- up[seq(1:diseases_set),]
names(up) <- c("Disease", "Frequency")
write.xlsx(up, file = "Top Diseases by Disgenet.xlsx", sheetName = "upreg", append = TRUE)

down <- as.data.frame(table(unlist(diseases_down)))
down <- down[order(down$Freq, decreasing = TRUE),]
down <- down[seq(1:diseases_set),]
names(down) <- c("Disease", "Frequency")
write.xlsx(down, file = "Top Diseases by Disgenet.xlsx", sheetName = "downreg", append = TRUE)
}




### PANTHER.DB

if (panther_analysis == TRUE){
pan_up <- et_annot_high
pan_down <- et_annot_low

pan_up$goslim <- mapIds(PANTHER.db, 
                        keys=et_annot_high$entrez, 
                        column="GOSLIM_ID", 
                        keytype="ENTREZ",
                        multiVals="first")

pan_up$pathway <- mapIds(PANTHER.db, 
                        keys=et_annot_high$entrez, 
                        column="PATHWAY_ID", 
                        keytype="ENTREZ",
                        multiVals="first")
pan_down$goslim <- mapIds(PANTHER.db, 
                        keys=et_annot_low$entrez, 
                        column="GOSLIM_ID", 
                        keytype="ENTREZ",
                        multiVals="first")

pan_down$pathway <- mapIds(PANTHER.db, 
                        keys=et_annot_low$entrez, 
                        column="PATHWAY_ID", 
                        keytype="ENTREZ",
                        multiVals="first")


go_pan_up <- as.data.frame(table(unlist(pan_up$goslim)))
go_pan_up <- go_pan_up[order(go_pan_up$Freq, decreasing = TRUE),]
go_pan_up <- go_pan_up[seq(1:go_terms_set),]
names(go_pan_up) <- c("GO term", "Frequency")

go_pan_down <- as.data.frame(table(unlist(pan_down$goslim)))
go_pan_down <- go_pan_down[order(go_pan_down$Freq, decreasing = TRUE),]
go_pan_down <- go_pan_down[seq(1:go_terms_set),]
names(go_pan_down) <- c("GO term", "Frequency")

pth_pan_up <- as.data.frame(table(unlist(pan_up$pathway)))
pth_pan_up <- pth_pan_up[order(pth_pan_up$Freq, decreasing = TRUE),]
pth_pan_up <- pth_pan_up[seq(1:pathways_set),]
names(pth_pan_up) <- c("Pathway ID", "Frequency")

pth_pan_down <- as.data.frame(table(unlist(pan_down$pathway)))
pth_pan_down <- pth_pan_down[order(pth_pan_down$Freq, decreasing = TRUE),]
pth_pan_down <- pth_pan_down[seq(1:pathways_set),]
names(pth_pan_down) <- c("Pathway ID", "Frequency")

write.xlsx(go_pan_up, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(go_pan_down, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)

write.xlsx(pth_pan_up, file = "Top Pathways by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(pth_pan_down, file = "Top Pathways by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)
}

### DESEQ2 PART (ADDITIONAL)
directory <- '~/GitHub/counts/ALS Mice/experimental/'
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res <- results(dds, tidy = FALSE )
rld<- rlogTransformation(dds, blind=TRUE)

setwd('~/GitHub/counts/ALS Mice/experimental/results/')
setwd(stattest)
png(file = "PCAPlot(DESeq2).png")
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

pdf(file = "MAplot(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotMA(dds,ylim=c(-10,10),main='DESeq2')
dev.off()

pdf(file = "Dispestimate(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotDispEsts(dds)
dev.off()

pdf(file = "Sparsity(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotSparsity(dds)
dev.off()

sink(type="message")
close(zz)
