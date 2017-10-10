library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(org.Dm.eg.db)
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
### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO
custom_heatmap <- FALSE
custom_genes_plots <- FALSE
fisherGO <- FALSE
analyze_all_samples <- FALSE
kegg_plots <- TRUE

pvalue_cutoff <- 0.05
logfchigh_cutoff <- 0.3
logfclow_cutoff <- -0.3
cpm_cutoff <- 0.5
gr_control <- c("K")
gr_case <- c("N")
gs_size <- 10
diseases_set <- 50
number_of_kegg_plots <- 50
### Statistical analysis
directory <- '~/GitHub/counts/dr_multimap/'
setwd(directory)

if (analyze_all_samples == TRUE){
        sampleFiles <- grep('fly',list.files(directory),value=TRUE)
        sampleCondition <- c('K', 'K', 
                             'N', 'N',
                             'F', 'F', 
                             'M', 'M')
  
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
directory <- '~/GitHub/counts/dr_multimap/results/'
setwd(directory)
dir.create(stattest)
setwd(stattest)
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

normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
CountsTable <- as.data.frame(y$counts)
raw_counts <- as.data.frame(y$counts)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
keep <- rowSums(cpm(y) > cpm_cutoff) >= ncol(sampleTable)
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )
cpm <- cpm(y)
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ]
cpm <- as.data.frame(cpm(y))
cpm$rowsum <- rowSums(cpm)
y <- y[keep, , keep.lib.sizes=FALSE]
normalized_lib_sizes <- calcNormFactors(y, method = "TMM")
logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
et <- exactTest(y)
top <- as.data.frame(topTags(et))
et_annot <- as.data.frame(et$table)
et_annot_non_filtered <- as.data.frame(et$table)



### ANNOTATE

y$genes$Symbol <- mapIds(org.Dm.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="SYMBOL", 
                         keytype="FLYBASE",
                         multiVals="first")
y$genes$Name <- mapIds(org.Dm.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
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


et_annot_non_filtered$Symbol <- mapIds(org.Dm.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="SYMBOL", 
                          keytype="FLYBASE",
                          multiVals="first")
et_annot_non_filtered$Symbol <- mapIds(org.Dm.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="ENTREZID", 
                          keytype="FLYBASE",
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




### FISHER GO TESTS
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

if (fisherGO == TRUE){

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
}

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

### REACTOME PART
foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "fly", minGSSize=gs_size, readable = TRUE )
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
x <- enrichPathway(gene=df_high, organism = "fly", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "fly", minGSSize=gs_size, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()



###KEGG expression profile (without lfc, but with generatio)

kk <- enrichKEGG(gene = df_high, organism = "dme", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

kk <- enrichKEGG(gene = df_low, organism = "dme", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()


# TOP 100 PVALUE GENES
if (heatmaps == TRUE){

logCPM <- as.data.frame(logCPM)
rownames(logCPM) <- make.names(y$genes$Symbol, unique = TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
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

}



### Multiple MDPlots
dir.create("mdplots")
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
setwd(directory)
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
data(kegg.sets.dm)
data(sigmet.idx.dm)
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
         species="dme", 
         new.signature=FALSE)
}
detach("package:dplyr", unload=TRUE)
dir.create("kegg")
setwd("kegg")
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="dme"))
}

setwd(directory)

pathview(gene.data=foldchanges, 
         pathway.id="dme03050", 
         species="dme", 
         new.signature=FALSE)


