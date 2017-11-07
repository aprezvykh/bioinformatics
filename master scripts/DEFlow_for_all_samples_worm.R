library(DOSE)
library("GenomicRanges")
library(checkmate)
library(org.Ce.eg.db)
library(ENCODExplorer)
library(DT)
library(reshape2)
library(htmlwidgets)
library(visNetwork)
library(beepr)
library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(gplots)
library(lazyeval)
library(ggplot2)
library(pheatmap)
library(xlsx)
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
library(biomaRt)
library(devtools)
library(enrichR)
library(zoo)
library(rvest)
library(plyr)
library(AnnotationDbi)
library(gage)
library(gageData)
library(DESeq2)
heatmaps <- TRUE
custom_heatmap <- FALSE
custom_genes_plots <- FALSE
fisherGO <- TRUE
analyze_all_samples <- FALSE
kegg_plots <- TRUE
panther_analysis <- TRUE
deseq2_part <- TRUE
qlm_test <- TRUE
logging <- FALSE

### CONSTANTS BLOCK

pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5
gs_size <- 10
go_terms_set <- 50
pathways_set <- 30
genes_in_term <- 3
filter_thresh <- 5
baseMean_cutoff <- 1.5


### GROUPS. FIRST GROUP WILL BE USED AS CONTROL!
gr_control <- c("control")
gr_case <- c("case")


### BUILDING A SPECIFIC DESIGN TABLE
if (logging == TRUE){
  zz <- file("error.log", open="wt")
  sink(zz, type="message")
}

col.pan <- colorpanel(100, "blue", "white", "red")
directory <- '~/counts/worm_test/'
setwd(directory)

if (analyze_all_samples == TRUE){
  sampleFiles <- grep('worm',list.files(directory),value=TRUE)
  sampleCondition <- c('control', 'conrol', 'exp', 'exp')
  
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
directory <- '~/counts/worm_test/results'
setwd(directory)

if (analyze_all_samples == FALSE){
  dir.create(stattest)
  setwd(stattest)
} else if (analyze_all_samples == TRUE){
  dir.create("all")
  setwd("all")
}



# PLOTTING HTSEQ QUALITY BARPLOTS
dir.create("htseq-count quality plots")
setwd("htseq-count quality plots")
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
setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}

if (qlm_test == TRUE){ 
  a <- DGEList(counts=y, group = sampleTable$condition) 
  CountsTable <- as.data.frame(y$counts)
  cpm <- cpm(y) 
  cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
  cpm <- as.data.frame(cpm(y)) 
  cpm$rowsum <- rowSums(cpm) 
  keep <- rowSums(cpm > cpm_cutoff) >= ncol(sampleTable) 
  a <- a[keep, , keep.lib.sizes=FALSE] 
  a <- calcNormFactors(a, method = "TMM") 
  design <- model.matrix(~sampleTable$condition) 
  a <- estimateDisp(a,design) 
  fit <- glmQLFit(a,design, robust = TRUE) 
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  et_annot <- as.data.frame(qlf$table) 
  et_annot_non_filtered <- as.data.frame(qlf$table)
  top <- as.data.frame(topTags(qlf))
  logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
  et <- qlf
} else if (qlm_test == FALSE){ 
  design <- model.matrix(~sampleTable$condition) 
  normalized_lib_sizes <- calcNormFactors(y, method = "TMM") 
  CountsTable <- as.data.frame(y$counts) 
  raw_counts <- as.data.frame(y$counts) 
  y <- estimateCommonDisp(y) 
  y <- estimateTagwiseDisp(y) 
  y <- estimateDisp(y, design = design) 
  nf <- exactTest(y) 
  no_filtered <- as.data.frame(nf$table) 
  keep <- rowSums(cpm(y) > cpm_cutoff) >= ncol(sampleTable) 
  cpm <- cpm(y) 
  cpm <- as.data.frame(cpm(y)) 
  cpm$rowsum <- rowSums(cpm) 
  y <- y[keep, , keep.lib.sizes=FALSE] 
  logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
  et <- exactTest(y) 
  top <- as.data.frame(topTags(et)) 
  et_annot <- as.data.frame(et$table) 
  et_annot_non_filtered <- as.data.frame(et$table) 
}



y$genes$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")
y$genes$Name <- mapIds(org.Ce.eg.db, 
                       keys=row.names(et_annot_non_filtered), 
                       column="GENENAME", 
                       keytype="WORMBASE",
                       multiVals="first")

et_annot$symbol <- mapIds(org.Ce.eg.db, 
                          keys=row.names(et_annot), 
                          column="SYMBOL", 
                          keytype="WORMBASE",
                          multiVals="first")

et_annot$name <- mapIds(org.Ce.eg.db, 
                        keys=row.names(et_annot), 
                        column="GENENAME", 
                        keytype="WORMBASE",
                        multiVals="first")


et_annot$entrez <- mapIds(org.Ce.eg.db, 
                          keys=row.names(et_annot), 
                          column="ENTREZID", 
                          keytype="WORMBASE",
                          multiVals="first")
et_annot$uniprot <-mapIds(org.Ce.eg.db, 
                       keys=row.names(et_annot), 
                       column="UNIPROT", 
                       keytype="WORMBASE",
                       multiVals="first")


et_annot$GOID <-   mapIds(org.Ce.eg.db, 
                            keys=row.names(et_annot), 
                            column="GO", 
                            keytype="WORMBASE",
                            multiVals="first")

et_annot$term <- mapIds(GO.db, 
                        keys=et_annot$GOID, 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")
et_annot$term <- as.character(et_annot$term)

top$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(top), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")
top$Name <- mapIds(org.Ce.eg.db, 
                       keys=row.names(top), 
                       column="GENENAME", 
                       keytype="WORMBASE",
                       multiVals="first")

et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]

### BIOTYPE ANNOT FIX IT!
taxon = 'Caenorhabditis elegans'
taxon = tolower(taxon)
tmp = unlist(strsplit(x = taxon, split = ' '))
dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
use.official.gene.symbol <- TRUE
needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene")
if(use.official.gene.symbol == TRUE){
  gmt_flt = getBM(attributes=needed.attributes,filters="external_gene_name",values=et_annot$symbol, mart=mart)
  gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"external_gene_name"])),]
  rownames(gmt_flt) = gmt_flt[,"external_gene_name"]
} else {
  gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(et_annot), mart=mart)
  gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
  rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]
}

gmt_flt[,"description"] = gsub(pattern = "\\[Source:.*", replacement = "", x = gmt[,"description"], ignore.case = T,perl = FALSE)
gmt_flt_freq <- as.data.frame(table(unlist(gmt_flt$gene_biotype)))
pdf(file = "Filtered genes biotype distribution.pdf", width = 12, height = 17, family = "Helvetica")
pie(gmt_flt_freq$Freq, labels = gmt_flt_freq$Var1, radius = 0.5, main = "Filtered genes biotype distribution")
dev.off()

#TESTING A HYPOTESIS
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
  print("Correction protocol executer! All LogFC have been inverted!")
  beep()
  et_annot$logFC <- et_annot$logFC*(-1)
}

### Distribution, with moda and median

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

dir.create("distributions")
setwd("distributions")

png(file = "logCPM distribution, filtered.png")
hist(et_annot$logCPM, main = "logCPM distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM")
abline(v=median(et_annot$logCPM), col = "red")
abline(v=getmode(et_annot$logCPM), col = "blue")
dev.off()
png(file = "p-value distribution, filtered.png")
hist(et_annot$PValue, main = "p-value distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value")
abline(v=median(et_annot$PValue), col = "red")
abline(v=getmode(et_annot$PValue), col = "blue")
dev.off()
png(file = "LogFC distribution, filtered.png")
hist(et_annot$logFC, main = "logFC distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC")
abline(v=median(et_annot$logFC), col = "red")
abline(v=getmode(et_annot$logFC), col = "blue")
dev.off()

png(file = "logCPM distribution, nonfiltered.png")
hist(et_annot_non_filtered$logCPM, main = "logCPM distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM")
abline(v=median(et_annot_non_filtered$logCPM), col = "red")
abline(v=getmode(et_annot_non_filtered$logCPM), col = "blue")
dev.off()
png(file = "p-value distribution, nonfiltered.png")
hist(et_annot_non_filtered$PValue, main = "p-value distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value")
abline(v=median(et_annot_non_filtered$PValue), col = "red")
abline(v=getmode(et_annot_non_filtered$PValue), col = "blue")
dev.off()
png(file = "LogFC distribution, nonfiltered.png")
hist(et_annot_non_filtered$logFC, main = "logFC distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC")
abline(v=median(et_annot_non_filtered$logFC), col = "red")
abline(v=getmode(et_annot_non_filtered$logFC), col = "blue")
dev.off()

setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}

### Simple summary
all <- nrow(CountsTable)
allpadj <- sum(et_annot$PValue < pvalue_cutoff, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > logfchigh_cutoff, na.rm=TRUE)
down <- sum(et_annot$logFC < logfclow_cutoff, na.rm=TRUE)
header <- c('all genes', 'mean of logCPM', 'padj<0,05', 'genes with > high', 'genes with < low')
meaning <- c(print(all), print(avg_cpm), print(allpadj), print(up), print(down))
df <- data.frame(header, meaning)


# WRITE RESULTS
write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(top, file = "Results edgeR.xlsx", sheetName = "Top Tags (with FDR)", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)


### MY GO TESTS
setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}

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
                     keys=row.names(my_go), 
                     column="TERM", 
                     keytype="GOID",
                     multiVals="first")
my_go$term <- as.character(my_go$term)
names(my_go) <- c("term", "logFC", "genes in term", "name of term")

for_hm <- my_go
my_go <- my_go[order(my_go$`genes in term`, decreasing = TRUE),]
write.xlsx(my_go, file = "My GO.xlsx", sheetName = "GO classes", append = TRUE)


### GO HEATMAP
## FIX IT MUTAFUKA
cpmfgh <- as.data.frame(logCPM)

cpmfgh$GOID <-            mapIds(org.Ce.eg.db, 
                                 keys=row.names(cpmfgh), 
                                 column="GO", 
                                 keytype="WORMBASE",
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
pdf(file = "GO heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(go_for_heatmap, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=0.7, cexCol=0.7, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,15))

dev.off()

### ANNOTATE COUNTS

CountsTable$symbol <- mapIds(org.Ce.eg.db, 
                             keys=row.names(CountsTable), 
                             column="SYMBOL", 
                             keytype="WORMBASE",
                             multiVals="first")

CountsTable$name <- mapIds(org.Ce.eg.db, 
                           keys=row.names(CountsTable), 
                           column="GENENAME", 
                           keytype="WORMBASE",
                           multiVals="first")



### SPLITTING DATA
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

for_fisher <- as.data.frame(subset(et_annot_non_filtered, PValue < pvalue_cutoff))

for_fisher_high <- as.data.frame(subset(et_annot_non_filtered, logFC > 0))
for_fisher_low <- as.data.frame(subset(et_annot_non_filtered, logFC < 0))



### GOANA

goana_up <- goana(de = et_annot_high$entrez, species = "Ce")
go_up_30 <- topGO(goana_up, n=30)
go_up_100 <- topGO(goana_up, n=100)
go_up_500 <- topGO(goana_up, n=500)

goana_down <- goana(de = et_annot_low$entrez, species = "Ce")
go_down_30 <- topGO(goana_down, n=30)
go_down_100 <- topGO(goana_down, n=100)
go_down_500 <- topGO(goana_down, n=500)


write.xlsx(go_up_30, file = "Goana GO tests, upreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_up_100, file = "Goana GO tests, upreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_up_500, file = "Goana GO tests, upreg.xlsx", sheetName = "top500", append = TRUE)

write.xlsx(go_down_30, file = "Goana GO tests, downreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_down_100, file = "Goana GO tests, downreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_down_500, file = "Goana GO tests, downreg.xlsx", sheetName = "top500", append = TRUE)




### FISHER GO TESTS
  
GOFisherBP <- function(df, nodes, nrows, p){
    all_genes <- c(df$logFC)
    names(all_genes) <- rownames(df)
    go_data <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
                     p, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", 
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
                     p, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", 
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
                     p, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", 
                   ID = "ENSEMBL", nodeSize = nodes)
    go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
    go_table <- GenTable(go_data, weightFisher = go_test,
                         orderBy = "weightFisher", ranksOf = "weightFisher",
                         topNodes = nrows)
    return(go_table)
  }

gobp_h_50 <- GOFisherBP(for_fisher_high, 5, 50, 0.05)
gomf_h_50 <- GOFisherMF(for_fisher_high, 5, 50, 0.05)
gocc_h_50 <- GOFisherCC(for_fisher_high, 5, 50, 0.05)
  
gobp_l_50 <- GOFisherBP(for_fisher_low, 5, 50, 0.05)
gomf_l_50 <- GOFisherMF(for_fisher_low, 5, 50, 0.05)
gocc_l_50 <- GOFisherCC(for_fisher_low, 5, 50, 0.05)
  
gobp_h_100 <- GOFisherBP(for_fisher_high, 5, 100, 0.05)
gomf_h_100 <- GOFisherMF(for_fisher_high, 5, 100, 0.05)
gocc_h_100 <- GOFisherCC(for_fisher_high, 5, 100, 0.05)
  
gobp_l_100 <- GOFisherBP(for_fisher_low, 5, 100, 0.05)
gomf_l_100 <- GOFisherMF(for_fisher_low, 5, 100, 0.05)
gocc_l_100 <- GOFisherCC(for_fisher_low, 5, 100, 0.05)
  
gobp_h_250 <- GOFisherBP(for_fisher_high, 5, 250, 0.05)
gomf_h_250 <- GOFisherMF(for_fisher_high, 5, 250, 0.05)
gocc_h_250 <- GOFisherCC(for_fisher_high, 5, 250, 0.05)
  
gobp_l_250 <- GOFisherBP(for_fisher_low, 5, 250, 0.05)
gomf_l_250 <- GOFisherMF(for_fisher_low, 5, 250, 0.05)
gocc_l_250 <- GOFisherCC(for_fisher_low, 5, 250, 0.05)

write.xlsx(gobp_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 50", append = TRUE)
write.xlsx(gomf_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 50", append = TRUE)
write.xlsx(gocc_h_50, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 50", append = TRUE) 
write.xlsx(gobp_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(gomf_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(gocc_h_100, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 100", append = TRUE)
write.xlsx(gobp_h_250, file = "GO_Fisher_upreg.xlsx", sheetName = "BP, top 250", append = TRUE)
write.xlsx(gomf_h_250, file = "GO_Fisher_upreg.xlsx", sheetName = "MF, top 250", append = TRUE)
write.xlsx(gocc_h_250, file = "GO_Fisher_upreg.xlsx", sheetName = "CC, top 250", append = TRUE)
    
write.xlsx(gobp_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 50", append = TRUE)
write.xlsx(gomf_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 50", append = TRUE)
write.xlsx(gocc_l_50, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 50", append = TRUE)
write.xlsx(gobp_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(gomf_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(gocc_l_100, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 100", append = TRUE)
write.xlsx(gobp_l_250, file = "GO_Fisher_downreg.xlsx", sheetName = "BP, top 250", append = TRUE)
write.xlsx(gomf_l_250, file = "GO_Fisher_downreg.xlsx", sheetName = "MF, top 250", append = TRUE)  
write.xlsx(gocc_l_250, file = "GO_Fisher_downreg.xlsx", sheetName = "CC, top 250", append = TRUE)




## CORELLATION MARIX
correl <- cpm(y)
x <- cor(correl)
pdf(file = "Corellation matrix.pdf", width = 10, height = 10)
heatmap.2(x, col=col.pan,Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,10), lhei=c(2,10), lwid=c(2,6), main = "Corellation Matrix")
dev.off()
#heatmap.2(x, margins = c(10,10))      
### MDS PLOT
pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
# pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
pdf(file = "MDSPlot.pdf", width = 10, height = 10)
plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition])
legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
dev.off()

### VOLCANO PLOT
allgenes <- nrow(et_annot_non_filtered)
et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > logfchigh_cutoff & et_annot_non_filtered$PValue < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 10, height = 10)
g = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=1, size=2) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw()
g
dev.off()



### REACTOME PART

dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "celegans", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "All reactome", append = TRUE)
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9, legend.text = "Top 30 Reactome pathways")
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges, categorySize="pvalue", showCategory = 10)
dev.off()

#HIGH
df_high <- et_annot_high$entrez
x <- enrichPathway(gene=df_high, organism = "celegans", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "High", append = TRUE)
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "celegans", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Low", append = TRUE)
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

###KEGG expression profile (without lfc, but with generatio)
df_high <- et_annot_high$uniprot
kk_up <- enrichKEGG(gene = df_high, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk_up, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk_up, showCategory=30,  font.size = 9)
dev.off()


df_low <- et_annot_low$uniprot
kk_down <- enrichKEGG(gene = df_low, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk_down, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk_down, showCategory=30,  font.size = 9)
dev.off()


### KEGGA
keg_com <- kegga(de = et_annot$entrez, species="Ce")
tk_common <- topKEGG(keg_com, n=100)
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Ce")
tk_up <- topKEGG(keg_up, n=30)
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Ce")
tk_down <- topKEGG(keg_down, n=30)
write.xlsx(tk_down, file = "kegga.xlsx", sheetName = "Downreg", append = TRUE)

rownames(tk_common) <- substring(rownames(tk_common), 6)


# TOP 100 PVALUE GENES
if (heatmaps == TRUE){
  logCPM <- as.data.frame(logCPM)
  rownames(logCPM) <- make.names(et$genes$Symbol, unique = TRUE)
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
  # colnames(topcpm) <- paste(y$samples$group, 1:2, sep="-")
  topcpm$Symbol <- mapIds(org.Ce.eg.db, 
                          keys=row.names(topcpm), 
                          column="SYMBOL", 
                          keytype="WORMBASE",
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

dir.create("mdplots")
setwd("mdplots")
for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}

setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


### SEARCH AND PLOT!
if (custom_heatmap == TRUE) {
  logCPM$Name <- y$genes$Name
  let <- c("CD")
  
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
dir.create("kegg")
setwd("kegg")
plot_pathway = function(pid){
  pathview(gene.data=foldchanges, 
           pathway.id=pid, 
           species="cel", 
           new.signature=FALSE)
}

for (f in kk_up$ID){plot_pathway(paste(f))}
for (f in kk_down$ID){plot_pathway(paste(f))}


### PANTHER.DB
setwd(directory)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

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


###DESEQ2
directory <- '~/counts/worm_test/'

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res <- results(dds, tidy = FALSE )
rld<- rlogTransformation(dds, blind=TRUE)

directory <- '~/counts/worm_test/results'
setwd(directory)

if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


pdf(file = "PCAPlot.pdf", width = 10, height = 10)
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

if (logging == TRUE){
  sink(type="message")
  close(zz)
}


res_df <- as.data.frame(res)
res_nf <- as.data.frame(res)
res_df <- res_df[complete.cases(res_df), ]
res_df <- as.data.frame(subset(res_df, res_df$baseMean > baseMean_cutoff))
res_df <- as.data.frame(subset(res_df, res_df$padj < pvalue_cutoff))
res_df <- as.data.frame(subset(res_df, res_df$log2FoldChange > logfchigh_cutoff | res_df$log2FoldChange < logfclow_cutoff))

res_df_high <- as.data.frame(subset(res_df, res_df$log2FoldChange > 0))
res_df_low <- as.data.frame(subset(res_df, res_df$log2FoldChange < 0))

com_deseq_edger_high <- as.data.frame(intersect(rownames(res_df_high), rownames(et_annot_high)))
com_deseq_edger_low <- as.data.frame(intersect(rownames(res_df_low), rownames(et_annot_low)))

proc_high<- (nrow(com_deseq_edger_high)/nrow(et_annot_high))*100
proc_low <- (nrow(com_deseq_edger_low)/nrow(et_annot_low))*100

sum <- NULL
sum <- data.frame(paste(proc_high), paste(proc_low))
names(sum) <- c("Upreg common genes", "Downreg common genes")

sign_proc_edger <- (nrow(et_annot)/nrow(et_annot_non_filtered))*100
sign_proc_deseq <- (nrow(res_df)/nrow(res_nf))*100

aproc <- data.frame(paste(sign_proc_edger), paste(sign_proc_deseq))

names(aproc) <- c("edgeR, significance %", "DESeq2, significance %")

write.xlsx(sum, file = "DESeq and edgeR comparsion.xlsx", sheetName = "Common genes of upreg and downreg", append = TRUE)
write.xlsx(aproc, file = "DESeq and edgeR comparsion.xlsx", sheetName = "Percent of significance, w.filtering", append = TRUE)

### COMPARE EDGER AND DESEQ2
com <- as.data.frame(intersect(rownames(res_nf), rownames(et_annot_non_filtered)))
names(com) <- c("sas")
edger_com <- data.frame()
deseq_com <- data.frame()

edger_com <- et_annot_non_filtered[(rownames(et_annot_non_filtered))%in% com$sas ,]
deseq_com <- res_nf[(rownames(res_nf))%in% com$sas ,]

common <- data.frame()
common <- data.frame(edger_com$logFC, deseq_com$log2FoldChange)
rownames(common) <- rownames(edger_com)
common <- transform(common, SD=apply(common,1, sd, na.rm = TRUE))

common <- common[order(common$SD, decreasing = TRUE),]

png(file = "Deviance between edgeR and deseq2.png")
plot(common$SD, type= "l", main= "Deviance between edgeR and deseq2(should be approximately zero)")
dev.off()

deviant <- common[seq(1:10),]


deviant$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(deviant), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
deviant$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(deviant), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

write.xlsx(deviant, file = "Top 10 deviant genes between deseq2 and edgeR.xlsx")