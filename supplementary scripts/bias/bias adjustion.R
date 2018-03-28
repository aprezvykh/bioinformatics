library(GOplot)
library(r2excel)
library(BatchJobs)
library(BiocParallel)
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
library(RColorBrewer)
library(Matrix)
library(PANTHER.db)
library(GO.db)
library(Hmisc)
library(DESeq2)
library(checkmate)
library(biomaRt)
library(RcisTarget)
library(devtools)
library(enrichR)
library(zoo)
library(rvest)
library(XML)
library(plyr)
library(AnnotationDbi)
library(ggsignif)
library(grid)
library(pcaExplorer)
library(ggbiplot)
library(cowplot)
library(gapminder)
library(gridExtra)

###BIOCBARALLEL SETTINGS
default <- registered()
register(BatchJobsParam(workers = 10), default = TRUE)
options(MulticoreParam=quote(MulticoreParam(workers=8)))
param <- SnowParam(workers = 2, type = "SOCK")


heatmaps <- TRUE
custom_genes_plots <- FALSE
analyze_all_samples <- FALSE
disease_association <- FALSE
kegg_plots <- TRUE
panther_analysis <- TRUE
deseq2_part <- FALSE
qlm_test <- TRUE
logging <- FALSE
motiv <- TRUE
boxplots <- TRUE
biotype <- FALSE
distrib <- TRUE
### CONSTANTS BLOCK

pvalue_cutoff <- 0.05
fdr_cutoff <- 0.5
logfchigh_cutoff <- 0
logfclow_cutoff <- -0
cpm_cutoff <- -1
gs_size <- 10
diseases_set <- 50
go_terms_set <- 50
genes_in_term <- 3
filter_thresh <- 5
baseMean_cutoff <- 1.5
significant_enriched_motif <- 5
go_heatmap_count <- 20
stattest_number <- 1

directory <- '~/counts/BAP/'
setwd(directory)
gr_control <- c("new_F")
gr_case <- c("new_M")

if (analyze_all_samples == TRUE){
  sampleFiles <- grep('fly',list.files(directory),value=TRUE)
  sampleFiles
  sampleCondition <- c("old_k", "old_k", "old_n", "old_n", 
                       "old_f", "old_f", "old_m", "old_m", 
                       "new_f", "new_f", "new_k", "new_k", 
                       "new_m", "new_m")
  
  sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
  col <- as.vector(sampleTable$sampleName)
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



setwd("results")
stattest <- paste(gr_control, gr_case, sep = "-")
results.dir <- paste(directory, "results", sep = "")
setwd(results.dir)


if (analyze_all_samples == FALSE){
  dir.create(stattest)
  setwd(stattest)
} else if (analyze_all_samples == TRUE){
  dir.create("all")
  setwd("all")
}
# PLOTTING HTSEQ QUALITY BARPLOTS
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )


### DIFFEXPRESSION STATISTICAL ANALYSIS - EXACT NEGATIVE-
### BINOMIAL OR QLM TEST

if (qlm_test == TRUE){ 
  a <- DGEList(counts=y, group = sampleTable$condition) 
  CountsTable <- as.data.frame(y$counts)
  cpm <- cpm(y) 
  cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
  cpm <- as.data.frame(cpm(y)) 
  cpm$rowsum <- rowSums(cpm) 
  keep <- rowSums(cpm > cpm_cutoff) >= ncol(sampleTable) 
  logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
  logCPM <- logCPM[!(row.names(logCPM) %in% row.names.remove), ]
  logCPM <- logCPM[keep,]
  a <- a[keep, , keep.lib.sizes=FALSE] 
  a <- calcNormFactors(a, method = "TMM", doWeighting = TRUE) 
  design <- model.matrix(~sampleTable$condition)
  a <- estimateDisp(a,design, tagwise = TRUE, trend.method = "loess") 
  fit <- glmQLFit(a,design, robust = TRUE, abundance.trend = TRUE) 
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  et_annot_non_filtered <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  top <- as.data.frame(topTags(qlf, n = 20))
  et <- exactTest(a)
  hist(a$common.dispersion)
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
  #et_annot <- as.data.frame(et$table)
  et_annot <- as.data.frame(topTags(et, n = nrow(logCPM), adjust.method = "BH", sort.by = "PValue"))
  et_annot_non_filtered <- as.data.frame(et$table)
  
  
}

counts_control <- CountsTable[,grep(gr_control, names(CountsTable))]
counts_case <- as.data.frame(CountsTable[,grep(gr_case, names(CountsTable))])
counts_control$rowsum.control <- rowSums(counts_control)
counts_case$rowsum.case <- rowSums(counts_case)
diff <- data.frame(counts_control$rowsum.control, counts_case$rowsum.case)
rownames(diff) <- rownames(CountsTable) 

fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
lfgene
c <- grep(paste(lfgene), rownames(diff))
dfc <- as.data.frame(diff[c,])
lfgenefc <- fc[4,1]
stat <- dfc$counts_control.rowsum.control > dfc$counts_case.rowsum.case
lfgene
if (stat == TRUE & lfgenefc < 0){
  print("No Correction Needed!")
} else {
  et_annot$logFC <- et_annot$logFC*(-1)
  print("Correction protocol executed, logFC have been inverted!")
}



taxon = 'Drosophila Melanogaster'
taxon = tolower(taxon)
tmp = unlist(strsplit(x = taxon, split = ' '))
dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene", "transcript_length")

gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(et_annot), mart=mart)
gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]

###transcript length-logFC
et_annot.trlen <- et_annot
trlen.common <- intersect(rownames(et_annot.trlen), gmt_flt$ensembl_gene_id)
et_annot.trlen.common <- et_annot.trlen[(rownames(et_annot.trlen) %in% trlen.common),]
gmt_flt.common <- gmt_flt[(gmt_flt$ensembl_gene_id %in% trlen.common),]


et_annot.trlen.common <- et_annot.trlen.common[order(rownames(et_annot.trlen.common), decreasing = TRUE),]
gmt_flt.common <- gmt_flt.common[order(gmt_flt.common$ensembl_gene_id, decreasing = TRUE),]


trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM, et_annot.trlen.common$PValue, et_annot.trlen.common$FDR)
names(trlen.distrib) <- c("fc", "tl", "cpm", "pvalue", "FDR")

##trlen - logFC
trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM, et_annot.trlen.common$PValue, et_annot.trlen.common$FDR)
names(trlen.distrib) <- c("fc", "tl", "cpm", "pvalue", "FDR")
trlen.distrib <- trlen.distrib[which(trlen.distrib$cpm > 5),]
trlen.distrib$tl <- log2(trlen.distrib$tl)
fit <- lm(fc~tl, data = trlen.distrib)
fit$coefficients
rmse <- round(sqrt(mean(resid(fit)^2)), 2)
coefs <- coef(fit)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(fit)$r.squared, 2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                r^2 == .(r2) * "," ~~ RMSE == .(rmse))
plot(fc ~ tl, data = trlen.distrib, pch = ".", xlab = "Log2(Transcript Length)", ylab = "Log2(Fold Change)", main = "Transcript Length vs LogFC (all genes, logCPM>2)")
abline(fit)
cor(trlen.distrib$fc, trlen.distrib$tl)

##trlen - CPM
trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM, et_annot.trlen.common$PValue, et_annot.trlen.common$FDR)
names(trlen.distrib) <- c("fc", "tl", "cpm", "pvalue", "FDR")
tr.len.median <- median(trlen.distrib$tl)
 

m.low <- trlen.distrib[which(trlen.distrib$tl < tr.len.median),]$cpm 
m.low <- mean(10^m.low)
m.high <- trlen.distrib[which(trlen.distrib$tl > tr.len.median),]$cpm
m.high <- mean(10^m.high)
assym.coef <- log10(m.low/m.high)
assym.coef



max(g$transcript_length)
### INDEPENDENT NORMALISATION
a <- DGEList(counts=y, group = sampleTable$condition) 
df <- as.data.frame(a$counts)
gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(df), mart=mart)
gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]

gmt_flt.ordered <- gmt_flt[order(gmt_flt$transcript_length),]

s <- seq(from = 1, to = nrow(gmt_flt.ordered), by = round(nrow(gmt_flt.ordered)/10))


g <- gmt_flt[which(gmt_flt.ordered$transcript_length %in% s[2]:s[3]),]
g$transcript_length

norm.factors <- calcNormFactors(df[which(rownames(df) %in% g$ensembl_gene_id),], method = "TMM", doWeighting = TRUE) 
a$samples$norm.factors <- norm.factors
a$counts <- as.matrix(df[which(rownames(df) %in% g$ensembl_gene_id),])
design <- model.matrix(~sampleTable$condition)
a <- estimateDisp(a,design, tagwise = TRUE, trend.method = "loess") 
fit <- glmQLFit(a,design, robust = TRUE, abundance.trend = TRUE) 
qlf <- glmQLFTest(fit,coef=ncol(fit$design))
et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))


