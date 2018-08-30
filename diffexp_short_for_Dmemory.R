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
qlm_test <- FALSE
logging <- FALSE
motiv <- TRUE
boxplots <- TRUE
biotype <- FALSE
distrib <- TRUE
### CONSTANTS BLOCK

pvalue_cutoff <- 0.05
fdr_cutoff <- 0.05
logfchigh_cutoff <- 0
logfclow_cutoff <- 0
cpm_cutoff <- -100
gs_size <- 10
diseases_set <- 50
go_terms_set <- 50
genes_in_term <- 3
filter_thresh <- 5
baseMean_cutoff <- 1.5
significant_enriched_motif <- 5
go_heatmap_count <- 20
stattest_number <- 1


directory <- '~/counts/Dmel.memory.29.06.18/HSP+/'
setwd(directory)
gr_control <- c("F")
gr_case <- c("mem")

### BUILDING A SPECIFIC DESIGN TABLE
if (logging == TRUE){
  zz <- file("error.log", open="wt")
  sink(zz, type="message")
}

col.pan <- colorpanel(100, "blue", "white", "red")
###DIRECTORY WHERE SAMPLES ARE LOCATED
library(dplyr)
if (analyze_all_samples == TRUE){
  sampleFiles <- grep('counts',list.files(directory),value=TRUE)
  sampleFiles
  sampleCondition <- c("F24-new", "F24-new", "F24-new",
                       "F24-old", "F24-old",
                       "F5-new", "F5-new", "F5-new", 
                       "F5-old", "F5-old", 
                       "K-new", "K-new", "K-new", 
                       "K-old", "K-old")
  
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

et_annot$GOID <-     mapIds(org.Dm.eg.db, 
                            keys=row.names(et_annot), 
                            column="GO", 
                            keytype="FLYBASE",
                            multiVals="first")

et_annot$term <- mapIds(GO.db, 
                        keys=et_annot$GOID, 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")

et_annot$term <- as.character(et_annot$term)

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


all <- nrow(y$counts)
allpadj <- sum(et_annot$PValue < pvalue_cutoff, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > logfchigh_cutoff, na.rm=TRUE)
down <- sum(et_annot$logFC < logfclow_cutoff, na.rm=TRUE)
header <- c('all genes', 'mean of logCPM', 'padj<0,05', 'genes with > high', 'genes with < low')
meaning <- c(print(all), print(avg_cpm), print(allpadj), print(up), print(down))
df <- data.frame(header, meaning)

et.full <- nrow(et_annot)
et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < fdr_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
#et_annot <- et_annot[complete.cases(et_annot), ]



# WRITE RESULTS
write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(top, file = "Results edgeR.xlsx", sheetName = "Top Tags (with FDR)", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)







####INTERSECTING#####
wt_kf <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP+/results/K-F/Results edgeR.xlsx", sheetIndex = 3)
wt_fmem <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP+/results/F-mem/Results edgeR.xlsx", sheetIndex = 3)
wt_kmem <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP+/results/K-mem/Results edgeR.xlsx", sheetIndex = 3)

hsp_kf <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP-/results/cont-F5/Results edgeR.xlsx", sheetIndex = 3)
hsp_fmem <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP-/results/F5-F24/Results edgeR.xlsx", sheetIndex = 3)
hsp_kmem <- read.xlsx("~/counts/Dmel.memory.29.06.18/HSP-/results/cont-F24/Results edgeR.xlsx", sheetIndex = 3)


intersect(wt_fmem$NA., hsp_fmem$NA.)
intersect(wt_kmem$NA., hsp_kmem$NA.)
intersect(wt_kf$NA., hsp_kf$NA.)






