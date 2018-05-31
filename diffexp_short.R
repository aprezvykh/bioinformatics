library(RDAVIDWebService)
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

### TFES
library(RcisTarget.mm9.motifDatabases.20k)
#data("mm9_10kbpAroundTss_motifRanking")
#data("mm9_direct_motifAnnotation")

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
deseq2_part <- TRUE
qlm_test <- FALSE
logging <- FALSE
motiv <- TRUE
boxplots <- TRUE
biotype <- FALSE
### CONSTANTS BLOCK

pvalue_cutoff <- 0.05
logfchigh_cutoff <- 0
logfclow_cutoff <- 0
cpm_cutoff <- 0.5
gs_size <- 10
diseases_set <- 50
number_of_kegg_plots <- 100
go_terms_set <- 50
pathways_set <- 30
genes_in_term <- 3
filter_thresh <- 5
baseMean_cutoff <- 1.5
significant_enriched_motif <- 5
go_heatmap_count <- 20
stattest_number <- 1


#a <- read.xlsx("tests.xlsx", sheetIndex = 1)
#a <- data.frame(a$control, a$case)
#gr_control <- as.character(a[1,1])
#gr_case <- as.character(a[1,2])

directory <- '~/counts/Early.markers/search-early/'
setwd(directory)
gr_control <- c("Control")
gr_case <- c("Tg")

### BUILDING A SPECIFIC DESIGN TABLE
if (logging == TRUE){
  zz <- file("error.log", open="wt")
  sink(zz, type="message")
}

col.pan <- colorpanel(100, "blue", "white", "red")
###DIRECTORY WHERE SAMPLES ARE LOCATED
library(dplyr)

if (analyze_all_samples == TRUE){
  sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
  sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                       'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                       'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                       'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                       'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')
  
  
  sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
  col <- as.vector(sampleTable$sampleName)
  names(col) <- c("yellow", "yellow", "yellow", "yellow", "yellow", 
                  "purple", "purple", "purple", "purple", "purple", 
                  "green", "green", "green", "green", "green", 
                  "red", "red", "red", "red", 
                  "blue", "blue", "blue", "blue", "blue")
  
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
  a <- calcNormFactors(a, method = "TMM") 
  design <- model.matrix(~sampleTable$condition) 
  a <- estimateDisp(a,design) 
  fit <- glmQLFit(a,design = design, robust = TRUE) 
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  et_annot_non_filtered <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  top <- as.data.frame(topTags(qlf, n = 20))
  et <- exactTest(a)
} else if (qlm_test == FALSE){ 
  design <- model.matrix(~sampleTable$condition) 
  normalized_lib_sizes <- calcNormFactors(y, method = "TMM") 
  CountsTable <- as.data.frame(y$counts) 
  raw_counts <- as.data.frame(y$counts) 
  y <- calcNormFactors(y, method = "TMM") 
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




et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]

et_annot <- as.data.frame(subset(et_annot, logFC > 1 | logFC < -1))
et_annot <- as.data.frame(subset(et_annot, FDR < 0.01))
et_annot <- as.data.frame(subset(et_annot, PValue < 0.01))


m <- et_annot
write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(top, file = "Results edgeR.xlsx", sheetName = "Top Tags (with FDR)", append = TRUE)
write.xlsx(m, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)

goana_straight <- goana(de = et_annot$entrez, species = "Mm")
go_all_1000 <- topGO(goana_straight, n=1000)
