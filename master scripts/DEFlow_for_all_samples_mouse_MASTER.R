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
### TFES
library(RcisTarget.mm9.motifDatabases.20k)
data("mm9_10kbpAroundTss_motifRanking")
data("mm9_direct_motifAnnotation")

heatmaps <- TRUE
custom_genes_plots <- FALSE
fisherGO <- TRUE
analyze_all_samples <- FALSE
disease_association <- TRUE
kegg_plots <- TRUE
panther_analysis <- TRUE
deseq2_part <- TRUE
qlm_test <- TRUE
logging <- FALSE
motiv <- FALSE
boxplots <- TRUE
### CONSTANTS BLOCK

pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
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

directory <- '~/counts/ALS Mice/experimental/'
setwd(directory)
gr_control <- c("Tg-1")
gr_case <- c("Tg-3")

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
  fit <- glmQLFit(a,design, robust = TRUE) 
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  et_annot_non_filtered <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
  top <- as.data.frame(topTags(qlf, n = 20))
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
nrow(logCPM)
write.csv(logCPM, file = "overall logCPM.csv")
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


top$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(top), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
top$Name<- mapIds(org.Mm.eg.db, 
                       keys=row.names(top), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")


top <- top[complete.cases(top),]
### CPM
cpm <- as.data.frame(cpm(y))

cpm$Symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$Name <- mapIds(org.Mm.eg.db, 
                   keys=row.names(cpm), 
                   column="GENENAME", 
                   keytype="ENSEMBL",
                   multiVals="first")

cpm$entrez <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$GOID <-     mapIds(org.Mm.eg.db, 
                       keys=row.names(cpm), 
                       column="GO", 
                       keytype="ENSEMBL",
                       multiVals="first")

cpm$term <- mapIds(GO.db, 
                   keys=cpm$GOID, 
                   column="TERM", 
                   keytype="GOID",
                   multiVals="first")
cpm$term <- as.character(cpm$term)
###TopTags Boxplots

if (boxplots == TRUE){
  dir.create("TopTags Boxplots")
  setwd("TopTags Boxplots")
    for (f in rownames(top)){
      r <- grep(paste(f), rownames(cpm), ignore.case = TRUE)
      thm <- cpm[r,]
      rownames(thm) <- thm$Symbol
      plot.name <- as.character(rownames(thm))
      plot.description <- as.character(thm$Name)
      plot.term <- as.character(thm$term)
      thm$Symbol <- NULL
      thm$Name <- NULL
      thm$GOID <- NULL
      thm$term <- NULL
      thm$entrez <- NULL
      colnames(thm) <- sampleCondition
      thm <-as.data.frame(t(thm))
      thm$Condition <- rownames(thm)
      thm$Group <- thm$Condition
      names(thm) <- c("gene", "Condition", "Group")
      #pdf(file = paste(plot.name, "pdf", sep = "."), width = 10, height = 10, family = "Helvetica")
      #png(paste(plot.name, "-", src, ".png", sep=""), height = 600, width = 800)
      g <- ggplot(thm, aes(x = Condition, y = gene)) + 
        geom_boxplot(data = thm, aes(fill = Group), alpha = 0.5) + 
        scale_x_discrete(name = "Experimental Groups") + 
        scale_y_continuous(name = "Counts per million") + 
        theme_bw() + 
        geom_signif(comparisons = list(c(paste(gr_control)), paste(gr_case)), map_signif_level = TRUE) 
        
      g <- g + ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
        theme(plot.title = element_text(hjust = 0.5))
                    
        ggsave(filename = paste(plot.name, "pdf", sep = "."), plot = g)

    }
}
setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}
getwd()
### FILTRATION
et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]

### BIOTYPE ANNOT FIX IT!
taxon = 'Mus Musculus'
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

# et_annot$biotype <- gmt_flt$gene_biotype
gmt_flt_freq <- as.data.frame(table(unlist(gmt_flt$gene_biotype)))

pdf(file = "Filtered genes biotype distribution.pdf", width = 12, height = 17, family = "Helvetica")
bp <- ggplot(gmt_flt_freq, aes(x="", y=Freq, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0) + theme_bw() +theme(axis.title.x=element_blank(),
                                                             axis.text.x=element_blank(),
                                                             axis.ticks.x=element_blank(),
                                                             axis.title.y=element_blank(),
                                                             axis.text.y=element_blank(),
                                                             axis.ticks.y=element_blank())

pie
dev.off()


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
  print("Correction protocol executed, logFC have been inverted!")
}



### Distribution, with moda and median

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

dir.create("distributions")
setwd("distributions")

png(file = "logCPM distribution, nonfiltered.png")
hist(et_annot_non_filtered$logCPM, main = "logCPM distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM", breaks = 50)
abline(v=median(et_annot_non_filtered$logCPM), col = "red")
abline(v=getmode(et_annot_non_filtered$logCPM), col = "blue")

dev.off()
png(file = "logCPM distribution, filtered.png")
hist(et_annot$logCPM, main = "logCPM distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "logCPM", breaks = 50)
abline(v=median(et_annot$logCPM), col = "red")
abline(v=getmode(et_annot$logCPM), col = "blue")

dev.off()
png(file = "p-value distribution, nonfiltered.png")
hist(et_annot_non_filtered$PValue, main = "p-value distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value", breaks = 50)
abline(v=median(et_annot_non_filtered$PValue), col = "red")
abline(v=getmode(et_annot_non_filtered$PValue), col = "blue")

dev.off()
png(file = "p-value distribution, filtered.png")
hist(et_annot$PValue, main = "p-value distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "p-value", breaks = 50)
abline(v=median(et_annot$PValue), col = "red")
abline(v=getmode(et_annot$PValue), col = "blue")

dev.off()
png(file = "LogFC distribution, nonfiltered.png")
hist(et_annot_non_filtered$logFC, main = "logFC distribution, non-filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC", breaks = 50)
abline(v=median(et_annot_non_filtered$logFC), col = "red")
abline(v=getmode(et_annot_non_filtered$logFC), col = "blue")

dev.off()
png(file = "LogFC distribution, filtered.png")
hist(et_annot$logFC, main = "logFC distribution, filtered", freq = TRUE, col = col.pan, labels = TRUE, xlab = "LogFC", breaks = 50)
abline(v=median(et_annot$logFC), col = "red")
abline(v=getmode(et_annot$logFC), col = "blue")
dev.off()



setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}

### Simple summary
all <- nrow(y$counts)
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
pdf(file = "GO heatmap.pdf", width = 12, height = 17, family = "Helvetica")
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


### GOANA
et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))
goana_up <- goana(de = et_annot_high$entrez, species = "Mm")
go_up_30 <- topGO(goana_up, n=30)
go_up_30$perc = (go_up_30$DE/go_up_30$N)*100
go_up_100 <- topGO(goana_up, n=100)
go_up_100$perc = (go_up_100$DE/go_up_100$N)*100
go_up_500 <- topGO(goana_up, n=500)
go_up_500$perc = (go_up_500$DE/go_up_500$N)*100
go_up_500$perc <- round(go_up_500$perc, digits = 4)

goana_down <- goana(de = et_annot_low$entrez, species = "Mm")
go_down_30 <- topGO(goana_down, n=30)
go_down_30$perc = (go_down_30$DE/go_down_30$N)*100
go_down_100 <- topGO(goana_down, n=100)
go_down_100$perc = (go_down_100$DE/go_down_100$N)*100
go_down_500 <- topGO(goana_down, n=500)
go_down_500$perc = (go_down_500$DE/go_down_500$N)*100
go_down_500$perc <- round(go_down_500$perc, digits = 4)




write.xlsx(go_up_30, file = "Goana GO tests, upreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_up_100, file = "Goana GO tests, upreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_up_500, file = "Goana GO tests, upreg.xlsx", sheetName = "top500", append = TRUE)

write.xlsx(go_down_30, file = "Goana GO tests, downreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_down_100, file = "Goana GO tests, downreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_down_500, file = "Goana GO tests, downreg.xlsx", sheetName = "top500", append = TRUE)


write.xlsx(go_up_100, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated GO terms", append = TRUE)
write.xlsx(go_down_100, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated GO terms", append = TRUE)

###ANNOTATE LOGCPM FOR GO HEATMAPS
for_go_heatmap <- logCPM

for_go_heatmap$entrez <- mapIds(org.Mm.eg.db, 
                                keys=row.names(for_go_heatmap), 
                                column="ENTREZID", 
                                keytype="ENSEMBL",
                                multiVals="first")

for_go_heatmap$symbol <- mapIds(org.Mm.eg.db, 
                                keys=row.names(for_go_heatmap), 
                                column="SYMBOL", 
                                keytype="ENSEMBL",
                                multiVals="first")



### GO HEATMAPS OF VARIOUS GENES

all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez
GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", nodeSize = 2)

allGO <- genesInTerm(GOdata)
##UP
remove <- c("MF", "CC")
dir.create("GO heatmap plots upregulated")
setwd("GO heatmap plots upregulated")

go_up_subset <- subset(go_up_500, go_up_500$DE > 40)
go_up_subset <- subset(go_up_subset, go_up_subset$DE < 60)
go_up_subset <- subset(go_up_subset, go_up_subset$N < 400)
go_up_subset <- go_up_subset[!(go_up_subset$Ont %in% remove), ] 
up_terms <- data.frame(go_up_subset$Term)
rownames(up_terms) <- rownames(go_up_subset)


for (f in seq(1:nrow(go_up_subset))){
      z <- allGO[paste(rownames(up_terms))[f]]
      z <- as.data.frame(z)
      goid <- names(z)
      goname <- go_up_subset[f,1]
      goperc <- go_up_subset[f,6]
      gopvalue <- go_up_subset[f,5]
      statperc <- paste("Percent of involved genes", goperc, sep = ":")
      name <-paste(goid, goname, sep = ":")
      name <- as.character(name)
      name <- gsub('[.]', "", name)
      name <- paste(name, statperc, sep = "\n")
      names(z) <- c("sus")
      mat <- for_go_heatmap[(for_go_heatmap$entrez %in% z$sus),]
      rownames(mat) <- mat$symbol
      mat$symbol <- NULL
      mat$entrez <- NULL
      mat <- t(scale(t(mat)))
      pdf(file = paste(goid,"upreg_GO.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
      heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none",
                trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
                margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = paste(name))
      dev.off()
      
}

##DOWN

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

dir.create("GO heatmap plots downregulated")
setwd("GO heatmap plots downregulated")
go_down_subset <- subset(go_down_500, go_down_500$DE > 40)
go_down_subset <- subset(go_down_subset, go_down_subset$DE < 60)
go_down_subset <- subset(go_down_subset, go_down_subset$N < 400)
go_down_subset <- go_down_subset[!(go_down_subset$Ont %in% remove), ] 
down_terms <- data.frame(go_down_subset$Term)
rownames(down_terms) <- rownames(go_down_subset)

for (f in seq(1:nrow(go_down_subset))){
  z <- allGO[paste(rownames(down_terms))[f]]
  z <- as.data.frame(z)
  goid <- names(z)
  goname <- go_down_subset[f,1]
  goperc <- go_down_subset[f,6]
  gopvalue <- go_down_subset[f,5]
  statperc <- paste("Percent of involved genes", goperc, sep = ":")
  name <-paste(goid, goname, sep = ":")
  name <- as.character(name)
  name <- gsub('[.]', "", name)
  name <- paste(name, statperc, sep = "\n")
  names(z) <- c("sus")
  mat <- for_go_heatmap[(for_go_heatmap$entrez %in% z$sus),]
  rownames(mat) <- mat$symbol
  mat$symbol <- NULL
  mat$entrez <- NULL
  mat <- t(scale(t(mat)))
  pdf(file = paste(goid,"downreg_GO.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
  heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none",
            trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = paste(name))
  dev.off()
  
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

getwd()
### FISHER GO TESTS
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

for_fisher <- as.data.frame(subset(et_annot_non_filtered, PValue < pvalue_cutoff))

for_fisher_high <- as.data.frame(subset(for_fisher, logFC > 0))
for_fisher_low <- as.data.frame(subset(for_fisher, logFC < 0))

###colors

## CORELLATION MARIX
correl <- logCPM
correl <- as.data.frame(correl)
correl <- correl[complete.cases(correl),]
x <- cor(correl, method = "pearson")
pdf(file = "Corellation matrix Spearman.pdf", width = 10, height = 10)
heatmap.2(x, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")
dev.off()
#legend("bottomright",
#       legend = c("Control-1", "Control-3", "Tg-1", "Tg-2", "Tg-3"),
#       col = c("yellow", "purple", "green", "red", "blue"),
#       lty= 1,
#       lwd = 10)

#heatmap.2(x, margins = c(10,10))      
### MDS PLOT
getwd()

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

#SMEAR PLOT
tt <- topTags(qlf, n=nrow(y), p.value=0.1)
pdf(file = "Smear plot.pdf", width = 10, height = 10)
plotSmear(qlf, de.tags = rownames(tt$table))
dev.off()


### REACTOME PART
dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "mouse", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "All reactome", append = TRUE)
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
write.xlsx(x, "Reactome.xlsx", sheetName = "High", append = TRUE)
write.xlsx(x, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated Reactome pathways", append = TRUE)
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "mouse", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Low", append = TRUE)
write.xlsx(x, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated Reactome pathways", append = TRUE)
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

###KEGGA
keg_com <- kegga(de = et_annot$entrez, species="Mm")
tk_common <- topKEGG(keg_com, n=100)
tk_common <- subset(tk_common, tk_common$P.DE < 0.05)
tk_common$perc <- (tk_common$DE/tk_common$N)*100
tk_common <- tk_common[order(tk_common$perc, decreasing = TRUE),]
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Mm")
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Mm")
tk_down <- topKEGG(keg_down, n=100)
tk_down <- subset(tk_down, tk_down$P.DE < 0.05)
tk_down$perc <- (tk_down$DE/tk_down$N)*100
tk_down <- tk_down[order(tk_down$perc, decreasing = TRUE),]
write.xlsx(tk_down, file = "kegga.xlsx", sheetName = "Downreg", append = TRUE)
getwd()
rownames(tk_common) <- substring(rownames(tk_common), 6)


write.xlsx(tk_up, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated KEGG pathways", append = TRUE)
write.xlsx(tk_down, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated KEGG pathways", append = TRUE)


# TOP 100 PVALUE GENES
logCPM <- as.data.frame(logCPM)
nrow(logCPM)
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
logCPMfc <- logCPMfc[complete.cases(logCPMfc),]
logCPMfc <- t(scale(t(logCPMfc)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 logFC Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
dev.off()


### High Expressed Genes
topcpm <- as.data.frame(cpm(y))
topcpm$rowsum <- rowSums(topcpm)
topcpm <- topcpm[order(topcpm$rowsum, decreasing = TRUE),]
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[!(row.names(topcpm) %in% row.names.remove), ] 
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

### heatmap with all samples
topcpm <- as.data.frame(cpm(y))
topcpm$rowsum <- rowSums(topcpm)
topcpm <- topcpm[!(row.names(topcpm) %in% row.names.remove), ] 
topcpm <- topcpm[order(topcpm$rowsum, decreasing = TRUE),]
topcpm$rowsum <- NULL
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[1:10000,]
topcpm <- t(scale(t(topcpm)))
pdf(file = "Top 10000 expressed genes.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(topcpm, col = col.pan, Rowv=TRUE, scale="column",
            trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Highest expressed genes",
            labRow = FALSE)
dev.off()

getwd()

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

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}



# KEGG PLOTS
setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}
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
         pathway.id="mmu04728", 
         species="mmu", 
         new.signature=FALSE)

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


getwd()
dir.create("another kegg plots")
setwd("another kegg plots")


for (f in rownames(tk_common)){
  plot_pathway(f)
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

### DISEASE ASSOCIATION

if (disease_association == TRUE){
et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

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

pth_pan_down <- pth_pan_down[order(pth_pan_down$Freq, decreasing = TRUE),]
pth_pan_down <- pth_pan_down[seq(1:pathways_set),]
names(pth_pan_down) <- c("Pathway ID", "Frequency")

write.xlsx(go_pan_up, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(go_pan_down, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)

write.xlsx(pth_pan_up, file = "Top Pathways by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(pth_pan_down, file = "Top Pathways by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)
}

### DESEQ2 PART (ADDITIONAL)


ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res <- results(dds, tidy = FALSE )
rld<- rlogTransformation(dds, blind=TRUE)
setwd(results.dir)

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
beep()
beep()
beep()

### COMPARE DESEQ2 AND EDGER

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

deviant$Symbol <- mapIds(org.Mm.eg.db, 
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
### MOTIV ENRICHMENT

beep()

if (motiv == TRUE){
      geneList1 <- et_annot$symbol
      geneLists <- list(geneListName=geneList1)
      motifRankings <- mm9_10kbpAroundTss_motifRanking
      
      motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                      motifAnnot_direct=mm9_direct_motifAnnotation)
      
      motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
      resultsSubset <- motifEnrichmentTable_wGenes_wLogo
      result <- datatable(resultsSubset[,-c("enrichedGenes", "TF_inferred"), with=FALSE], 
                          escape = FALSE, # To show the logo
                          filter="top", options=list(pageLength=5))
      
      
      saveWidget(result, file="TF enrichment.html")
      
      motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
      motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=mm9_direct_motifAnnotation)
      
      signifMotifNames <- motifEnrichmentTable$motif[1:3]
      incidenceMatrix <- getSignificantGenes(geneLists, 
                                             motifRankings,
                                             signifRankingNames=signifMotifNames,
                                             plotCurve=TRUE, maxRank=5000, 
                                             genesFormat="incidMatrix",
                                             method="aprox")$incidMatrix
      
      edges <- melt(incidenceMatrix)
      edges <- edges[which(edges[,3]==1),1:2]
      colnames(edges) <- c("from","to")
      
      
      motifs <- unique(as.character(edges[,1]))
      genes <- unique(as.character(edges[,2]))
      nodes <- data.frame(id=c(motifs, genes),   
                          label=c(motifs, genes),    
                          title=c(motifs, genes), # tooltip 
                          shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                          color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
      
      web <- visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                                     nodesIdSelection = TRUE)
      
      saveWidget(web, file="TF enrichment web.html")
      write.xlsx(motifEnrichmentTable_wGenes, file = "Transcription factor binding motif enrichment.xlsx")
}



### GO BOXPLOTS!
### SOME NUMBERS SHOULD BE FIXED dependent do dataset
dir.create("GO boxplots")
setwd("GO boxplots")

goana_straight <- goana(de = et_annot$entrez, species = "Mm")
go_all_1000 <- topGO(goana_straight, n=1000)
go_all_1000 <- go_all_1000[!(go_all_1000$Ont %in% remove), ] 
go_all_1000$perc = (go_all_1000$DE/go_all_1000$N)*100
go_all_1000 <- subset(go_all_1000, go_all_1000$P.DE < 0.05)
go_all_1000 <- subset(go_all_1000, go_all_1000$perc > 20)
go_all_1000 <- subset(go_all_1000, go_all_1000$N < 100)
go_all_1000 <- subset(go_all_1000, go_all_1000$DE <= 10)
go_all_1000 <- subset(go_all_1000, go_all_1000$DE > 5)
go_all_1000 <- go_all_1000[order(go_all_1000$P.DE, decreasing = TRUE),]
go_all_1000 <- go_all_1000[seq(1:20),]

for (i in 1:nrow(go_all_1000)){
  str <- go_all_1000[i,]
  f <- rownames(str)
  plot.name <- str$Term
  plot.ont <- str$Ont
  plot.perc <- round(str$perc, digits = 2)
  plot.pval <- round(str$P.DE, digits = 6)
  z <- allGO[[f]]
  z <- as.data.frame(z)
  thm <- cpm[(cpm$entrez %in% z$z),]
  rownames(thm) <- thm$Symbol
  thm$Symbol <- NULL
  thm$Name <- NULL
  thm$GOID <- NULL
  thm$term <- NULL
  thm$entrez <- NULL
  colnames(thm) <- sampleTable$condition
  thm <- log2(thm)
  thm <-as.data.frame(t(thm))
  thm$Ñondition <- rownames(thm)
  thm <- melt(thm, id.vars = "Ñondition")
  
  g <- ggplot(thm, aes(x = variable, y = value)) +
    geom_boxplot(data = thm, aes(fill = Ñondition), alpha = 0.5) + 
    scale_x_discrete(name = "Experimental Groups") + 
    scale_y_continuous(name = "logCPM") + 
    theme_bw() +
    geom_signif(comparisons = list(c("tg_2", "tg_3")), map_signif_level = TRUE) +
    #facet_wrap(~variable, scales="free") + 
    theme(axis.title.x=element_blank())
  g <- g + ggtitle(paste(f, ":", plot.name, "\n", "Ontology:", plot.ont, "\n", "Percent of significant genes: ", plot.perc, "%", "\n", "Term p-value:", plot.pval)) + theme(plot.title = element_text(hjust = 0.5, size = 10)) 
  ggsave(filename = paste(i, "_", plot.name, ".pdf", sep = ""), plot = g, width = 20, height = 20, units = "cm")
  
}  


###EDITABLE
### GREP SOME GENES
pattern <- c("lymphocyte antigen")
a <- grep(paste(pattern, collapse = "|"), cpm$Symbol, ignore.case = TRUE)
thm <- cpm[a,]
### REMOVE SOME GENES
#remove <- c("Cdr1")
#cpm <- cpm[!(cpm$Symbol %in% remove), ] 
rownames(thm) <- make.names(thm$Symbol, unique = TRUE)
thm$Symbol <- NULL
thm$Name <- NULL
thm <- t(scale(t(thm)))
thm <- thm[complete.cases(thm),]

pdf(file = "Ly.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(thm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(15,15), lhei=c(2,10), lwid=c(2,6), main = "", ColSideColors = names(col),
          labCol = FALSE)

legend("topright",
       legend = c("Control-1", "Control-3", "Tg-1", "Tg-2", "Tg-3"),
       col = c("yellow", "purple", "green", "red", "blue"),
       lty= 1,
       lwd = 10)
dev.off()

### BOX PLOT
