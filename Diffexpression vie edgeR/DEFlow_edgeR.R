source('http://bioconductor.org/biocLite.R')
install.packages("statmod")
library("edgeR")
library("topGO")
library("org.Dm.eg.db")
library("xlsx")
library("pheatmap")
library("biomaRt")
library("tidyr")
library("knitr")
library("ggplot2")
library(clusterProfiler)
library(reactome.db)
library(ReactomePA)

cpm_tr <- 1
low_logfc <- -1
high_logfc <- 1
pval_tr <- 0.05

### Experimental design
setwd("~/Fly memory project/K_vs_F24/")
files <- c("fly_K1.counts", "fly_K2.counts",
           "fly_F24.counts","fly_F24_2.counts")

class <- c("cont", "cont", "case", "case")
label <- c("K1", "K2", "F24_1", "F24_2")
ExpTable <- cbind(files, class, label)

y <- readDGE(files = files, group = class, labels = label)
y <- DGEList(y, group=class)
normalized_lib_sizes <- calcNormFactors(y)
log_cpm <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
CountsTable <- as.data.frame(y$counts)
raw_counts <- as.data.frame(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
top <- as.data.frame(topTags(et))
et_annot <- as.data.frame(et$table)
et_annot_non_filtered <- as.data.frame(et$table)



### ANNOTATE LogFC

columns(org.Dm.eg.db)
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

et_annot$logFC <- et_annot$logFC*(-1)

### Simple summary
all <- nrow(raw_counts)
allpadj <- sum(et_annot$PValue < pval_tr, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > high_logfc, na.rm=TRUE)
down <- sum(et_annot$logFC < low_logfc, na.rm=TRUE)
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


et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_tr))
et_annot <- as.data.frame(subset(et_annot, PValue < pval_tr))
et_annot <- as.data.frame(subset(et_annot, logFC > high_logfc | logFC < low_logfc))

write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)
write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table, logCPM>1", append = TRUE)

