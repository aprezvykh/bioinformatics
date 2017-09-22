source('http://bioconductor.org/biocLite.R')
install.packages("statmod")
install.packages("reshape")
library("edgeR")
library("topGO")
library("org.Mm.eg.db")
library("xlsx")
library("pheatmap")
library("biomaRt")
library("tidyr")
library("knitr")
library("ggplot2")
library(matrixStats)
library("reshape")
### Experimental design
setwd("~/counts_ens/2_late_tg_vs_ctrl_tg/")
files <- c("mouse_tg_1_182.counts", "mouse_tg_1_183.counts", "mouse_tg_1_184.counts", "mouse_tg_1_185.counts", "mouse_tg_1_212.counts", 
           "mouse_tg_3_171.counts", "mouse_tg_3_172.counts", "mouse_tg_3_173.counts", "mouse_tg_3_174.counts", "mouse_tg_3_175.counts")
class <- c("cont", "cont", "cont", "cont", "cont", "case", "case", "case", "case", "case")
label <- c("cont1", "cont2", "cont3", "cont4", "cont5", "case1", "case2", "case3", "case4", "case5")

ExpTable <- cbind(files, class, label)
y <- readDGE(files = files, group = class, labels = label)
y <- DGEList(y, group=class)
keep <- rowSums(cpm(y) > 20) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
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

et_annot$logFC <- et_annot$logFC*(-1)

et_annot <- as.data.frame(subset(et_annot, PValue < 0.05))

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

write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes", append = TRUE)
write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table", append = TRUE)

### ANNOTATE RAW COUNTS
raw_counts$symbol <- mapIds(org.Mm.eg.db, 
                             keys=row.names(raw_counts), 
                             column="SYMBOL", 
                             keytype="ENSEMBL",
                             multiVals="first")

raw_counts$name <- mapIds(org.Mm.eg.db, 
                           keys=row.names(raw_counts), 
                           column="GENENAME", 
                           keytype="ENSEMBL",
                           multiVals="first")
CountsTable <- CountsTable[rownames(et_annot),]

### HEATMAP
select <- et_annot[order(et_annot$logFC),]
select <- et_annot[1:20,]
b <- rownames(select)
z <- cpm(y, log=TRUE, prior.count = 1)
pheatmap(scale(z, center = TRUE)[b,])
dev.off()


### GO
et_annot_subset <- as.data.frame(subset(et_annot, logFC > 1 | logFC < -1))
all_genes <- c(et_annot_subset$logFC)
names(all_genes) <- rownames(et_annot_subset)
go_data <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                0.01, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", 
              ID = "Ensembl", nodeSize = 5)
go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
go_table <- GenTable(go_data, weightFisher = go_test,
                     orderBy = "weightFisher", ranksOf = "weightFisher",
                     topNodes = sum(score(go_test) < .1))


# ENTER HERE YOU INTERESTING PATHWAY ID! 
term <- c("cell maturation")

go_table[grep(paste(term), go_table$Term), ]


go_id_ra <- go_table[grep(paste(term), go_table$Term), "GO.ID"]
go_genes_ra <- genesInTerm(go_data, go_id_ra)[[1]]
ensembl <- useMart(host = "sep2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl")

gene_info_ra <- getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                     "external_gene_name", "description"),
                      filter = "ensembl_gene_id",
                      values = go_genes_ra,
                      mart = ensembl)
stopifnot(go_genes_ra == gene_info_ra$ensembl_gene_id)
gene_info_ra

log_cpm_ra <- log_cpm[go_genes_ra, ]
log_cpm_ra <- data.frame(gene = rownames(log_cpm_ra), log_cpm_ra,
                         stringsAsFactors = FALSE)

log_cpm_ra$symbol <- mapIds(org.Mm.eg.db, 
                            keys=row.names(log_cpm_ra), 
                            column="SYMBOL", 
                            keytype="ENSEMBL",
                            multiVals="first")
sel <- rownames(log_cpm_ra)
pheatmap(scale(log_cpm, scale = FALSE, center = TRUE)[sel,])

rownames(log_cpm_ra) <- log_cpm_ra$symbol

log_cpm_ra$cont_avg <- rowSums(log_cpm_ra[,2:6])/5
log_cpm_ra$case_avg <- rowSums(log_cpm_ra[,7:11])/5
log_cpm_ra <- transform(log_cpm_ra, SD_control=apply(log_cpm_ra[,2:6],1, sd, na.rm = TRUE))
log_cpm_ra <- transform(log_cpm_ra, SD_case=apply(log_cpm_ra[,7:11],1, sd, na.rm = TRUE))
log_cpm_ra$logFC <- log2(log_cpm_ra$case_avg/log_cpm_ra$cont_avg)

for_plot <- data.frame(log_cpm_ra$symbol, log_cpm_ra$cont_avg, log_cpm_ra$case_avg)

g <- ggplot() + geom_boxplot(data = log_cpm_ra, aes(x = rownames(log_cpm_ra), y=cont_avg, fill=SD_control), color = "blue", fill = "white") +
      
                geom_boxplot(data = log_cpm_ra, aes(x = rownames(log_cpm_ra), y = case_avg, fill=SD_control), color = "red", fill = "white")
g

p10 <- ggplot(for_plot, aes(x = log_cpm_ra.symbol, y = log_cpm_ra.cont_avg)) +
  geom_boxplot()
p10
