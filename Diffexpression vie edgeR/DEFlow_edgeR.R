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

### Experimental design
setwd("~/Fly memory project/ALL")
files <- c("fly_K1.counts", "fly_K2.counts",
            "fly_N1.counts", "fly_N2.counts",
           "fly_F1.counts", "fly_F2.counts",
           "fly_F24.counts","fly_F24_2.counts")
           
           
class <- c("cont", "cont", "cont", "cont", "case", "case", "case", "case")
label <- c("K1", "K2", "N1", "N2", "F1", "F2", "F24_1", "F24_2")


ExpTable <- cbind(files, class, label)
y <- readDGE(files = files, group = class, labels = label)
y <- DGEList(y, group=class)
# keep <- rowSums(cpm(y) > 20) >= 2
# y <- y[keep, , keep.lib.sizes=FALSE]
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

et_annot$logFC <- et_annot$logFC*(-1)

et_annot <- as.data.frame(subset(et_annot, PValue < 0.05))

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

write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes", append = TRUE)
write.xlsx(CountsTable, file = "Results edgeR.xlsx", sheetName = "Counts Table", append = TRUE)

### ANNOTATE RAW COUNTS
raw_counts$symbol <- mapIds(org.Dm.eg.db, 
                             keys=row.names(raw_counts), 
                             column="SYMBOL", 
                             keytype="FLYBASE",
                             multiVals="first")

raw_counts$name <- mapIds(org.Dm.eg.db, 
                           keys=row.names(raw_counts), 
                           column="GENENAME", 
                           keytype="FLYBASE",
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
all_genes <- c(et_annot$logFC)
names(all_genes) <- rownames(et_annot)
go_data <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
              ID = "Ensembl", nodeSize = 5)
go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
go_table <- GenTable(go_data, weightFisher = go_test,
                     orderBy = "weightFisher", ranksOf = "weightFisher",
                     topNodes = sum(score(go_test) < .05))
go_table[grep("neu", go_table$Term), ]


go_id_ra <- go_table[grep("oogenesis", go_table$Term), "GO.ID"]
go_genes_ra <- genesInTerm(go_data, go_id_ra)[[1]]
ensembl <- useMart(host = "sep2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "dmelanogaster_gene_ensembl")

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

log_cpm_ra$symbol <- mapIds(org.Dm.eg.db, 
                            keys=row.names(log_cpm_ra), 
                            column="SYMBOL", 
                            keytype="FLYBASE",
                            multiVals="first")
sel <- rownames(log_cpm_ra)

rownames(log_cpm_ra) <- log_cpm_ra$symbol 
log_cpm_ra$mean_cont <- (log_cpm_ra$cont1+log_cpm_ra$cont2)/2
log_cpm_ra$mean_case <- (log_cpm_ra$case1+log_cpm_ra$case2)/2

g <- ggplot() + geom_boxplot(data = log_cpm_ra, aes(x = rownames(log_cpm_ra), y = log_cpm_ra$mean_case), color = "red") + 
                geom_boxplot(data = log_cpm_ra, aes(x = rownames(log_cpm_ra), y = log_cpm_ra$mean_cont), color = "blue")
g
dev.off()       


