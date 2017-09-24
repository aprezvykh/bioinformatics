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


### Reactome
foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "fly", minGSSize=70, readable = TRUE )
head(as.data.frame(x))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=10,  font.size = 9)
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges, categorySize="pvalue", showCategory = 10, fixed = FALSE)
dev.off()

