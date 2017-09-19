source('http://bioconductor.org/biocLite.R')
biocLite("topGO")
library("topGO")
library(biomaRt)
library(org.Dm.eg.db)
df <- read.csv(file = "~/res_GO.csv")
df <- df[complete.cases(df), ]

all_genes <- c(df$log2FoldChange)

names(all_genes) <- df$X

GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
              ID = "Ensembl")


resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
GenTable(GOdata, classicFisher = resultFisher, topNodes = 100)

