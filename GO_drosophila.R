biocLite("topGO")
biocLite("AffyLib")
install.packages("AffyLib")
library(package = affyLib, character.only = TRUE)

library("topGO")
df <- read.csv(file = "~/diffexp_reports/res_GO.csv")
df <- df[complete.cases(df), ]

changes <- data.frame(df$X, df$log2FoldChange, df$padj)
geneList <- as.character(df$X)

GOdata <- new("topGOdata", description = "Simple session",
                    ontology = "BP", allGenes = geneList,
                    geneSel = topDiffGenes,
                    nodeSize = 10, annot = annFUN.db,
                    affyLib = affyLib)
