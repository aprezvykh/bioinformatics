source('http://bioconductor.org/biocLite.R')
biocLite("topGO")
biocLite("AffyLib")
install.packages("AffyLib")
library(package = affyLib, character.only = TRUE)
library("topGO")
df <- read.csv(file = "~/diffexpr_reports/res_GO.csv")
df <- df[complete.cases(df), ]

vec <- c(df$padj)
names(vec) <- as.character(df$X)
tdg <- topDiffGenes(vec)

dros_BP <- new("topGOdata",
                    description = "Simple", ontology = "BP",
                    allGenes = vec, geneSel = tdg,
                    nodeSize = 100,
                    annot = annFUN.db,
                    mapping = "Org.Dm.eg.db")
data(geneList)

typeof(data)
