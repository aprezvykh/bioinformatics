## GO ENRICHMENT WITHOUT CONSIDERING FOLD CHANGES! ##
   ## REQUIRES DATASET WITH ENSEMBL/FLYBASE ID! ##
source('http://bioconductor.org/biocLite.R')
biocLite("topGO")
library("topGO")
library(biomaRt)
library(org.Dm.eg.db)

logfc_low <- -1
logfc_high <- 1
basemean <- 100
topNodes <- 1000

df <- read.csv(file = "~/res_GO.csv")
df <- df[complete.cases(df), ]
sub_low <- as.data.frame(subset(df, log2FoldChange < logfc_low))
sub_high <- as.data.frame(subset(df, log2FoldChange > logfc_high))

myGOBP <- function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  gores <- as.data.frame(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
  gores <- gores[order(gores$Annotated, decreasing = TRUE), ]
  return(gores)
  
}
myGOMF <- function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p < 
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  gores <- as.data.frame(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
  gores <- gores[order(gores$Annotated, decreasing = TRUE), ]
  return(gores)
  
}
myGOCC <- function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(p) p < 
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  gores <- as.data.frame(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
  gores <- gores[order(gores$Annotated, decreasing = TRUE), ]
  return(gores)
  
}



write.xlsx(myGOBP(sub_low, topNodes), file = "GO_low.xlsx", sheetName = "low_bp", append = TRUE)
write.xlsx(myGOMF(sub_low, topNodes), file = "GO_low.xlsx", sheetName = "low_mf", append = TRUE)
write.xlsx(myGOCC(sub_low, topNodes), file = "GO_low.xlsx", sheetName = "low_cc", append = TRUE)
write.xlsx(myGOBP(sub_high, topNodes), file = "GO_high.xlsx", sheetName = "high_bp", append = TRUE)
write.xlsx(myGOMF(sub_high, topNodes), file = "GO_high.xlsx", sheetName = "high_mf", append = TRUE)
write.xlsx(myGOCC(sub_high, topNodes), file = "GO_high.xlsx", sheetName = "high_cc", append = TRUE)

write.xlsx(myGOBP(df, topNodes), file = "GO_all.xlsx", sheetName = "all_bp", append = TRUE)


help(topGO)



