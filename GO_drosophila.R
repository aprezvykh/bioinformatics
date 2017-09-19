## GO ENRICHMENT WITHOUT CONSIDERING FOLD CHANGES! ##
   ## REQUIRES DATASET WITH ENSEMBL/FLYBASE ID! ##

logfc_low <- -0.5
logfc_high <- 0.5
basemean <- 100
source('http://bioconductor.org/biocLite.R')
biocLite("topGO")
library("topGO")
library(biomaRt)
library(org.Dm.eg.db)
df <- read.csv(file = "~/res_GO.csv")
df <- df[complete.cases(df), ]
sub_low <- as.data.frame(subset(df, log2FoldChange < logfc_low))
sub_high <- as.data.frame(subset(df, log2FoldChange > logfc_high))

myGOBP <-   function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  return(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
}

myGOMF <-   function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p < 
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  return(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
}

myGOÑÑ <-   function(mydf, topNodes){
  all_genes <- c(mydf$log2FoldChange)
  names(all_genes) <- mydf$X
  GOdata <- new("topGOdata", ontology = "ÑÑ", allGenes = all_genes, geneSel = function(p) p < 
                  0.01, description = "Test", annot = annFUN.org, mapping = "org.Dm.eg.db", 
                ID = "Ensembl")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  return(GenTable(GOdata, classicFisher = resultFisher, topNodes = topNodes))
}



