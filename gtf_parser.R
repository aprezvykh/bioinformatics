  library(rtracklayer)
  library(RCurl)
  library(matrixStats)
  library(data.table)
  setwd("~/refgenome/")
  
  
  
gtf <- rtracklayer::import('Caenorhabditis_elegans.WBcel235.90.gtf')
df <- as.data.frame(gtf)


genes <- unique(df$gene_id)
g <- data.frame()
g <- df[which(df$gene_id == genes[12]),]
g <- g[which(g$type != "CDS"),]
gene.gene <- g[which(g$type == "gene"),]
gene.length <- gene.gene$width
gene.length.adjusted <- round(gene.length*0.7, digits = 0)
exons <- g[which(g$type == "exon"),]
exon.length <- sum(exons$width)
exon.length.adjusted <- round(exon.length*0.7, digits = 0)
gene.element.count <- nrow(g)



rtracklayer::export(df_parse, "cel_modified_3-hatch_30.gtf",format = "GTF")


