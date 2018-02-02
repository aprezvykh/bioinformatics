  library(rtracklayer)
  library(RCurl)
  library(matrixStats)
  library(data.table)
  setwd("~/refgenome/")
  
  
  
gtf <- rtracklayer::import('Caenorhabditis_elegans.WBcel235.90.gtf')
df <- as.data.frame(gtf)


genes <- unique(df$gene_id)
gene.count <- length(genes)
new.gtf <- data.frame()

for(f in genes){
    print(f)
    g <- df[which(df$gene_id == f),]
    df.sub <- g[which(g$type == "gene"),]
    gene.length <- df.sub$width
    if (gene.length >= 4000){
      print("Gene is too long! Skipping")
      
    } else if (gene.length < 4000){
      print("Gene is short! Rbinding!")  
      new.gtf <- rbind(g, new.gtf)

    }
}

      
      #gene.gene <- g[which(g$type == "gene"),]
#gene.length <- gene.gene$width
#gene.length.adjusted <- round(gene.length*0.7, digits = 0)
#exons <- g[which(g$type == "exon"),]
#exon.length <- sum(exons$width)
#exon.length.adjusted <- round(exon.length*0.7, digits = 0)
#gene.element.count <- nrow(g)



rtracklayer::export(new.gtf, "cel_modified_1000nc.gtf",format = "GTF")


