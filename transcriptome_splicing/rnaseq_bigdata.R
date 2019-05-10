source("https://bioconductor.org/biocLite.R")
devtools::install_github("stephenturner/annotables")

library(annotables)
library(ggplot2)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(plyr)
library(dplyr)
library(tidyr)
library(org.Mm.eg.db)
library(GO.db)

annotate.rsem.rownames <- function(rsem){
  trans <- rownames(rsem)
  trans <- data.frame(trans)
  trans <- data.frame(separate(data = trans, col = trans, into = c("enstrans", "gene.name"), sep = "_"))
  trans <- trans %>% 
    dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))
  
  trans$symbol <- mapIds(org.Mm.eg.db, 
                         keys=as.character(trans$ensgene), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
  
  trans$name <- mapIds(org.Mm.eg.db, 
                       keys=as.character(trans$ensgene), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")
  trans$entrez <- mapIds(org.Mm.eg.db, 
                         keys=as.character(trans$ensgene), 
                         column="ENTREZID", 
                         keytype="ENSEMBL",
                         multiVals="first")
  
  trans$GOID <- mapIds(org.Mm.eg.db, 
                       keys=as.character(trans$ensgene), 
                       column="GO", 
                       keytype="ENSEMBL",
                       multiVals="first")
  
  trans$term <- mapIds(GO.db, 
                       keys=trans$GOID, 
                       column="TERM", 
                       keytype="GOID",
                       multiVals="first")
  
  trans$term <- as.character(trans$term)
  return(trans)
}

exact.test.and.cutoff <- function(p,edger.obj){
  et <- exactTest(edger.obj, pair = p)
  tt <- data.frame(topTags(et, n = nrow(df)))
  tt <- tt[tt$FDR<0.05,]
  tt.ann <- annotate.rsem.rownames(tt)
  annotated.et <- data.frame(tt, tt.ann)
  return(annotated.et)
}

col.pan <- colorpanel(100, "blue", "white", "red")
df <- read.csv("~/transcriptomes/reads/intron_retention/all.rsem/all+pfn.csv")
rownames(df) <- df$X
df$X <- NULL
model <- c("pfn", "pfn", "pfn", "pfn", "pfn", 
           "pfn", "pfn", "pfn", 
           "pfn", "pfn", "pfn", "pfn", "pfn", 
           "pfn", "pfn", "pfn", 
           "fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus",
           "tdp", "tdp", "tdp", "tdp", 
           "tdp", "tdp", "tdp", "tdp", 
           "sod_moto", "sod_moto",
           "sod_moto", "sod_moto",
           "sod_glia", "sod_glia", 
           "sod_glia", "sod_glia", "sod_glia", 
           "sod_glia", "sod_glia", "sod_glia", "sod_glia", 
           "sod_glia", "sod_glia")

age <- c(190, 190, 250, 250, 200,
         50, 50, 50,
         190, 190, 250, 250, 200,
         50, 50, 50,
         60, 60, 60, 60, 60, 
         120, 120, 120, 120, 120, 
         60, 60, 60, 60, 60, 
         90, 90, 90, 90, 
         120, 120, 120, 120, 120,
         38, 38, 38, 38, 
         38, 38, 38, 38, 
         90, 90, 
         90, 90, 
         65, 65, 
         100, 100, 100, 
         150, 150, 150, 150, 
         65, 65)
transgenity <- c("tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", 
                 "wt", "wt", 
                 "tg", "tg", 
                 "tg", "tg", 
                 "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "wt", "wt") 

sample.id <- c(1,2,3,4,5,
               6,7,8,
               9,10,11,12,13,
               14,15,16,
               1,2,3,4,5,
               6,7,8,9,10,
               11,12,13,14,15,
               16,17,18,19,
               20,21,22,23,25,
               1,2,3,4,
               5,6,7,8,
               1,2,
               3,4,
               1,2,
               3,4,5,
               6,7,8,9,
               10,11)

age_fac <- vector()

for(f in age){
  if(f<70){
    age_fac <- append("young", age_fac)
  } else if(f>70 & f<100){
    age_fac <- append("middle", age_fac)
  } else {
    age_fac <- append("old", age_fac)
  }
}

age_fac <- rev(age_fac)
metadata <- data.frame(model = model,
                        age = age,
                        transgenity = transgenity,
                        index = sample.id,
                        age_fac = age_fac)

metadata$trans.index <- seq(1:nrow(metadata))

metadata$final <- paste(model, "_", transgenity, "_", age, "_", sample.id, sep = "")
metadata$expgroup <- paste(model, "_",transgenity, "_", age_fac ,sep = "")
names(df) <- metadata$final


rsem.anotation <- annotate.rsem.rownames(df)
rsem <- df



rownames(rsem) <- unlist(lapply(strsplit(rownames(df), "_"), function(x)x[[1]]))

plot.gene <- function(x,name= "", comp = NULL, height = NULL, msl = "", xl = "", yl = "", fs = 12){
    dat <- rsem[rownames(rsem) %in% rsem.anotation[grep(x, rsem.anotation$ensgene),]$enstrans,]
    dat <- data.frame(t(dat))
    da <- data.frame(rowSums(dat))
    names(da) <- c("gene")
    da <- bind_cols(metadata,da)
    ggplot(data=da, aes(x = expgroup, y = gene, fill = transgenity)) + geom_boxplot(aes(fill = model, fatten = transgenity), coef = 6) + 
      ylab("Counts per million (CPM)") + 
      geom_signif(comparisons = comp, map_signif_level = F, y_position = height,annotations = msl) + 
      xlab(xl) + 
      ylab(yl) + 
      theme(axis.title.x = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black")) + 
      theme(text=element_text(family="Liberation Serif", face = "bold", size=fs)) + 
      theme(legend.position="none") + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

}



plot.gene("ENSMUSG00000026864")
