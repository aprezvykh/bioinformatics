library(gplots)
library(pheatmap)
library(xlsx)
control3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Control-3/kegg_down.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-1/kegg_down.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-2/kegg_down.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-3/kegg_down.csv")
sel <- read.xlsx("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/kegga.xlsx", sheetIndex = 3)

sel1 <- control3[(control3$Pathway %in% sel$Pathway),]
sel2 <- tg1[(tg1$Pathway %in% sel$Pathway),]
sel3 <- tg2[(tg2$Pathway %in% sel$Pathway),]
sel4 <- tg3[(tg3$Pathway %in% sel$Pathway),]
df <- data.frame(sel1$Pathway, sel2$DE, sel3$DE, sel4$DE)
names(df) <- c("term", "Tg-1", "Tg-2", "Tg-3")
rownames(df) <- df$term
df$term <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Significantly downregulated KEGG terms in motoneurons.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
         trace="none", dendrogram="none", cexRow=5, cexCol=5, density.info="none",
         margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "Significantly downregulated KEGG terms in motoneurons", border_color = NA,
         cluster_rows=T, cluster_cols=F, fontsize_row = 12, fontsize_col = 12)
dev.off()

sel
