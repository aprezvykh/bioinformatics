col.pan <- colorpanel(100, "blue", "white", "red")
tg12 <- read.xlsx("~/counts/ALS Mice/experimental/results/fg1_w_kegga_w_FDR(final)/fc0.5/Tg-1-Tg-2/Goana GO tests, upreg.xlsx", sheetIndex = 3)
tg23 <- read.xlsx("~/counts/ALS Mice/experimental/results/fg1_w_kegga_w_FDR(final)/fc0.5/Tg-2-Tg-3/Goana GO tests, upreg.xlsx", sheetIndex = 3)
control13 <- read.xlsx("~/counts/ALS Mice/experimental/results/fg1_w_kegga_w_FDR(final)/fc0.5/Control-1-Control-3/Goana GO tests, upreg.xlsx", sheetIndex = 3)
i <- intersect(tg12$NA., tg23$NA.)
x <- intersect(i, control13$NA.)
i <- x
com_tg12 <- tg12[(tg12$NA. %in% i),]
com_tg23 <- tg23[(tg23$NA. %in% i),]
com_cnt13 <- control13[(control13$NA. %in% i),]

df <- data.frame(com_tg12$NA., com_tg12$Term, com_tg12$Ont, com_cnt13$DE ,com_tg12$DE, com_tg23$DE)

names(df) <- c("id", "term", "ont","control1-3", "tg1-2", "tg2-3")

df <- df[which(df$ont == "BP"),]
rownames(df) <- df$term
df$id <- NULL
df$term <- NULL
df$ont <- NULL
df <- as.matrix(df)
df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
pdf(file = "GO HEATMAP COUNT.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(df, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(20,20), lhei=c(2,10), lwid=c(1,1), main = "GO terms")
dev.off()