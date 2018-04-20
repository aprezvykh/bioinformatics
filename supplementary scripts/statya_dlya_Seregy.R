library(xlsx)
library(gplots)
library(RColorBrewer)
library(ggbiplot)
library(Rtsne)
col.pan <- colorpanel(100, "blue", "white", "red")
df <- read.xlsx("~/cpm_for_correlation.xlsx", sheetIndex = 3, header = TRUE)
rownames(df) <- df$NA.
df$NA. <- NULL


c <- cor(t(df), method = "pearson")
pdf(file = "3_serg_pearson.pdf", width = 10, height = 10, family = "Helvetica")
heatmap.2(c, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Pearson corellation")
dev.off()

c <- cor(t(df), method = "spearman")
pdf(file = "3_serg_spearman.pdf", width = 10, height = 10, family = "Helvetica")
heatmap.2(c, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "spearman corellation")
dev.off()









heatmap.2(t(scale(dt)), col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")


p <- prcomp(t(df), center = TRUE, scale. = TRUE)

ggbiplot(p, var.axes = FALSE, labels = colnames(df))



t <- Rtsne(df, verbose = TRUE, check_duplicates = FALSE)

df <- as.data.frame(t$Y)
fit <- kmeans(df, 8)
clusplot(df, fit$cluster, color=TRUE, shade = TRUE, lines = 2)
