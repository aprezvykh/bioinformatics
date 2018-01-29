int <- read.xlsx("~/counts/AIKAR.final/filtered diffexpression ET_tt_LFC0,5.xlsx", sheetIndex = 1)
cpm <- read.csv("~/counts/AIKAR.elim.2/results/control_late-aikar_late/cpm.csv")

int <- int[order(int$FDR),]
sub <- int[1:50,]
sub <- int



hm <- cpm[(cpm$X %in% sub$NA.),]
rownames(hm) <- hm$Symbol
hm$Symbol <- NULL
hm$X <- NULL
hm$Name <- NULL
hm$entrez <- NULL
hm$GOID <- NULL
hm$term <- NULL

setwd("~/counts/AIKAR.final/")
mat <- t(scale(t(hm)))

pdf(file = "Heatmap aikar early all.pdf", height = 10, width = 10, family = "Helvetica")
heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=0.9, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), labRow = FALSE)
dev.off()



correl <- cpm[,2:5]
correl <- as.data.frame(correl)
correl <- correl[complete.cases(correl),]
x <- cor(correl, method = "spearman")
pdf(file = "Corellation matrix Spearman LATE.pdf", width = 10, height = 10)
heatmap.2(x, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")

dev.off()



hist(cpm$X10_worm_aikar_late_10.counts)
