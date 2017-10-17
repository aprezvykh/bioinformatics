library("ggplot2")
library("xlsx")
et_annot_non_filtered <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg2/manifestation_deg_12.csv")
et_annot_non_filtered <- et_annot_non_filtered[complete.cases(et_annot_non_filtered),]
allgenes <- nrow(et_annot_non_filtered)
# et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > 1 & et_annot_non_filtered$PValue < 0.05/allgenes)
setwd("~/GitHub/counts/ALS Mice/new filtering/volcano/")
pdf(file = "Volcano plot12.pdf", width = 10, height = 10)
g = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour="Red")) +
  geom_point(alpha=1, size=1) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
g
dev.off()
