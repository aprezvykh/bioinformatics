library("ggplot2")
library("xlsx")
et_annot_non_filtered <- read.xlsx("~/GitHub/counts/ALS Mice/new filtering/tg2-tg3/manifestation_deg.xlsx", sheetIndex = 2)
et_annot_non_filtered <- et_annot_non_filtered[complete.cases(et_annot_non_filtered),]
allgenes <- nrow(et_annot_non_filtered)

moto <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg2-tg3/moto/deg_moto_23.csv")
glia <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg2-tg3/glia/deg_glia_23.csv")
other <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg2-tg3/other/deg_out_23.csv")

et_annot_non_filtered$symbol <- NULL
et_annot_non_filtered$name <- NULL
et_annot_non_filtered$entrez <- NULL

df <- data.frame()
et_annot_non_filtered$keep <- NULL

et_annot_non_filtered$keep_moto <- ifelse(et_annot_non_filtered$NA. %in% moto$NA., "Moto", "")
et_annot_non_filtered$keep_glia <- ifelse(et_annot_non_filtered$NA. %in% glia$NA., "Glia", "")  
et_annot_non_filtered$keep_other <- ifelse(et_annot_non_filtered$NA. %in% other$NA., "Other", "")


et_annot_non_filtered$keep <- paste(et_annot_non_filtered$keep_moto
                                    , et_annot_non_filtered$keep_glia
                                    , et_annot_non_filtered$keep_other)

et_annot_non_filtered$keep_moto <- NULL
et_annot_non_filtered$keep_glia <- NULL
et_annot_non_filtered$keep_other <- NULL
names(et_annot_non_filtered) <- c("NA", "logFC", "logCPM", "PValue", "Tissue")
# et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > 1 & et_annot_non_filtered$PValue < 0.05/allgenes)
setwd("~/GitHub/counts/ALS Mice/new filtering/volcano/")
pdf(file = "Volcano plot23.pdf", width = 10, height = 10)
g = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour=Tissue)) +
  geom_point(alpha=1, size=2) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw()
g
dev.off()
g