library(xlsx)
setwd("~/transcriptomes/reads/intron_retention/fus/cash/")
ref <- read.xlsx("as.events.xlsx",sheetIndex = 1)
ref <- ref[which(ref$padj.fus1.fus2 < 0.05),]

df <- read.delim("tg3vstg2.alldiff.txt")
df <- df[which(df$FDR < 0.05),]
df <- df[which(df$P.Value < 0.05),]
ggplot(data=df) + geom_boxplot(aes(x = SplicingType, y = delta_PSI, fill = SplicingType)) + theme_bw()