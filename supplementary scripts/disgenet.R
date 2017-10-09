library("xlsx")
library("dplyr")
library(data.table)
moto <- read.xlsx('~/GitHub/counts/ALS Mice/motoneurons marker.xlsx', header = TRUE, sheetIndex = 1)
moto$NA. <- NULL
rownames(moto) <- moto$X
up <- subset(moto, moto$logFC > 0)
down <- subset(moto, moto$logFC < 0)

table <- read.delim(file = '~/GitHub/data/curated_gene_disease_associations.tsv')
table <- as.data.frame(table)
test <- up[seq(1:3),]

newdf <- data.frame()

for (f in moto$symbol){
sub <- NULL
sub <- table[grepl(paste(f), table$geneSymbol, ignore.case = TRUE),]
sub <- sub[order(sub$score, decreasing = TRUE),]
sub <- sub[seq(1:3),]
sub <- as.data.frame(sub$diseaseName)
sub <- transpose(sub)
sub$gene <- paste(f)
newdf <- rbind(sub, newdf)
}

