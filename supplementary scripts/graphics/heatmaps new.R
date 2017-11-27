library("xlsx")
library("gplots")
library(org.Mm.eg.db)
library(pheatmap)
library(GO.db)
library(ggplot2)
row.names.remove <- c("NA.1")
col.pan <- colorpanel(100, "blue", "white", "red")
cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
cpm <- cpm[complete.cases(cpm),]
rownames(cpm) <- cpm$X
cpm$X <- NULL
#a <- grep("tg", colnames(cpm))
#cpm <- cpm[,a]

tg_12_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/glia/tg12-glia.csv")
tg_12_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/moto/tg12-moto.csv")
tg_12_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/other/tg12-other.csv")

tg_13_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/glia/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/tg13-other.csv")


setwd("~/counts/ALS Mice/experimental/results/all/")
#all
common <- data.frame()
common <- tg_13_glia
common <- rbind(tg_13_moto, common)
common <- rbind(tg_13_others, common)
common$X <- NULL

rownames(common) <- common$NA.
common$NA. <- NULL

colors <- c("yellow", "yellow", "yellow", "yellow", "yellow", 
            "purple", "purple", "purple", "purple", "purple", 
            "green", "green", "green", "green", "green", 
            "red", "red", "red", "red", 
            "blue", "blue", "blue", "blue", "blue")
sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                     'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                     'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                     'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')
### 23
hm <- cpm[(rownames(cpm) %in% tg_13_moto$X),]
hm$rowsum <- rowSums(hm)
hm <- hm[order(hm$rowsum, decreasing = TRUE),]
hm <- hm[seq(1:50),]
hm$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(hm), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
rownames(hm) <- hm$Symbol

hm$Symbol <- NULL
hm$rowsum <- NULL
colnames(hm) <- sampleCondition
hm <- t(scale(t(hm)))
names(hm)
#pdf(file = "Top 50 Tg1-Tg3 Others.pdf", width = 12, height = 17, family = "Helvetica")

heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.5, cexCol=1.5, density.info="none",
          margin=c(15,11), lhei=c(2,10), lwid=c(2,6), main = "Top 50 Tg1-Tg3 Others", ColSideColors = colors)
legend("topright",
       legend = c("Tg-1", "Tg-2", "Tg-3"),
       col = c("yellow","green", "red"),
       lty= 1,
       lwd = 10)
dev.off()





getwd()




##12
tg_12_glia$type <- paste("yellow")
tg_12_moto$type <- paste("green")
tg_12_others$type <- paste("red")

common <- data.frame()
common <- tg_12_glia
common <- rbind(tg_12_moto, common)
common <- rbind(tg_12_others, common)
common$X <- NULL

rownames(common) <- common$X

hm <- cpm[(rownames(cpm) %in% rownames(common)),]
hm$sum <- rowSums(hm)
hm <- hm[order(hm$sum, decreasing = TRUE),]
hm <- hm[seq(1:50),]

hm$Symbol <- mapIds(org.Mm.eg.db, 
                    keys=row.names(hm), 
                    column="SYMBOL", 
                    keytype="ENSEMBL",
                    multiVals="first")
rownames(hm) <- hm$Symbol
hm$Symbol <- NULL
hm$sum <- NULL
hm <- as.matrix(hm)
hm <- t(scale(t(hm)))
pdf(file = "Top 50 Tg1-Tg2.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=0.9, cexCol=1.4, density.info="none",
          margin=c(15,11), lhei=c(2,10), lwid=c(2,6), RowSideColors = common$type)

legend("topright",
       legend = c("glia", "Moto", "Others"),
       col = c("red", "green", "blue"),
       lty= 1,
       lwd = 10
)

dev.off()


### GO HEATMAP
row.names.remove <- c("NA.1")
cpm <- read.csv("~/counts/ALS Mice/new filtering/experimental overall logCPM.csv")
cpm <- cpm[complete.cases(cpm),]
cpm <- cpm[!(cpm$X %in% row.names.remove),] 
rownames(cpm) <- cpm$X
cpm$X <- NULL

cpm$Symbol <- mapIds(org.Mm.eg.db, 
                    keys=row.names(cpm), 
                    column="SYMBOL", 
                    keytype="ENSEMBL",
                    multiVals="first")

cpm$Name <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="GENENAME", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$GOID <-     mapIds(org.Mm.eg.db, 
                            keys=row.names(cpm), 
                            column="GO", 
                            keytype="ENSEMBL",
                            multiVals="first")

cpm$term <- mapIds(GO.db, 
                        keys=cpm$GOID, 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")

cpm$term <- as.character(cpm$term)

freq <- as.data.frame(table(unlist(cpm$term)))
freq <- freq[order(freq$Freq, decreasing = TRUE),]

a <- grep("mRNA splicing, via spliceosome", cpm$term)
hm <- cpm[a,]
rownames(hm) <- hm$Symbol
hm$Symbol <- NULL
hm$Name <- NULL
hm$term <- NULL
hm$GOID <- NULL
hm <- t(scale(t(hm)))
pdf(file = "mRNA splicing, via spliceosome.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(11,11), lhei=c(2,10), lwid=c(2,6), main = "mRNA splicing, via spliceosome")
dev.off()
cpm$GOID <- NULL
cpm$term <- NULL

##for CD
a <- grep("glutamate", cpm$Name, ignore.case = TRUE)
hm <- cpm[a,]
#a <- grep("antigen", hm$Name)
#hm <- hm[a,]
#remove <- c("Cdr1")
#hm <- hm[!(hm$Symbol %in% remove), ] 
rownames(hm) <- hm$Symbol
hm$Symbol <- NULL
hm$Name <- NULL

hm <- t(scale(t(hm)))


pdf(file = "glutamate.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(15,11), lhei=c(2,10), lwid=c(2,6), main = "glutamate")
dev.off()
getwd()

##VOLCANO
pdf(file = "Volcano 1-3.pdf", width = 10, height = 10, family = "Helvetica")
g = ggplot(data=common, aes(x=logFC, y=-log10(PValue), colour=type)) +
  geom_point(alpha=1, size=2) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw()
g

dev.off()
grep("Ras", cpm$Symbol, value = TRUE)
