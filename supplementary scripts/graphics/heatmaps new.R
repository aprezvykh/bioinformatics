library("xlsx")
library("gplots")
library(org.Mm.eg.db)
library(pheatmap)
library(GO.db)
library(ggplot2)
row.names.remove <- c("NA.1")
col.pan <- colorpanel(100, "blue", "white", "red")
cpm <- read.csv("~/counts/ALS Mice/new filtering/experimental overall logCPM.csv")
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
cpm <- cpm[complete.cases(cpm),]
rownames(cpm) <- cpm$X
cpm$X <- NULL
a <- grep("tg", colnames(cpm))
cpm <- cpm[,a]

tg_12_glia <- read.csv("~/counts/ALS Mice/new filtering/tg1-tg2/glia/deg_glia_12.csv")
tg_12_moto <- read.csv("~/counts/ALS Mice/new filtering/tg1-tg2/moto/deg_moto_12.csv")
tg_12_others <- read.csv("~/counts/ALS Mice/new filtering/tg1-tg2/other/deg_out_12.csv")

tg_23_glia <- read.csv("~/counts/ALS Mice/new filtering/tg2-tg3/glia/deg_glia_23.csv")
tg_23_moto <- read.csv("~/counts/ALS Mice/new filtering/tg2-tg3/moto/deg_moto_23.csv")
tg_23_others <- read.csv("~/counts/ALS Mice/new filtering/tg2-tg3/other/deg_out_23.csv")

tg_23_glia$type <- paste("glia")
tg_23_moto$type <- paste("moto")
tg_23_others$type <- paste("others")

setwd("~/counts/ALS Mice/experimental/results/all/")
#all
common <- data.frame()
common <- tg_23_glia
common <- rbind(tg_23_moto, common)
common <- rbind(tg_23_others, common)
common$X <- NULL

rownames(common) <- common$NA.
common$NA. <- NULL


### 23
hm <- cpm[(rownames(cpm) %in% tg_23_others$NA.),]
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
hm <- t(scale(t(hm)))
names(hm)
pdf(file = "Top 50 Tg2-Tg3 Others.pdf", width = 12, height = 17, family = "Helvetica")

heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.5, cexCol=1.5, density.info="none",
          margin=c(15,11), lhei=c(2,10), lwid=c(2,6), main = "Top 50 Tg2-Tg3 Others")
dev.off()





getwd()




##12
tg_12_glia$type <- paste("red")
tg_12_moto$type <- paste("green")
tg_12_others$type <- paste("blue")

common <- data.frame()
common <- tg_12_glia
common <- rbind(tg_12_moto, common)
common <- rbind(tg_12_others, common)
common$X <- NULL

rownames(common) <- common$NA.
common$NA. <- NULL

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
pdf(file = "Volcano 2-3.pdf", width = 10, height = 10, family = "Helvetica")
g = ggplot(data=common, aes(x=logFC, y=-log10(PValue), colour=type)) +
  geom_point(alpha=1, size=2) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw()
g

dev.off()


