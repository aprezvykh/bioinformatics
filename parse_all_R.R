#install.packages("factoextra")
library(rtracklayer)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(factoextra)
col.pan <- brewer.pal(n = 9, name = "Set3")
scf.dir <- "~/genomes/comp.transposones/fa"
gtf.dir <- "~/genomes/12_genomes/all_gtf/"
out.dir <- "~/genomes/comp.transposones/out/"
gtf.list <- grep("gtf", list.files(gtf.dir), value = TRUE)
scf.list <- grep("fa", list.files(scf.dir), value = TRUE)
out.list <- grep("out", list.files(out.dir), value = TRUE)

coords <- read.delim("~/genomes/comp.transposones/coords")
coords <- coords[order(coords$spec),]
ir.length <- coords$end-coords$start

unique(gtf$type)
compart <- data.frame(scf.list, gtf.list)
genes.count <- data.frame()
percentages.df <- data.frame()
for (i in seq(1:9)){
    setwd(gtf.dir)
    gtf <- rtracklayer::import(gtf.list[i])
    gtf <- as.data.frame(gtf)
    gtf$seqnames <- as.character(gtf$seqnames)
    sub.test <- gtf[which(gtf$seqnames == coords$chr[i]),]
    sub.test <- sub.test[which(sub.test$start > coords$start[i]),]
    sub.test <- sub.test[which(sub.test$end < coords$end[i]),]
    all.genes <- sub.test[which(sub.test$type == "exon"),]
    sum.length <- sum(all.genes$end-all.genes$start)
    setwd(out.dir)
    out.df <- read.table(out.list[i])
    names(out.df) <- c("score", "div", "del", "ins", 
                   "qseq", "qbegin", "qend", "qleft",
                   "strand", "repeat_id", "class",
                   "sbegin","send", "sleft")
    
    sum.trans.length <- sum(out.df$qend - out.df$qbegin)
    df <- data.frame(coords$spec[i], nrow(sub.test[which(sub.test$type == "gene"),]), sum.length, sum.trans.length)
    print(df)
    genes.count <- rbind(df, genes.count)
}
names(genes.count) <- c("stain", "genes", "sum.gene.length", "sum.trans.length")
genes.count <- genes.count[order(genes.count$stain),]
genes.count$ir.length <- coords$end-coords$start
View(genes.count)

genes.count$trans.idx <- 100*(genes.count$sum.trans.length/(coords$end-coords$start))
genes.count$gene.idx <-100*(genes.count$sum.gene.length/(coords$end-coords$start))


idx <- data.frame(genes.count$gene.idx, genes.count$trans.idx)
idx$names <- genes.count$stain
names(idx) <- c("genes", "trans", "names")
idx <- melt(idx, id.vars = "names")

ggplot(data = idx) + geom_bar(aes(x = names, y = value, color = variable, fill = variable), stat = "identity") + theme_bw()



setwd(out.dir)
tr.for.pca <- data.frame()
for (i in seq(1:9)){
  df <- read.table(out.list[i])
  names(df) <- c("score", "div", "del", "ins", 
                     "qseq", "qbegin", "qend", "qleft",
                     "strand", "repeat_id", "class",
                     "sbegin","send", "sleft")
  df <- as.data.frame(df)
  df$spec <- paste(out.list[i])
  df$qbegin <- df$qbegin/ir.length[i]
  df$qend <- df$qend/ir.length[i]
  sub.tr.for.pca.df <- data.frame(df$class, df$spec, df$score, df$div, df$del, df$ins,
                                  df$qbegin, df$qend,
                                  df$sbegin,
                                  df$send)
  tr.for.pca <- rbind(sub.tr.for.pca.df, tr.for.pca)
}

tr.for.pca$df.spec <- strtrim(tr.for.pca$df.spec, 4)
tr.for.pca$df.sbegin <- gsub("\\(|\\)", "", tr.for.pca$df.sbegin)
tr.for.pca$df.sbegin <- as.numeric(tr.for.pca$df.sbegin)
View(tr.for.pca)
p <- prcomp(tr.for.pca[,3:8],center = TRUE)
ggbiplot(p, groups = tr.for.pca$df.spec, var.axes = TRUE, ellipse = TRUE, choices = c(1,3)) + theme_bw()


plot(p$x, col = as.factor(tr.for.pca$df.spec), lwd = as.factor(tr.for.pca$df.class))

legend("topright",unique(tr.for.pca$df.spec),col=1:length(tr.for.pca$df.spec),pch=1)



dev.off()

dev.off()
ggsave(g, filename = "trans_PCA.pdf",device = "pdf",width = 32, height = 22)


fviz_pca_ind(p,
             col.ind = tr.for.pca$df.spec, # color by groups
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Drosophila species",
             repel = F,
             geom.ind = c("point"),
             habillage = tr.for.pca$df.spec) + 
             scale_shape_manual(values=1:9)




plot(tr.for.pca$df.score)
