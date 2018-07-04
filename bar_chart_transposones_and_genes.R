library(rtracklayer)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(factoextra)
library(stringr)
###parse coords
coords <- read.delim("~/genomes/Drosophila.IR.flanks.no.genes/fa/script_normalized_length.sh", header = FALSE)
coords$V1 <- as.character(coords$V1)

new.table <- data.frame()
for (f in 1:nrow(coords)){
    fs <- strsplit(coords$V1[f], " ")
    spec <- strsplit(fs[[1]][3], "/")[[1]][8]
    spec <- tolower(spec)
    scaf <- strsplit(fs[[1]][4], ":")[[1]][1]
    start <- strsplit(strsplit(fs[[1]][4], "-")[[1]][1],":")[[1]][2]
    end <- strsplit(fs[[1]][4],"-")[[1]][2]
    coords.1 <- data.frame(spec,scaf,start,end)
    names(coords.1) <- c("spec", "reg", "start", "end")
    new.table <- rbind(coords.1, new.table)
}

col.pan <- brewer.pal(n = 9, name = "Set3")
scf.dir <- "~/genomes/Drosophila.IR.flanks.no.genes/fa_only/"
gtf.dir <- "~/genomes/12_genomes/all_gtf/"
out.dir <- "~/genomes/Drosophila.IR.flanks.no.genes/out/"
gtf.list <- grep(".gtf", list.files(gtf.dir), value = TRUE)
scf.list <- grep(".fa", list.files(scf.dir), value = TRUE)
out.list <- grep(".out", list.files(out.dir), value = TRUE)

fa.length <- read.table("~/genomes/Drosophila.IR.flanks.no.genes/fa/fa.length", header = F)
fa.length <- fa.length$V2
ir.length <- fa.length


setwd(out.dir)
tr.for.pca <- data.frame(stringsAsFactors = FALSE)
for (i in seq(1:12)){
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
tr.for.pca$df.class <- as.character(tr.for.pca$df.class)

tr.for.pca[grepl("DNA", tr.for.pca$df.class),]$df.class <- "DNA"
tr.for.pca[grepl("LINE", tr.for.pca$df.class),]$df.class <- "LINE"
tr.for.pca[grepl("LTR", tr.for.pca$df.class),]$df.class <- "LTR"
tr.for.pca[grepl("Simple", tr.for.pca$df.class),]$df.class <- "Satellite"
tr.for.pca[grepl("Low", tr.for.pca$df.class),]$df.class <- "Satellite"


tr.for.pca
tr.for.pca <- tr.for.pca[!grepl("Unknown", tr.for.pca$df.class),]
#tr.for.pca <- tr.for.pca[!grepl("Satellite", tr.for.pca$df.class),]



#new.table <- read.table("~/genomes/Drosophila.IR.flanks.no.genes/fa/coords")
#names(new.table) <- c("spec", "reg", "start", "end")

setwd(gtf.dir)
big.sus <- data.frame()
for(f in gtf.list){
  print(f)
  gtf <- as.data.frame(rtracklayer::import(f))
  spec <- strtrim(f, 4)
  regs <- new.table[which(new.table$spec== spec),]
  a <- 0
  sgm <- 0
  c <- 0
  for (z in 1:nrow(regs)){
    sub.gtf <- as.data.frame(gtf[which(gtf$seqnames == as.character(regs[z,2])),])
    sub.gtf <- sub.gtf[which(sub.gtf$start > as.numeric(as.character(regs[z,3]))),]
    sub.gtf <- sub.gtf[which(sub.gtf$start < as.numeric(as.character(regs[z,4]))),]
    b <- sub.gtf[which(sub.gtf$type == "gene"),]
    sgm <- abs(sum(b$end-b$start))
    a <- a + nrow(sub.gtf[which(sub.gtf$type == "gene"),])
    c <- sgm + c
    print(paste(nrow(b), "is not final!", sep = " "))
  }
  print(paste(a, "is final!", sep = " "))
  sus <- data.frame(spec, a, c)
  big.sus <- rbind(sus, big.sus)
}

 big.sus <- big.sus[order(big.sus$spec, decreasing = TRUE),]
big.sus$ir.length <- ir.length
big.sus$gene.perc <- 100*(big.sus$c/big.sus$ir.length)


sub_big <- data.frame()
for (f in unique(tr.for.pca$df.spec)){
  print(f)
  sub <- tr.for.pca[which(tr.for.pca$df.spec == f),]
  frq.sub <- data.frame()
  for(z in unique(sub$df.class)){
    df.sub <- sub[which(sub$df.class == z),]
    ss <- sum(df.sub$df.qend - df.sub$df.qbegin)
    freqs <- data.frame(z,ss)
    freqs$ss <- 100*freqs$ss
    frq.sub <- rbind(freqs, frq.sub)
    print("additional row Binded!")
  }
  dft <- data.frame()
  dft <- data.frame("Genes",big.sus[which(big.sus$spec == f),]$gene.perc, f)
  names(dft) <- c("z","ss", "spec")
  frq.sub$spec <- f
  frq.sub <- rbind(dft, frq.sub)
  sub_big <- rbind(frq.sub, sub_big)
  print("Binded!")
}

#add unannotated 

names(sub_big) <- c("Feature", "freq", "spec")

for(f in unique(sub_big$spec)){
  zz <- sub_big[which(sub_big$spec == f),]
  unc <- 100-sum(zz$freq)
  row <- data.frame("Uncharacterized", unc, as.character(f))
  names(row) <- c("rep", "freq", "spec")
  sub_big <- rbind(row, sub_big)
}

setwd("~/genomes/Paper.Pictures/")
pdf("pic_2_feature_barplot.pdf", width = 16, height = 12, family = "Helvetica")
ggplot(data=sub_big) + 
  geom_bar(aes(x = spec, y = freq, color = Feature, fill = Feature), stat = "identity") + 
  theme_bw() + 
  scale_x_discrete("Species") + 
  scale_y_continuous("Genomic feature")
dev.off()


fa.length <- read.table("~/genomes/Drosophila.IR.flanks.no.genes/fa/fa.length", header = F)
fa.length$V1 <- strtrim(fa.length$V1, 4)

pdf("pic_2_intergenic_length.pdf", width = 16, height = 12, family = "Helvetica")
ggplot(data = fa.length) + geom_bar(aes(x = V1, y = V2), stat = "identity", color = "blue", fill = "blue") + 
  theme_bw() +    
  scale_x_discrete("Species") + 
  scale_y_continuous(breaks = c(0,50000,100000,150000), labels = c("0kb","50kb","100kb","150kb"),name = "Intergenic region length")
dev.off()

sss <- as.matrix(sub_big)
xtabs(formula = rep~spec, data = sss)


##diversity barplot
big.s <- data.frame()
for(f in unique(tr.for.pca$df.spec)){
  sub <- tr.for.pca[which(tr.for.pca$df.spec == f),]
#  sub <- sub[!grepl("Satellite", sub$df.class),]
  s <- data.frame(mean(sub$df.div), f, nrow(sub))
  names(s) <- c("div", "spec", "num")
  big.s <- rbind(s, big.s)
}



ggplot(data = big.s) + 
  geom_bar(aes(x = spec, y = div),stat="identity",color="blue",fill="blue") + 
  theme_bw()
