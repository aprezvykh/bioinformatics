library(rtracklayer)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(factoextra)
library(stringr)
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


new.table <- read.table("~/genomes/Drosophila.IR.flanks.no.genes/fa/coords")
names(new.table) <- c("spec", "reg", "start", "end")

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
    sub.gtf <- sub.gtf[which(sub.gtf$start > regs[z,3]),]
    sub.gtf <- sub.gtf[which(sub.gtf$start < regs[z,4]),]
    b <- sub.gtf[which(sub.gtf$type == "gene"),]
    sgm <- sum(b$width)
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
  dft <- data.frame("Genes",big.sus[which(big.sus$spec == f),]$gene.perc/100, f)
  names(dft) <- c("z","ss", "spec")
  frq.sub$spec <- f
  frq.sub <- rbind(dft, frq.sub)
  sub_big <- rbind(frq.sub, sub_big)
  print("Binded!")
}

#add unannotated 
for(f in unique(sub_big$spec)){
  oo <- sub_big[which(sub_big$spec == f),]
  print(sum(oo$freq))  
}



names(sub_big) <- c("rep", "freq", "spec")

ggplot(data=sub_big) + geom_bar(aes(x = spec, y = freq, color = rep, fill = rep), stat = "identity") + theme_bw()



