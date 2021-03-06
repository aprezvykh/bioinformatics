print("sas")
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(dplyr)
library(gridExtra)
library(GenomicAlignments)
library(lattice)
library(ggbiplot)
library(mgcv)
library(data.table)
setwd("~/genomes/IR.w.flanks.20kb/")
df <- read.table("rmask.out")

names(df) <- c("score", "div", "del", "ins", 
               "qseq", "qbegin", "qend", "qleft",
               "strand", "repeat_id", "class",
               "sbegin","send", "sleft")


df$sbegin <- gsub("\\(|\\)", "", df$sbegin)
df$sleft <- gsub("\\(|\\)", "", df$sleft)


df$sbegin <- as.numeric(df$sbegin)
df$send <- as.numeric(df$send)
df$sleft <- as.numeric(df$sleft)
df$sum <- df$sbegin + df$send + df$sleft
df$perc <- 100*((df$send - df$sbegin)/df$sum)

df <- data.frame(df)
df$class <- as.character(df$class)

unique(df$class)

df[which(df$class == 'Unknown'),] <- NULL


df[grepl("DNA", df$class),]$class <- "DNA"
df[grepl("LINE", df$class),]$class <- "LINE"
df[grepl("LTR", df$class),]$class <- "LTR"
df[grepl("Simple", df$class),]$class <- "Satellite"
df[grepl("Low", df$class),]$class <- "Satellite"

df$sum.inv <- ifelse(test = df$perc > 0, "Straight", "Inverted")
g.start <- 1740222
g.end <- 1953040
gff <- readGFF("dvir_IR_flanks.fasta.out.gff")
df <- cbind(gff,df)
param <- ScanBamParam(which=GRanges("scaffold_13050", IRanges(start = g.start,end = g.end)))
read.bam.coverage <- function(x){
  bam <- readGAlignments(x, param = param)
  cov <- coverage(bam)
  xcov <- as.numeric(cov$`scaffold_13050`)
  xcov <- as.data.frame(xcov)
  xcov$idx <- seq(1:nrow(xcov))
  xcov <- xcov[which(xcov$idx < 1953041),]
  xcov <- xcov[which(xcov$idx > 1740221),]
  xcov$idx <- seq(1:nrow(xcov))
  return(xcov)
}
read.bam.coverage.no.param <- function(x){
  bam <- readGAlignments(x)
  cov <- coverage(bam)
  xcov <- as.numeric(cov$`scaffold_13050:1740222-1953040`)
  xcov <- as.data.frame(xcov)
  xcov$idx <- seq(1:nrow(xcov))
  return(xcov)
}


transcripts <- read.bam.coverage.no.param("transcripts.sorted.bam")
transcripts$xcov <- ifelse(transcripts$xcov>0, 1, 0)

h3k9 <- read.bam.coverage("/mnt/raid/illumina/ftp_data/Rezvykh.Chip/H3K9_1.sorted.bam")
pol2 <- read.bam.coverage("/mnt/raid/illumina/ftp_data/Rezvykh.Chip/Pol_1.sorted.bam")
pi <- read.bam.coverage("full_bam/160_ova_23_29_genome.sorted.bam")
#tr.reads <- read.bam.coverage("transcripts_reads.sorted.bam")

#h3k9[which(h3k9$idx == 0)] <- 1
#pol2[which(pol2$idx == 0)] <- 1
#pi[which(pi$idx == 0)] <- 1

h3k9.bed <- read.table("/mnt/raid/illumina/ftp_data/Rezvykh.Chip/H3K9_1.sorted_peaks.bed",header = FALSE)
pol2.bed <- read.table("/mnt/raid/illumina/ftp_data/Rezvykh.Chip/Pol_1.sorted_peaks.bed", header = FALSE)

parse.bed <- function(x){
  bed.src <- x
  bed.src <- bed.src[which(bed.src$V1 == "scaffold_13050"),]
  bed.src <- bed.src[which(bed.src$V2 > g.start),]
  bed.src <- bed.src[which(bed.src$V2 < g.end),]
  bed.src$V2 <- bed.src$V2 - g.start
  bed.src$V3 <- bed.src$V3 - g.start
  return(bed.src)
}

h3k9.bed <- parse.bed(h3k9.bed)
pol2.bed <- parse.bed(pol2.bed)
df <- as.data.frame(df)

rg <- rtracklayer::import("~/genomes/12_genomes/Dvir/dvir-all-r1.06.gtf")
rg=as.data.frame(rg)

rg <- rg[which(rg$seqnames == "scaffold_13050"),]
rg <- rg[which(rg$start > g.start),]
rg <- rg[which(rg$start < g.end),]
rg$start <- rg$start - g.start
rg$end <- rg$end - g.start

rg$type_num <- revalue(rg$type, c("gene"="0", "mRNA"="1", "5UTR"="2", "exon"="3", "start_codon"="4", "CDS"="5", "stop_codon"="6","3UTR"="7"))

myb.3 <- 20001
myb.5 <- 22879
ran.3 <- 143218
ran.5 <- 192819


#pdf("PLOTZ.pdf")
ggplot(data = df) + 
  geom_segment(aes(x=qbegin, 
                   xend=qend,
                   y=abs(perc), 
                   yend=abs(perc),
                   color = class), 
               size = 0.1 ,
               arrow = arrow(length = unit(0.2, "cm"), 
                             ends = ifelse(df$sum.inv == "Straight", "first", "last"), type = "closed")) + 
  scale_shape_manual(values=1:nlevels(df$class)) + 
  theme_bw() + 
  geom_segment(aes(x = start, xend = end, y = -84, yend = -84),size=1.5, color = ifelse(df$class == "Satellite", "red", "black")) + 
  scale_x_continuous(breaks=c(0,50000,100000,150000, 200000),
                     labels=c("0kb", "50kb", "100kb", "150kb", "200kb"), 
                     name = "Genomic region") + 
  scale_y_continuous(breaks=c(-150, -137,-125, -113, -100,-84,-63,-50,0,50),
                     labels=c("Pol2", "Pol2-peaks", "H3K9me3","H3K9me3-peaks", "Transcripts","Transposones", "Reference", "-100%", "0%","100%"),
                     name = "") +
  geom_hline(yintercept = 0, size=0.5) + 
  geom_line(data = transcripts, aes(x = idx, y = xcov-100),
            color=ifelse(transcripts$xcov > 0, "red", "white"),size = 1) + 
  geom_line(data = h3k9, aes(x = idx,
                             y = ifelse(-is.infinite(log2(xcov)), -125, log2(xcov)-125)),
            size = 0.1, color = ifelse(h3k9$xcov>0, "blue", "white")) + 
  geom_line(data = pol2, aes(x = idx,
                             y = ifelse(-is.infinite(log2(xcov)), -150, log2(xcov)-150)),
            size = 0.1, color = ifelse(pol2$xcov>0, "red", "white")) + 
  #geom_line(data = pi, aes(x = idx,
  #                         y = ifelse(-is.infinite(log2(xcov)), -175, log2(xcov)-175)),
  #          size = 0.1, color = ifelse(pi$xcov>0, "green", "white")) 
  geom_segment(data = h3k9.bed, aes(x = V2, xend = V3, y = -113, yend = -113), color = "blue", size = 1) + 
  geom_curve(data = h3k9.bed, aes(x = V2, xend = V3, y = -113, yend = -113), color = "blue", size = 0.1, curvature = 0.5) + 
  geom_segment(data = pol2.bed, aes(x = V2, xend = V3, y = -137, yend = -137), color = "red", size = 1) +
  geom_curve(data = pol2.bed, aes(x = V2, xend = V3, y = -137, yend = -137), color = "red", size = 0.1, curvature = 1)
  #geom_segment(aes(x = myb.3, xend = myb.5, y = -63, yend = -63)) + 
 # geom_text(aes(x = (myb.5+myb.3)/2, y = -55), label = "Myb") + 
  #geom_segment(aes(x = ran.3, xend = ran.5, y = -63, yend = -63)) + 
#  geom_text(aes(x = (ran.5+ran.3)/2, y = -55), label = "Ranbp16") +
  #geom_line(data = rin, aes(x = idx, y = rin-62), color = ifelse(rin$rin>0, "red", "white"))
#  geom_segment(data=rg, aes(x = start, xend = end, y = type_num, yend = type_num))
                                       
                                       
                                       

for_clust <- data.frame(df$div, df$del, df$ins, df$perc, df$repeat_id, df$class)
names(for_clust) <- c("div", "del", "ins", "perc","repeat_id", "class")
for_clust$perc <- as.numeric(for_clust$perc)
p <- prcomp(for_clust[,1:4], center = TRUE, scale. = TRUE)

plot(p$x, col = for_clust$class)


ggplot() + geom_segment(data=rg, aes(x = start, xend = end, y = type_num, yend = type_num)) + ylim(c(0,100))




####standalone_graph





h3k9.start <- h3k9[(h3k9$idx %in% h3k9.bed$V2),]
h3k9.stop <- h3k9[(h3k9$idx %in% h3k9.bed$V3),]
h3k9.col.df <- data.frame(h3k9.start$idx, h3k9.stop$idx)

pol2.start <- pol2[(pol2$idx %in% pol2.bed$V2),]
pol2.stop <- pol2[(pol2$idx %in% pol2.bed$V3),]
pol2.col.df <- data.frame(pol2.start$idx, pol2.stop$idx)





transcripts.bed <- read.table("~/genomes/virilis_160/tr.bed", header = FALSE)
transcripts.bed
transcripts.bed <-transcripts.bed[which(transcripts.bed$V2 > 1740221),]
transcripts.bed <-transcripts.bed[which(transcripts.bed$V3 > 1953041),]

transcripts.bed$V2 <- transcripts.bed$V2 - 1740221
transcripts.bed$V3 <- transcripts.bed$V3 - 1740221
transcripts.bed

g1 <- ggplot() + 
  geom_segment(data = transcripts.bed, aes(x = V2, xend = V3, y = 5, yend = 5)) + 
  geom_line(data = h3k9, aes(x = idx,
                             y = ifelse(-is.infinite(log2(xcov)), 10, log2(xcov)+10)),
            size = 0.1, color = ifelse(inrange(x = h3k9$idx, lower = h3k9.start$idx, upper = h3k9.stop$idx), "red", "blue"),
            alpha = ifelse(inrange(x = h3k9$idx, lower = h3k9.start$idx, upper = h3k9.stop$idx), 1, 0.2)) + 
  geom_line(data = pol2, aes(x = idx,
                             y = ifelse(-is.infinite(log2(xcov)), 20, log2(xcov)+20)),
            size = 0.1, color = ifelse(inrange(x = pol2$idx, lower = pol2.start$idx, upper = pol2.stop$idx), "red", "blue"),
            alpha = ifelse(inrange(x = pol2$idx, lower = pol2.start$idx, upper = pol2.stop$idx), 1, 0.2)) + 
  scale_y_continuous(breaks=c(0,10,20),
                     labels=c("Transcripts", "H3K9me3", "Pol2"),
                     name = "") + 
  scale_x_continuous(breaks=c(0,50000,100000,150000, 200000),
                     labels=c("0kb", "50kb", "100kb", "150kb", "200kb"), 
                     name = "Genomic region") + 
  theme_bw()
  

g1
g2 <- ggplot(data = df) + 
geom_segment(aes(x=qbegin, 
                   xend=qend,
                   y=perc, 
                   yend=perc,
                   color = class), 
               size = 0.1,
               arrow = arrow(length = unit(0.2, "cm"), 
                             ends = ifelse(df$sum.inv == "Straight", "first", "last"), type = "closed")) + 
  scale_x_continuous(breaks=c(0,50000,100000,150000, 200000),
                     labels=c("0kb", "50kb", "100kb", "150kb", "200kb"), 
                     name = "Genomic region") + 
  theme_bw()

g2
grid.arrange(g1, g2, ncol= 2)
dev.off()


ggplot() + geom_smooth(data = h3k9, aes(x = idx,
                           y = ifelse(-is.infinite(log2(xcov)), 10, log2(xcov)+10)),
          size = 0.1, color = ifelse(inrange(x = h3k9$idx, lower = h3k9.start$idx, upper = h3k9.stop$idx), "red", "blue"),
          alpha = ifelse(inrange(x = h3k9$idx, lower = h3k9.start$idx, upper = h3k9.stop$idx), 1, 0.2), method = "gam")





ggplot(h3k9) + 
  geom_line(aes(y=impressions, x=, color=cvr)) +
  stat_smooth(aes(y=impressions, x=hour), method = lm, formula = y ~ poly(x, 10), se = FALSE)  
                       
plot(x = h3k9$idx, y = h3k9$xcov, pch = "", xlab = "Genomic region, bp", ylab = "Coverage")
lines(lowess(h3k9$xcov~h3k9$idx,f = 0.001), col = "blue")
lines(lowess(pol2$xcov~pol2$idx,f = 0.001), col = "red")



par(mfrow=c(3,1))
plot(x = h3k9$idx, 
     y = h3k9$xcov, 
     type = "l", 
     xlab = "Genomic region, bp", 
     ylab = "Coverage",
     main = "H3K9 coverage",
     xlim=c(10000,30000))
lines(lowess(h3k9$xcov~h3k9$idx,f = 0.00001), col = "blue")


plot(x = pol2$idx, 
     y = pol2$xcov, 
     type = "l", 
     xlab = "Genomic region, bp", 
     ylab = "Coverage",
     main = "Pol2 coverage",
     xlim=c(10000,30000))
lines(lowess(pol2$xcov~pol2$idx,f = 0.001), col = "red")


plot(x = transcripts$idx, 
     y = transcripts$xcov, 
     type = "l", 
     xlab = "Genomic region, bp", 
     ylab = "Coverage",
     main = "Transcripts coverage",
     xlim=c(10000,30000),
     ylim = c(0,30))


dev.off()
       
