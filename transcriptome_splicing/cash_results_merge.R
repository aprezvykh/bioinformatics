options(scipen=999)
library("Rcpp")
library("Vennerable")
library(ggplot2)
library(gridExtra)
library(digest)
library(org.Mm.eg.db)
library("edgeR")
library("limma")
library("topGO")
library(stringr)
library(reshape)
library(ggbiplot)
library(xlsx)
library(dplyr)
library(ggrepel)
library(rtracklayer)
library(rebus)
mycols <- colors()[c(8, 5, 30, 53, 118, 72)]

#gtf <- import("~/transcriptomes/ref/gff/Mus_musculus.GRCm38.90.gtf")
#gtf <- as.data.frame(gtf)
###commonize
setwd("~/transcriptomes/reads/intron_retention/fus/cash/")
tg1tg2 <- read.delim("tg2vstg1.alldiff.txt")
tg2tg3 <- read.delim("tg3vstg2.alldiff.txt")
wt1tg1 <- read.delim("tg1vswt1.alldiff.txt")
wt1wt3 <- read.delim("wt3vswt1.alldiff.txt")


tg1tg2$st <- c("tg1tg2")
tg2tg3$st <- c("tg2tg3")
wt1tg1$st <- c("wt1tg1")
wt1wt3$st <- c("wt1wt3")
sod <- read.delim("~/transcriptomes/reads/intron_retention/sod_moto/cash/tgvswt.alldiff.txt")
tdp <- read.delim("~/transcriptomes/reads/intron_retention/tdp43/cash/tgvswt.alldiff.txt")
sod$st <- c("sod")
tdp$st <- c("tdp")

reftg1tg2 <- read.xlsx("~/ALS.results/final results/Tg-1-Tg-2/Results edgeR.xlsx", sheetIndex = 3)
reftg2tg3 <- read.xlsx("~/ALS.results/final results/Tg-2-Tg-3/Results edgeR-1.xlsx", sheetIndex = 3)
refwt1tg1 <- read.xlsx("~/ALS.results/final results/Control-1-Tg-1/Results edgeR.xlsx", sheetIndex = 3)
refwt1wt3 <- read.xlsx("~/ALS.results/final results/Control-1-Control-3/Results edgeR.xlsx", sheetIndex = 3)


reftg1tg2$st <- c("tg1tg2")
reftg2tg3$st <- c("tg2tg3")
refwt1tg1$st <- c("wt1tg1")
refwt1wt3$st <- c("wt1wt3")

glia <- read.csv("~/ALS.results/final results/sort/glial genes.csv")
moto <- read.csv("~/ALS.results/final results/sort/motoneuron genes.csv")

names(glia) <- c("n", "gene")
names(moto) <- c("n", "gene")

glia$tissue <- c("glia")
moto$tissue <- c("moto")

glia$name <- mapIds(org.Mm.eg.db, 
                        keys=as.character(glia$gene), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")

moto$name <- mapIds(org.Mm.eg.db, 
                        keys=as.character(moto$gene), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")



comm <- bind_rows(tg1tg2, tg2tg3, wt1tg1, wt1wt3, sod, tdp)
exp <- bind_rows(reftg1tg2, reftg2tg3, refwt1tg1, refwt1wt3)

parse_coords <- function(df){
  s <- str_split_fixed(df$Location, ":", 2)
  chr <- s[,1]
  start <- str_split_fixed(s[,2],"-",2)[,1]
  stop <- str_split_fixed(s[,2],"-",2)[,2]
  df$chr <- chr
  df$start <- as.numeric(start)
  df$stop <- as.numeric(stop)
  df$div <- df$stop - df$start
  return(df)
}
cutoff_sign <- function(df){
  zz <- df[which(df$FDR<0.05),]
  zz <- zz[which(abs(zz$delta_PSI)>0.1),]
  return(zz)
}

tdp <- cutoff_sign(tdp)
sod <- cutoff_sign(sod)


comm <- cutoff_sign(comm)
comm <- parse_coords(comm)
comm$micro <- ifelse(comm$div < 27, "micro", "normal")
comm <- comm[,!grepl("^[tg]*[wt]", names(comm))]


tg1tg2 <- cutoff_sign(tg1tg2)
tg2tg3 <- cutoff_sign(tg2tg3)


write.csv(comm, "all_significant_events.csv")
write.csv(exp, "all_significant_diff.csv")


#pdf("Report.pdf", family = "Helvetica")
ggplot(data=comm) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("All events, p < 0.05")
ggplot(data=comm[comm$micro == "micro",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("Micro events, p < 0.05")
ggplot(data=comm[comm$micro == "normal",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("Non-micro events, p < 0.05")

c_sub <- comm[which(comm$FDR > 1.00e-30),]

ggplot(data=c_sub, aes(x = delta_PSI, y = -log10(FDR))) + 
  geom_point(aes(shape = st,color = SplicingType, size = micro)) + 
  geom_text(aes(label=ifelse(c_sub$FDR<0.000000001 & abs(c_sub$delta_PSI) > 0.3, AccID, "")),hjust=1, vjust=1,check_overlap = F) +
  theme_bw()

microexons.moto.genes <- intersect(comm[which(comm$micro == "micro" & comm$SplicingType == "Cassette"),]$AccID, moto$name)

g1 <- ggplot(data=comm[comm$micro == "micro" & comm$SplicingType == "Cassette",],aes(x = st, y = delta_PSI, fill = st)) + 
  geom_boxplot(show.legend = F) + 
  geom_point(size = 0.5, position=position_jitterdodge(),show.legend = F) + 
  theme_bw() + 
  ggtitle("Microexons only > 27 n.t")
g2 <- ggplot(data=comm,aes(x = st, y = delta_PSI, fill = st)) + 
  geom_boxplot(show.legend = F) + 
  geom_point(size = 0.5, position=position_jitterdodge(),show.legend = F) + 
  theme_bw() + 
  ggtitle("All splicing events, p < 0.05")
g3 <- ggplot(data=comm[comm$SplicingType == "Cassette",],aes(x = st, y = delta_PSI, fill = st)) + 
  geom_boxplot(show.legend = F) + 
  geom_point(size = 0.5, position=position_jitterdodge(),show.legend = F) + 
  theme_bw() + 
  ggtitle("Cassete exons")
g4 <- ggplot(data=comm[comm$SplicingType == "IR",],aes(x = st, y = delta_PSI, fill = st)) + 
  geom_boxplot(show.legend = F) + 
  geom_point(size = 0.5, position=position_jitterdodge(),show.legend = F) + 
  theme_bw() + 
  ggtitle("Retained introns")

grid.arrange(g1,g2,g3,g4)


####plot unique gene
sub <- comm[which(comm$AccID == "Ubc"),]
ggplot(data=sub) + 
  geom_point(aes(x = st, y = delta_PSI, shape = SplicingType, fill = FDR, alpha = div),size = 5) + 
  ylim(c(-1,1)) +
  theme_bw()

#GTF_PARSER
get_splice_regions <- function(gene){
  gene_name <- gene
  coords_common <- data.frame(data[grep(gene_name, comm$AccID),]$Location,
                              data[grep(gene_name, comm$AccID),]$SplicingType,
                              data[grep(gene_name, comm$AccID),]$div)
  names(coords_common) <- c('coords','type', 'div')
  gstart <- vector()
  gstop <- vector()
  coords <- vector()
  spl_type <- vector()
  for(i in 1:nrow(coords_common)){
    coords <- append(as.character(coords_common$coords[i]),coords)
    spl_type <- append(as.character(coords_common$type[i]),spl_type)
    #structure.feature.list <- c("exon", "cds", "start_codon", "stop_codon", "gene")
    chr <- strsplit(coords,":")[[1]][1]
    gstart <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][1], gstart)
    gstop <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][2], gstop)
  }
  s_lim <- min(gstart)
  e_lim <- max(gstop)
  sub_gtf <- gtf[which(gtf$seqnames == chr),]
  structure <- sub_gtf[grep(paste("\\b",gene_name, "\\b", sep = ""), sub_gtf$gene_name),]
  structure <- as.data.frame(structure)
  
  data.frame(gstart,gstop,spl_type)
  gstart <- as.numeric(gstart)
  gstop <- as.numeric(gstop)
  rx <- number_range(gstart, gstop)
  return(as.character(structure[grep(rx, structure$start),]$type))
}

plot_alt <- function(data,gene){
    gene_name <- gene
    coords_common <- data.frame(data[grep(gene_name, comm$AccID),]$Location,
                                data[grep(gene_name, comm$AccID),]$SplicingType,
                                data[grep(gene_name, comm$AccID),]$div)
    names(coords_common) <- c('coords','type', 'div')
    gstart <- vector()
    gstop <- vector()
    coords <- vector()
    spl_type <- vector()
    structure.feature.list <- c("exon")  
      for(i in 1:nrow(coords_common)){
            coords <- append(as.character(coords_common$coords[i]),coords)
            spl_type <- append(as.character(coords_common$type[i]),spl_type)

            chr <- strsplit(coords,":")[[1]][1]
            gstart <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][1], gstart)
            gstop <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][2], gstop)
        }
    s_lim <- min(gstart)
    e_lim <- max(gstop)
    sub_gtf <- gtf[which(gtf$seqnames == chr),]
    structure <- sub_gtf[grep(paste("\\b",gene_name, "\\b", sep = ""), sub_gtf$gene_name),]
    structure <- as.data.frame(structure)
    #structure <- structure[structure$type %in% structure.feature.list,]
    g <- ggplot(data=structure) + geom_segment(aes(y = type, yend = type, x = start, xend = end, col = type),size = 5)
        for(i in 1:length(gstart)){                     
            g <- g + geom_segment(aes(x = as.numeric(gstart[i]) ,xend = as.numeric(gstop[i]), y = -1, yend = -1), size = 5) + 
                    geom_vline(aes(xintercept = as.numeric(gstart[i])),color = mycols[i]) + 
                    geom_vline(aes(xintercept = as.numeric(gstop[i])),color = mycols[i])
        }
    
    g1 <- ggplot(data=structure) + geom_segment(aes(y = type, yend = type, x = start, xend = end, col = type),size = 5) + xlim(c(s_lim,e_lim))
    for(i in 1:length(gstart)){                     
      g1 <- g1 + geom_segment(aes(x = as.numeric(gstart[i]) ,xend = as.numeric(gstop[i]), y = -1, yend = -1), size = 5) + 
        geom_vline(aes(xintercept = as.numeric(gstart[i])),color = mycols[i]) + 
        geom_vline(aes(xintercept = as.numeric(gstop[i])),color = mycols[i])
    }
    
    g <- g + theme_bw()
    g1 <- g1 + theme_bw()
    ga <- grid.arrange(g,g1)
    return(ga)

}


plot_alt(comm,"Mpc1")

get_splice_regions("Gprasp1")
View(comm)






data <- comm
gene_name <- "Mpc1"
coords_common <- data.frame(data[grep(gene_name, comm$AccID),]$Location,
                            data[grep(gene_name, comm$AccID),]$SplicingType,
                            data[grep(gene_name, comm$AccID),]$div)
names(coords_common) <- c('coords','type', 'div')
gstart <- vector()
gstop <- vector()
coords <- vector()
spl_type <- vector()
structure.feature.list <- c("exon")  
for(i in 1:nrow(coords_common)){
  coords <- append(as.character(coords_common$coords[i]),coords)
  spl_type <- append(as.character(coords_common$type[i]),spl_type)
  
  chr <- strsplit(coords,":")[[1]][1]
  gstart <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][1], gstart)
  gstop <- append(strsplit(strsplit(coords, ":")[[1]][2],"-")[[1]][2], gstop)
}
s_lim <- min(gstart)
e_lim <- max(gstop)
sub_gtf <- gtf[which(gtf$seqnames == chr),]
structure <- sub_gtf[grep(paste("\\b",gene_name, "\\b", sep = ""), sub_gtf$gene_name),]
structure <- as.data.frame(structure)
#structure <- structure[structure$type %in% structure.feature.list,]

g <- ggplot(data=structure) + geom_segment(aes(y = type, yend = type, x = start, xend = end, col = type),size = 5)
g
for(i in 1:length(gstart)){                     
  g <- g + geom_segment(aes(x = as.numeric(gstart[i]) ,xend = as.numeric(gstop[i]), y = -1, yend = -1), size = 5) + 
    geom_vline(aes(xintercept = as.numeric(gstart[i])),color = mycols[i]) + 
    geom_vline(aes(xintercept = as.numeric(gstop[i])),color = mycols[i])
}

g1 <- ggplot(data=structure) + geom_segment(aes(y = type, yend = type, x = start, xend = end, col = type),size = 5) + xlim(c(s_lim,e_lim))
for(i in 1:length(gstart)){                     
  g1 <- g1 + geom_segment(aes(x = as.numeric(gstart[i]) ,xend = as.numeric(gstop[i]), y = -1, yend = -1), size = 5) + 
    geom_vline(aes(xintercept = as.numeric(gstart[i])),color = mycols[i]) + 
    geom_vline(aes(xintercept = as.numeric(gstop[i])),color = mycols[i])
}

g <- g + theme_bw()
g1 <- g1 + theme_bw()
ga <- grid.arrange(g,g1)
return(ga)











