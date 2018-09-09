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

glia <- read.csv("~/counts/ALS Mice/glial genes.csv")
moto <- read.csv("~/counts/ALS Mice/motoneuron genes.csv")

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

comm <- cutoff_sign(comm)
comm <- parse_coords(comm)
comm$micro <- ifelse(comm$div < 27, "micro", "normal")
comm <- comm[,!grepl("^[tg]*[wt]", names(comm))]


tg1tg2 <- cutoff_sign(tg1tg2)
tg2tg3 <- cutoff_sign(tg2tg3)


pdf("Report.pdf", family = "Helvetica")
ggplot(data=comm) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("All events, p < 0.05")
ggplot(data=comm[comm$micro == "micro",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("Micro events, p < 0.05")
ggplot(data=comm[comm$micro == "normal",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = SplicingType)) + theme_bw() + ggtitle("Non-micro events, p < 0.05")

comm$FDR
c_sub <- comm[which(comm$FDR > 1.00e-30),]

ggplot(data=c_sub, aes(x = delta_PSI, y = -log10(FDR))) + 
  geom_point(aes(shape = st,color = SplicingType, size = micro)) + 
  geom_text(aes(label=ifelse(c_sub$FDR<0.000000001 & abs(c_sub$delta_PSI) > 0.3, AccID, "")),hjust=1, vjust=1,check_overlap = F) +
  theme_bw()

microexons.moto.genes <- intersect(comm[which(comm$micro == "micro" & comm$SplicingType == "Cassette"),]$AccID, moto$name)
microexons.moto.genes


g1 <- ggplot(data=comm[comm$micro == "micro" & comm$SplicingType == "Cassette",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + theme_bw() + ggtitle("Microexons only > 27 n.t")
g2 <- ggplot(data=comm) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + theme_bw() + ggtitle("All splicing events, p < 0.05")
g3 <- ggplot(data=comm[comm$SplicingType == "Cassette",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + theme_bw() + ggtitle("Cassete exons")
g4 <- ggplot(data=comm[comm$SplicingType == "IR",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + theme_bw() + ggtitle("Retained introns")

grid.arrange(g1,g2,g3,g4)


tdp <- cutoff_sign(tdp)
sod <- cutoff_sign(sod)

intersect(tdp$AccID, comm$AccID[which(comm$st == "tg2tg3")])

ggplot(data=exp) + geom_point(aes( x = PValue, y = FDR, col = st))
ggplot(data=comm) + geom_point(aes( x = P.Value, y = FDR, col = st))

Reduce(intersect, list(tdp$AccID,sod$AccID,tg1tg2$AccID))




####plot unique gene
sub <- comm[which(comm$AccID == "Bhlhb9"),]
ggplot(data=sub) + geom_point(aes(x = st, y = delta_PSI, shape = SplicingType))
