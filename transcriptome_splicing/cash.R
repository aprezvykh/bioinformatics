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
###commonize
setwd("~/SAJR/cash/")
tg1tg2 <- read.delim("tg2vstg1.alldiff.txt")
tg2tg3 <- read.delim("tg3vstg2.alldiff.txt")
wt1tg1 <- read.delim("tg1vswt1.alldiff.txt")
wt1wt3 <- read.delim("wt3vswt1.alldiff.txt")

reftg1tg2 <- read.xlsx("~/counts/ALS Mice/experimental/results/final results/Tg-1-Tg-2/Results edgeR.xlsx", sheetIndex = 3)
reftg2tg3 <- read.xlsx("~/counts/ALS Mice/experimental/results/final results/Tg-2-Tg-3/Results edgeR-1.xlsx", sheetIndex = 3)
refwt1tg1 <- read.xlsx("~/counts/ALS Mice/experimental/results/final results/Control-1-Tg-1/Results edgeR.xlsx", sheetIndex = 3)
refwt1wt3 <- read.xlsx("~/counts/ALS Mice/experimental/results/final results/Control-1-Control-3/Results edgeR.xlsx", sheetIndex = 3)


tg1tg2$Location <- as.character(tg1tg2$Location)
tg2tg3$Location <- as.character(tg2tg3$Location)
wt1tg1$Location <- as.character(wt1tg1$Location)
wt1wt3$Location <- as.character(wt1wt3$Location)


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


tg1tg2 <- parse_coords(tg1tg2) 
tg2tg3 <- parse_coords(tg2tg3) 
wt1tg1 <- parse_coords(wt1tg1)
wt1wt3 <- parse_coords(wt1wt3)

tg1tg2$AccID <- as.character(tg1tg2$AccID)
tg2tg3$AccID <- as.character(tg2tg3$AccID)
wt1tg1$AccID <- as.character(wt1tg1$AccID)
wt1wt3$AccID <- as.character(wt1wt3$AccID)


cutoff_sign <- function(df){
    zz <- df[which(df$FDR<0.05),]
    zz <- zz[which(abs(zz$delta_PSI)>0.1),]
  return(zz)
}

tg1tg2 <- cutoff_sign(tg1tg2)
tg2tg3 <- cutoff_sign(tg2tg3)
wt1tg1 <- cutoff_sign(wt1tg1)
wt1wt3 <- cutoff_sign(wt1wt3)


ii <- intersect(tg1tg2$AccID, tg2tg3$AccID)
iii <- intersect(wt1tg1$AccID, wt1wt3$AccID)

i <- intersect(ii,iii)


tg1tg2$entrez <- mapIds(org.Mm.eg.db, 
                     keys=as.character(tg1tg2$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

tg2tg3$entrez <- mapIds(org.Mm.eg.db, 
                     keys=as.character(tg2tg3$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

wt1tg1$entrez <- mapIds(org.Mm.eg.db, 
                     keys=as.character(wt1tg1$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

wt1wt3$entrez <- mapIds(org.Mm.eg.db, 
                     keys=as.character(wt1wt3$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")





tg1tg2_c <- tg1tg2[which(tg1tg2$AccID %in% ii),]
tg2tg3_c <- tg1tg2[which(tg2tg3$AccID %in% ii),]

tg1tg2_c <- tg1tg2_c[order(tg1tg2_c$AccID),]
tg2tg3_c <- tg2tg3_c[order(tg2tg3_c$AccID),]



###analysis pt

df <- tg2tg3

df$micro <- ifelse(df$div<28,"micro","normal")

g1 <- ggplot(data=df) + geom_boxplot(aes(x = SplicingType, y = delta_PSI, fill = SplicingType)) + 
  theme_bw() + 
  ggtitle("dPSI for events")

g2 <- ggplot(data=df) + geom_point(aes(y = -log10(P.Value), x = delta_PSI, col = SplicingType, size = micro)) + 
  theme_bw() + 
  ggtitle("Volcano plot")


ts <- table(t(df$SplicingType))
ts <- as.data.frame(ts)
pie(x = ts$Freq, labels = ts$Var1)

g3 <- ggplot(data=ts) + geom_bar(aes(x = Var1, y = Freq), stat = "identity") + 
  theme_bw() + 
  ggtitle("Splicing events frequency, FDR<0.05, dPSI > 0.1")




dev.off()
g1
g2
g3


micro <- df[df$micro == "micro",]

###microexons

go <- goana(sub$entrez,species = "Mm")

topgo <- topGO(go,ontology = "BP")
topgo$perc <- 100*(topgo$DE/topgo$N)
topgo <- head(topgo,10)
ggplot(data=topgo) + geom_bar(aes(x = Term, y = perc, fill = P.DE), stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top-10 enriched GO terms, BP")

topgo <- topGO(go,ontology = "MF")
topgo$perc <- 100*(topgo$DE/topgo$N)
topgo <- head(topgo,10)
ggplot(data=topgo) + geom_bar(aes(x = Term, y = perc, fill = P.DE), stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top-10 enriched GO terms, MF")

topgo <- topGO(go,ontology = "CC")
topgo$perc <- 100*(topgo$DE/topgo$N)
topgo <- head(topgo,10)
ggplot(data=topgo) + geom_bar(aes(x = Term, y = perc, fill = P.DE), stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top-10 enriched GO terms, CC")

kk <- kegga(sub$entrez, species.KEGG = "mmu")
tk <- topKEGG(kk,number = 10)
ggplot(data=tk) + geom_bar(aes(x = Pathway, y = DE, fill = P.DE), stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Top-10 enriched KEGG-pathways")



###intersecting
i1 <- intersect(tg1tg2$AccID, reftg1tg2$symbol)
i2 <- intersect(tg2tg3$AccID, reftg2tg3$symbol)
i3 <- intersect(wt1tg1$AccID, refwt1tg1$symbol)
i4 <- intersect(wt1wt3$AccID, refwt1wt3$symbol)


c_spl <- data.frame(tg2tg3[tg2tg3$AccID %in% i2,])
c_diff <- data.frame(reftg2tg3[reftg2tg3$symbol %in% i2,])


c_spl <- c_spl[order(c_spl$AccID),]
c_diff <- c_diff[order(c_diff$symbol),]


table(t(as.character(c_spl$S)))

venn.plot <- venn.diagram(list(i1,i2,i3,i4), NULL, fill=c("red", "green", "yellow", "blue"), cex = 2, cat.fontface=4, category.names=c("Tg2 vs Tg1", "Tg3 vs Tg2", "Tg1 vs Wt1", "Wt3 vs Wt1"))
grid.draw(venn.plot)
dev.off()
