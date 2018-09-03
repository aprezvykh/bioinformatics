source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("topGO")
library(ggplot2)
library(gridExtra)
library(digest)
library(org.Mm.eg.db)
library("edgeR")
library("limma")
library("topGO")
df <- read.delim("~/Documents/intron_retention/cash/tg3vstg2.alldiff.txt")

###commonize
setwd("~/Documents/intron_retention/cash/")
tg1tg2 <- read.delim("tg2vstg1.alldiff.txt")
tg2tg3 <- read.delim("tg3vstg2.alldiff.txt")
wt1tg1 <- read.delim("tg1vswt1.alldiff.txt")
wt1wt3 <- read.delim("wt3vswt1.alldiff.txt")

tg1tg2$ENS <- mapIds(org.Mm.eg.db, 
                     keys=as.character(tg1tg2$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

tg2tg3$ENS <- mapIds(org.Mm.eg.db, 
                     keys=as.character(tg2tg3$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

wt1tg1$ENS <- mapIds(org.Mm.eg.db, 
                     keys=as.character(wt1tg1$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

wt1wt3$ENS <- mapIds(org.Mm.eg.db, 
                     keys=as.character(wt1wt3$AccID), 
                     column="ENTREZID", 
                     keytype="SYMBOL",
                     multiVals="first")

i <- Reduce(intersect,list(tg1tg2$ENS, tg2tg3$ENS, wt1tg1$ENS, wt1wt3$ENS))
which(tg1tg2$ENS %in% i)


tg1tg2 <- tg1tg2[(tg1tg2$ENS %in% i),]
tg2tg3 <- tg2tg3[(tg1tg2$ENS %in% i),]
wt1tg1 <- wt1tg1[(wt1tg1$ENS %in% i),]
wt1wt3 <- wt1wt3[(wt1wt3$ENS %in% i),]

tg1tg2 <- tg1tg2[order(tg1tg2$ENS),]
tg2tg3 <- tg2tg3[order(tg2tg3$ENS),]
wt1tg1 <- wt1tg1[order(wt1tg1$ENS),]
wt1wt3 <- wt1wt3[order(wt1wt3$ENS),]


i <- Reduce(intersect,list(tg1tg2$ENS, tg2tg3$ENS, wt1tg1$ENS, wt1wt3$ENS))

comm <- data.frame(wt1wt3$AccID, wt1wt3$delta_PSI, wt1wt3$FDR,
                   wt1tg1$delta_PSI, wt1tg1$FDR,
                   tg1tg2$delta_PSI, tg1tg2$FDR,
                   tg2tg3$delta_PSI, tg2tg3$FDR)




###analysis pt

df <- tg2tg3
df <- df[which(df$FDR < 0.05),]
df <- df[which(abs(df$delta_PSI) > 0.1 ),]

g1 <- ggplot(data=df) + geom_boxplot(aes(x = SplicingType, y = delta_PSI, fill = SplicingType)) + 
  theme_bw() + 
  ggtitle("dPSI for events")

g2 <- ggplot(data=df) + geom_point(aes(y = -log10(P.Value), x = delta_PSI, col = "red")) + 
  theme_bw() + 
  ggtitle("Volcano plot")
ts <- table(t(df$SplicingType))
ts <- as.data.frame(ts)

g3 <- ggplot(data=ts) + geom_bar(aes(x = Var1, y = Freq), stat = "identity") + 
  theme_bw() + 
  ggtitle("Splicing events frequency, FDR<0.05, dPSI > 0.1, Tg3 vs Tg2")

g1
g2
g3

big <- data.frame()
for(f in unique(df$SplicingType)){
  sub <- df[which(df$SplicingType ==f ),]
  gt <- data.frame(f, mean(sub$delta_PSI), mean(sub$FDR))
  big <- rbind(gt,big)
}
names(big) <- c("alt", "dpsi", "padj")
ggplot(data=big) + geom_bar(aes(x = alt, y = dpsi, fill = padj), stat = "identity") + theme_bw()

###GO
    sub <- df[which(df$SplicingType == f),]
    go <- goana(sub$ENS,universe = df$ENS,species = "Mm")
    topgo <- topGO(go,ontology = "BP")
    topgo$perc <- 100*(topgo$DE/topgo$N)
    topgo <- head(topgo,10)
    g <- ggplot(data=topgo) + geom_bar(aes(x = Term, y = perc, fill = P.DE), stat = "identity") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    g
#    ggsave(paste(f,".png",sep = ""),g)
