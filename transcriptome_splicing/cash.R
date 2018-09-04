library(ggplot2)
library(gridExtra)
library(digest)
library(org.Mm.eg.db)
library("edgeR")
library("limma")
library("topGO")

###commonize
setwd("~/transcriptomes/reads/intron_retention/fus/cash/")
tg1tg2 <- read.delim("tg2vstg1.alldiff.txt")
tg2tg3 <- read.delim("tg3vstg2.alldiff.txt")
wt1tg1 <- read.delim("tg1vswt1.alldiff.txt")
wt1wt3 <- read.delim("wt3vswt1.alldiff.txt")


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


sub <- df[which(df$SplicingType == f),]
go <- goana(sub$entrez,universe = df$entrez,species = "Mm")
topgo <- topGO(go,ontology = "BP")
topgo$perc <- 100*(topgo$DE/topgo$N)
topgo <- head(topgo,10)
g <- ggplot(data=topgo) + geom_bar(aes(x = Term, y = perc, fill = P.DE), stat = "identity") + 
theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
#    ggsave(paste(f,".png",sep = ""),g)

kk <- kegga(sub$entrez, species.KEGG = "mmu")
tk <- topKEGG(kk)
