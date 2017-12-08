library(xlsx)
library(ggbiplot)
library(org.Mm.eg.db)
tg_12_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/glia/tg12-glia.csv")
tg_12_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/moto/tg12-moto.csv")
tg_12_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/other/tg12-other.csv")

tg_13_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/glia/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/tg13-other.csv")

cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
sod <- read.xlsx("~/counts/SOD/results/Tg-1-Tg-3/Results edgeR.xlsx", sheetIndex = 3)
sod_early <- read.xlsx("~/counts/SOD/results/Tg-1-Tg-2/Results edgeR.xlsx", sheetIndex = 3)
tdp <- read.xlsx("~/counts/ALS Mice/TDP43/results/control-tg/Results edgeR.xlsx", sheetIndex = 3)
nrow(tdp)

names(sod) <- c("X", "logFC", "logCPM", "F", "PValue", "FDR", "symbol",
                "name", "entrez", "GOID", "term")

names(tdp) <- c("X", "logFC", "logCPM", "F", "PValue", "FDR", "symbol",
                "name", "entrez", "GOID", "term")

b <- NULL
a <- cpm[(cpm$X %in% tg_13_glia$X),]
rownames(a) <- a$X
a$X <- NULL
a$type <- c("FUS-model(Microglia)")

b <- cpm[(cpm$X %in% tg_13_moto$X),]
rownames(b) <- b$X
b$X <- NULL
b$type <- c("FUS-model(Motoneurons)")

c <- cpm[(cpm$X %in% tg_13_others$X),]
rownames(c) <- c$X
c$X <- NULL
c$type <- c("FUS-model(Non specific)")

d <- cpm[(cpm$X %in% sod$X),]
rownames(d) <- d$X
d$X <- NULL
d$type <- c("SOD-model")
e <- cpm[(cpm$X %in% tdp$X),]
rownames(e) <- e$X
e$X <- NULL
e$type <- c("TDP-model")

b <- rbind(a, b)
b <- rbind(c, b)
b <- rbind(d, b)
b <- rbind(e, b)

x <- data.frame()
x <- rbind(tg_13_glia, x)
x <- rbind(tg_13_moto, x)
x <- rbind(tg_13_others, x)
x <- rbind(sod, x)
x <- rbind(tdp, x)
x$entrez <- NULL

b$symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(b), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

f <- b[,seq(1:24)]

p <- prcomp(f, scale. = TRUE, center = TRUE)
#pdf(file = "PCA Models Compare with tissues and gene names.pdf", height = 10, width = 10, family = "Helvetica")



ggbiplot(p, ellipse = TRUE, groups = b$type, var.axes = FALSE, circle = TRUE) + theme_bw() +
        ggtitle("PCA analysis of ALS mice models") + 
        theme(plot.title = element_text(hjust = 0.5))


dev.off()



plot(x$logFC, x$ )
