setwd("~/counts/ALS Mice/experimental/results/all/")
library(xlsx)
library(ggbiplot)
library(org.Mm.eg.db)
library(gplots)
library(dplyr)
tg_13_glia <- read.csv("~/counts/ALS Mice/filtering + FDR/csv/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering + FDR/csv/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering + FDR/csv/tg13-other.csv")

cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
sod <- read.xlsx("~/counts/ALS Mice/SOD/Tg-1-Tg-3/Results edgeR.xlsx", sheetIndex = 3)
tdp <- read.xlsx("~/counts/ALS Mice/TDP43/results/control-tg/Results edgeR.xlsx", sheetIndex = 3)
mot <- read.xlsx("~/counts/ALS Mice/SOD moto/results/Control-Tg/Results edgeR.xlsx", sheetIndex = 3)

tg_13_glia$type <- c("FUS-Microglia")
tg_13_moto$type <- c("FUS-Motoneurons")
tg_13_others$type <- c("FUS-No Specific")
sod$type <- c("SOD-Microglia")
tdp$type <- c("TDP")
mot$type <- c("SOD-Motoneurons")

names(sod) <- c("X", "logFC", "logCPM", "F", "PValue", "FDR", "symbol",
                "name", "entrez", "GOID", "term", "type")

names(tdp) <- c("X", "logFC", "logCPM", "F", "PValue", "FDR", "symbol",
                "name", "entrez", "GOID", "term", "type")

names(mot) <- c("X", "logFC", "logCPM", "F", "PValue", "FDR", "symbol",
                "name", "entrez", "GOID", "term", "type")

b <- NULL
a <- cpm[(cpm$X %in% tg_13_glia$X),]
rownames(a) <- a$X
a$X <- NULL
a$type <- c("Microglia")

b <- cpm[(cpm$X %in% tg_13_moto$X),]
rownames(b) <- b$X
b$X <- NULL
b$type <- c("Motoneurons")

c <- cpm[(cpm$X %in% tg_13_others$X),]
rownames(c) <- c$X
c$X <- NULL
c$type <- c("Non specific)")

b <- rbind(a, b)
b <- rbind(c, b)

## PCA BETWEEN CELL TYPES
b$symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(b), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

f <- b[,seq(1:24)]

p <- prcomp(f, scale. = TRUE, center = TRUE)
pdf(file = "PCA Plot between cell types in FUS.pdf", height = 10, width = 10, family = "Helvetica")

ggbiplot(p,
         ellipse = TRUE,
         var.axes = FALSE,
         circle = TRUE,
         groups = b$type
        ) +
         theme_bw()
dev.off()


### ANALYSIS BETWEEN MODELS
nrow(tg_13_moto)

x <- data.frame(stringsAsFactors = FALSE)
x <- rbind(tg_13_glia, x)
x <- rbind(tg_13_moto, x)
x <- rbind(tg_13_others, x)
x <- rbind(sod, x)
x$entrez <- NULL

x$lfc <- ifelse(x$logFC > 0, paste("Positive"), paste("Negative"))
x$expr <- ifelse(x$logCPM > 1, paste("High"), paste("Low"))
x$index <- seq(1:nrow(x))
z <- x[,2:6]
p <- prcomp(z, scale. = TRUE, center = TRUE)

pdf(file = "PCA Plot between models NEW.pdf", height = 10, width = 10, family = "Helvetica")
ggbiplot(p, data = z,
         ellipse = TRUE,
         groups = x$type,
         var.axes = FALSE, alpha = 1) + 
         ggtitle("PCA plot between FUS and SOD") + 
         theme(plot.title = element_text(hjust = 0.5)) + 
         theme_bw()

dev.off()
ggplot(x) + geom_point(aes(x = PValue, y = FDR, color = type))

fus <- bind_rows(tg_13_glia, tg_13_moto, tg_13_others)

i <- intersect(fus$X, sod$X)

c_fus <- fus[(fus$X %in% i),]
c_sod <- sod[(sod$X %in% i),]

s <- data.frame(c_fus$logFC, c_sod$logFC, c_fus$type, c_fus$logCPM, c_sod$logCPM)
c <- apply(s[,4:5], 1, mean)
s <- cbind(c, s)
names(s) <- c( "LogCPM-value", "fus", "sod", "type", "logCPM-FUS", "logCPM-SOD")
pdf(file = "LogFC scatter w CPM.pdf", height = 10, width = 10, family = "Helvetica")
ggplot(s) + geom_point(aes(x = fus, y = sod, color = type, size = `LogCPM-value`), alpha = 0.7) + 
                 theme_bw() + 
                 xlim(-7,7) + 
                 ylim(-7,7) + 
                 geom_vline(xintercept = 0) + 
                 geom_hline(yintercept = 0) + 
                 scale_x_continuous(name = ("Log2FoldChange in delta-FUS")) + 
                 scale_y_continuous(name = ("Log2FoldChange in SOD"))

dev.off()
