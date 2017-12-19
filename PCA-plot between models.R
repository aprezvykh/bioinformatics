setwd("~/counts/ALS Mice/experimental/results/all/")
library(xlsx)
library(ggbiplot)
library(org.Mm.eg.db)
library(gplots)
library(dplyr)
library(calibrate)
library(limma)
tg_13_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/csv/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/csv/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/csv/tg13-other.csv")

cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
sod <- read.xlsx("~/counts/SOD/results/Tg-1-Tg-3/Results edgeR.xlsx", sheetIndex = 3)
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
c$type <- c("Non specific")
x <- bind_rows(tg_13_glia, tg_13_moto, tg_13_others)
x$fc <- round(x$FDR, digits = 2)
b <- rbind(a, b)
b <- rbind(c, b)
## PCA BETWEEN CELL TYPES
b$symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(b), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

f <- b[,seq(1:24)]
#f <- scale(f)
p <- prcomp(f, scale. = TRUE, center = TRUE)
pdf(file = "PCA Plot between cell types in FUS - TRUE PCA - by cell types.pdf", height = 10, width = 10, family = "Helvetica")

ggbiplot(p,
         ellipse = TRUE,
         var.axes = FALSE,
         circle = TRUE,
         groups = b$type) + 
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

pdf(file = "PCA Plot between models glia + sod.pdf", height = 10, width = 10, family = "Helvetica")
ggbiplot(p, data = z,
         ellipse = TRUE,
         groups = x$type,
         var.axes = FALSE, alpha = 1) + 
         ggtitle("PCA plot between FUS and SOD") + 
         theme(plot.title = element_text(hjust = 0.5)) + 
         theme_bw() + 
         geom_abline(intercept = 0.5, slope = -2.4)

dev.off()
getwd()
ggplot(x) + geom_point(aes(x = PValue, y = FDR, color = type))

pdf(file = "Volcano between SOD and FUS.pdf", height = 10, width = 10, family = "Helvetica")

ggplot(x) + geom_point(aes(x = logFC, y = -log10(FDR), color = type), alpha = 1) + 
                 theme_bw() + 
                 scale_x_continuous(name = ("Log2FoldChange")) + 
                 scale_y_continuous(name = ("-log10(FDR)")) + 
                 geom_text(aes(x = logFC, y = -log10(FDR), label = ifelse(abs(logFC) > 4 & -log10(FDR) > 2.5, as.character(symbol), ''), color = type), hjust = 0, vjust = 0, size = 3)
dev.off()


pdf(file = "Volcano between SOD and FUS w tick labels.pdf", height = 10, width = 10, family = "Helvetica")
ggplot(x) + geom_point(aes(x = logFC, y = -log10(PValue), color = type)) + 
            theme_bw() + 
            scale_x_continuous(breaks = c(-6,-3,0,3,6), name = "Log2FoldChange") + 
            ggtitle("Volcano plot") + 
            theme(plot.title = element_text(hjust = 0.5))
dev.off()



x <- read.xlsx(file = "~/counts/SOD/results/compare.xlsx", sheetIndex = 1)
x$type <- ifelse(x$logCPM.FUS > x$logCPM.SOD, "FUS", "SOD")
ggplot(x) + geom_point(aes(x = logFC.SOD, y = logFC.FUS)) + 
            theme_bw() + 
            geom_text(aes(x = logFC.SOD, y = logFC.FUS, label = symbol),
            hjust = 0, vjust = 0, size = 3) + 
            scale_x_continuous(name = ("LogFC-SOD")) + 
            scale_y_continuous(name = ("LogFC-FUS")) + 
            geom_smooth(aes(x = logFC.SOD, y = logFC.FUS), method = "lm", se = FALSE)
