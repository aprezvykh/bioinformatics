lym <- read.csv("~/counts/ALS Mice/immune cells compare/CD22/cd22.csv")
lym <- lym[complete.cases(lym),]
lym <- lym[order(lym$avg, decreasing = TRUE),]
lym <- lym[seq(1:1000),]


glia <- read.csv("~/counts/ALS Mice/immune cells compare/microglia.csv")
glia <- glia[complete.cases(glia),]
glia <- glia[order(glia$avg, decreasing = TRUE),]
glia <- glia[seq(1:1000),]


moto <- read.csv("~/counts/ALS Mice/immune cells compare/moto.csv")
moto <- moto[complete.cases(moto),]
moto <- moto[order(moto$avg, decreasing = TRUE),]
moto <- moto[seq(1:1000),]

glia.13 <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/glia/tg13-glia.csv")
moto.13 <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/tg13-moto.csv")
other.13 <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/tg13-other.csv")

in1 <- intersect(lym$X, moto$X)
in2 <- intersect(glia$X, moto$X)
i <- data.frame(intersect(in1, in2))
rownames(i) <- i$intersect.in1..in2.


i$symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(i), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
i$name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(i), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")
names(i) <- c("ens", "symbol", "name")


intersect(other.13$X, i$ens)


write.xlsx(i, "Housekeeping genes.xlsx", sheetName = "HKG")
