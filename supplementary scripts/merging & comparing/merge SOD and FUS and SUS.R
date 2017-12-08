tg_12_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/glia/tg12-glia.csv")
tg_12_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/moto/tg12-moto.csv")
tg_12_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/other/tg12-other.csv")

tg_13_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/glia/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/tg13-other.csv")


glia_12 <- read.xlsx("~/counts/SOD/results/Tg-1-Tg-2/Results edgeR.xlsx", sheetIndex = 3)
glia_23 <- read.xlsx("~/counts/SOD/results/Tg-2-Tg-3/Results edgeR.xlsx", sheetIndex = 3)
glia_13 <- read.xlsx("~/counts/SOD/results/Tg-1-Tg-3/Results edgeR.xlsx", sheetIndex = 3)

moto <- read.xlsx("")

in_glia <- as.data.frame(intersect(tg_13_glia$X, glia_13$NA.))
in_moto <- as.data.frame(intersect(tg_13_moto$X, glia_13$NA.))
in_others <- as.data.frame(intersect(tg_13_others$X, glia_13$NA.))

length(intersect(tg_13_glia$X, glia_13$NA.))/nrow(tg_13_glia)*100
length(intersect(tg_13_moto$X, glia_13$NA.))/nrow(tg_13_glia)*100
length(intersect(tg_13_others$X, glia_13$NA.))/nrow(tg_13_glia)*100

intersect(glia_12$NA., tg_13_glia$X)
