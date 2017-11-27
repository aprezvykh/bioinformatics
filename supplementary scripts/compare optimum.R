library(dplyr)
glia <- read.csv("~/counts/ALS Mice/new filtering/microglia_expression_profile.csv", header = TRUE)
moto <- read.csv("~/counts/ALS Mice/new filtering/moto_expression_profile.csv", header = TRUE)
tg2 <- read.xlsx("~/counts/ALS Mice/filtering fc=1 + FDR/DEG_with_FDR.xlsx", sheetIndex = 9)
tg3 <- read.xlsx("~/counts/ALS Mice/filtering fc=1 + FDR/DEG_with_FDR.xlsx", sheetIndex = 10)
comp <- data.frame(glia$X, glia$sum, moto$sum)
rownames(comp) <- comp$glia.X
comp$glia.X <- NULL
names(comp) <- c("Glia", "Moto")

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

aff.glia <- comp[which(comp$Glia/comp$Moto >= 5),]
aff.moto <- comp[which(comp$Moto/comp$Glia >= 5),]
aff.glia$X <- rownames(aff.glia)
aff.moto$X <- rownames(aff.moto)
com <- bind_rows(aff.glia, aff.moto)
oth <- outersect(com$X, rownames(comp))
aff.other <- comp[(rownames(comp) %in% oth),]


deg.12 <- read.xlsx("~/counts/ALS Mice/filtering fc=1 + FDR/Deg 1-2.xlsx", sheetIndex = 1)
rownames(deg.12) <- deg.12$NA.
deg.12$NA. <- NULL
deg.13 <- read.xlsx("~/counts/ALS Mice/filtering fc=1 + FDR/Deg 1-3.xlsx", sheetIndex = 1)
rownames(deg.13) <- deg.13$NA.
deg.13$NA. <- NULL

tg12.glia <- deg.12[(intersect(tg2$head, rownames(aff.glia))),]
tg12.moto <- deg.12[(intersect(tg2$head, rownames(aff.moto))),]
tg12.other <- deg.12[(intersect(tg2$head, rownames(aff.other))),]


tg13.glia <- deg.13[(intersect(tg3$head, rownames(aff.glia))),]
tg13.moto <- deg.13[(intersect(tg3$head, rownames(aff.moto))),]
tg13.other <- deg.13[(intersect(tg3$head, rownames(aff.other))),]


setwd("~/counts/ALS Mice/filtering fc=1 + FDR/")
write.csv(tg12.glia, file = "tg12-glia.csv")
write.csv(tg12.moto, file = "tg12-moto.csv")
write.csv(tg12.other, file = "tg12-other.csv")

write.csv(tg13.glia, file = "tg13-glia.csv")
write.csv(tg13.moto, file = "tg13-moto.csv")
write.csv(tg13.other, file = "tg13-other.csv")
