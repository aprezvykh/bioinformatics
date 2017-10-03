library(edgeR)
library(org.Mm.eg.db)
library(xlsx)

### Regression coefficients!
up_reg <- 1.5
down_reg <- 0.5
up_reg_n <- -0.5
down_reg_n <- -1.5

directory <- '~/counts_ens/'
setwd(directory)
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control_early', 'control_early', 'control_early', 'control_early', 'control_early', 
                     'control_late', 'control_late', 'control_late', 'control_late', 'control_late', 
                     'tg_early', 'tg_early', 'tg_early', 'tg_early', 'tg_early', 
                     'tg_mid', 'tg_mid', 'tg_mid', 'tg_mid', 
                     'tg_late', 'tg_late', 'tg_late', 'tg_late', 'tg_late')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleFiles, group = sampleCondition, labels = sampleFiles)
cpm <- cpm(y)
readqual <- as.data.frame(tail(y$counts, 3))

y <- calcNormFactors(y, method = "TMM")
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature")
cpm <- cpm(y)
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ]
cpm <- as.data.frame(cpm(y))


cpm$ntg1 <- (cpm$mouse_ntg_1_176.counts + 
               cpm$mouse_ntg_1_177.counts + 
               cpm$mouse_ntg_1_178.counts +
               cpm$mouse_ntg_1_213.counts + 
               cpm$mouse_ntg_1_225.counts)/5

cpm$ntg3 <- (cpm$mouse_ntg_3_222.counts +
               cpm$mouse_ntg_3_223.counts +
               cpm$mouse_ntg_3_224.counts +
               cpm$mouse_ntg_3_230.counts +
               cpm$mouse_ntg_3_245.counts)/5

cpm$tg1 <- (cpm$mouse_tg_1_182.counts +
              cpm$mouse_tg_1_183.counts +
              cpm$mouse_tg_1_184.counts +
              cpm$mouse_tg_1_185.counts +
              cpm$mouse_tg_1_212.counts)/5
cpm$tg2 <- (cpm$mouse_tg_2_204.counts +
              cpm$mouse_tg_2_205.counts +
              cpm$mouse_tg_2_226.counts +
              cpm$mouse_tg_2_227.counts)/4
cpm$tg3 <- (cpm$mouse_tg_3_171.counts +
              cpm$mouse_tg_3_172.counts +
              cpm$mouse_tg_3_173.counts +
              cpm$mouse_tg_3_174.counts +
              cpm$mouse_tg_3_175.counts)/5

### INCREASING TRENDS
avgcpm <- data.frame(cpm$ntg1, cpm$ntg3, cpm$tg1, cpm$tg2, cpm$tg3)
rownames(avgcpm) <- rownames(cpm)

trends <- data.frame()
for (n in seq(1:nrow(avgcpm))){
  row <- avgcpm[n,]
  if (row[,5] > row[,4] & row[,4] > row[,3] & row[,1] >= row[,2]) { 
    trends <- rbind(trends, row)
  } else {
    print("Trend not exists!")
  }
}

coords <- data.frame(trends$cpm.tg1, trends$cpm.tg2, trends$cpm.tg3)
rownames(coords) <- rownames(trends)
regress <- data.frame()
for (f in seq(1:nrow(coords))){
  x <- as.data.frame(coords[f,])
  v1 <- c(1, 2)
  v2 <- c(x$trends.cpm.tg1, x$trends.cpm.tg2)
  v3 <- c(2, 3)
  v4 <- c(x$trends.cpm.tg2, x$trends.cpm.tg3)
  fit1 <- lm(v1~v2)
  co1 <- as.data.frame(coef(fit1))
  co1 <- co1[2,]
  fit2 <- lm(v3~v4)
  co2 <- as.data.frame(coef(fit2))
  co2 <- co2[2,]
  coff <- data.frame(co1, co2)
  rownames(coff) <- rownames(x)
  regress <- rbind(coff, regress)
}


subset1 <- as.data.frame(subset(regress, co1 > down_reg))
subset2 <- as.data.frame(subset(subset1, co1 < up_reg))
subset3 <- as.data.frame(subset(subset2, co2 > down_reg))
subset4 <- as.data.frame(subset(subset3, co2 < up_reg))
subset_incr <- as.data.frame(subset4)

trgenes_up <- avgcpm[rownames(subset_incr),]

### DECREASING TRENDS
avgcpm <- data.frame(cpm$ntg1, cpm$ntg3, cpm$tg1, cpm$tg2, cpm$tg3)
rownames(avgcpm) <- rownames(cpm)

trends <- data.frame()
for (n in seq(1:nrow(avgcpm))){
  row <- avgcpm[n,]
  if (row[,5] < row[,4] & row[,4] < row[,3] & row[,2] >= row[,1]) { 
    trends <- rbind(trends, row)
  } else {
    print("Trend not exists!")
  }
}

coords <- data.frame(trends$cpm.tg1, trends$cpm.tg2, trends$cpm.tg3)
rownames(coords) <- rownames(trends)
regress <- data.frame()
for (f in seq(1:nrow(coords))){
  x <- as.data.frame(coords[f,])
  v1 <- c(1, 2)
  v2 <- c(x$trends.cpm.tg1, x$trends.cpm.tg2)
  v3 <- c(2, 3)
  v4 <- c(x$trends.cpm.tg2, x$trends.cpm.tg3)
  fit1 <- lm(v1~v2)
  co1 <- as.data.frame(coef(fit1))
  co1 <- co1[2,]
  fit2 <- lm(v3~v4)
  co2 <- as.data.frame(coef(fit2))
  co2 <- co2[2,]
  coff <- data.frame(co1, co2)
  rownames(coff) <- rownames(x)
  regress <- rbind(coff, regress)
}


subset1 <- as.data.frame(subset(regress, co1 > down_reg_n))
subset2 <- as.data.frame(subset(subset1, co1 < up_reg_n))
subset3 <- as.data.frame(subset(subset2, co2 > down_reg_n))
subset4 <- as.data.frame(subset(subset3, co2 < up_reg_n))

subset_decr <- as.data.frame(subset4)
trgenes_down <- avgcpm[rownames(subset_decr),]

trgenes_down <- trgenes_down[complete.cases(trgenes_down), ]
trgenes_up <- trgenes_up[complete.cases(trgenes_up), ]

### ANNOT
subset_incr$symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(subset_incr), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

subset_incr$name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(subset_incr), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")

subset_decr$symbol <- mapIds(org.Mm.eg.db, 
                             keys=row.names(subset_decr), 
                             column="SYMBOL", 
                             keytype="ENSEMBL",
                             multiVals="first")

subset_decr$name <- mapIds(org.Mm.eg.db, 
                           keys=row.names(subset_incr), 
                           column="GENENAME", 
                           keytype="ENSEMBL",
                           multiVals="first")

write.xlsx(subset_incr, file = "Trends.xlsx", sheetName = "Trends up", append = TRUE)
write.xlsx(subset_decr, file = "Trends.xlsx", sheetName = "Trends down", append = TRUE)
