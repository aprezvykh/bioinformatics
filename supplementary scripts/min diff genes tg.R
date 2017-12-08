sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                     'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                     'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                     'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')    


row.names.remove <- c("NA.1")
cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
cpm <- cpm[complete.cases(cpm),]
rownames(cpm) <- cpm$X
cpm$X <- NULL
colnames(cpm) <- sampleCondition


#a <- grep("Tg", colnames(cpm))
tg <- cpm
tg$var <- rowSdDiffs(tg)

tg <- tg[order(tg$var, decreasing = FALSE),]
tg <- tg[seq(1:500),]
tg$rowsum <- rowSums(tg)
tg$symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(tg), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

tg$name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(tg), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")
tg <- tg[order(tg$rowsum, decreasing = TRUE),]

write.xlsx(tg, "Minimal difference genes.xlsx", sheetName = "MDG")
x <- read.xlsx("~/Housekeeping genes.xlsx", sheetIndex = 1)


hkg <- as.data.frame(intersect(x$ens, rownames(tg)))
rownames(hkg) <- hkg$`intersect(x$ens, rownames(tg))`

hkg$symbol <- mapIds(org.Mm.eg.db, 
                    keys=row.names(hkg), 
                    column="SYMBOL", 
                    keytype="ENSEMBL",
                    multiVals="first")

hkg$name <- mapIds(org.Mm.eg.db, 
                  keys=row.names(hkg), 
                  column="GENENAME", 
                  keytype="ENSEMBL",
                  multiVals="first")

write.xlsx(hkg, "for rtPCR.xlsx", sheetName = "PCR")
