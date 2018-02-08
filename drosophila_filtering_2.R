library(xlsx)
new.kf <- read.xlsx("~/counts/dmel_memory_filtering/results/K-F/Results edgeR.xlsx", sheetIndex = 3)
new.fm <- read.xlsx("~/counts/dmel_memory_filtering/results/F-mem/Results edgeR.xlsx", sheetIndex = 3)
new.km <- read.xlsx("~/counts/dmel_memory_filtering/results/K-mem/Results edgeR.xlsx", sheetIndex = 3)

intersect(new.km$NA., rownames(et_annot))

new.kf <- subset(new.kf, new.kf$PValue < 0.01)
new.fm <- subset(new.fm, new.fm$PValue < 0.01) 
new.km <- subset(new.km, new.km$PValue < 0.01) 



outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

#subsetting data // 

new.kf.up <- subset(new.kf, new.kf$logFC > 0)
new.kf.down <- subset(new.kf, new.kf$logFC < 0)
new.fm.up <- subset(new.fm, new.fm$logFC > 0)
new.fm.down <- subset(new.fm, new.fm$logFC < 0)
new.km.up <- subset(new.km, new.km$logFC > 0)
new.km.down <- subset(new.km, new.km$logFC < 0)

new.kf.fm.up.common <- intersect(new.kf.up$NA., new.fm.up$NA.)
new.kf.fm.down.common <- intersect(new.kf.down$NA., new.fm.down$NA.)

new.kf.km.up.common <- intersect(new.kf.up$NA., new.km.up$NA.)
new.kf.km.down.common <- intersect(new.kf.down$NA., new.km.down$NA.)

new.fm.km.up.common <- intersect(new.fm.up$NA., new.km.up$NA.)
new.fm.km.down.common <- intersect(new.fm.down$NA., new.km.down$NA.)



fc.new.kf.fm.up <- new.kf.up[(new.kf.up$NA. %in% new.kf.fm.up.common),]
fc.new.kf.fm.down <- new.kf.down[(new.kf.down$NA. %in% new.kf.fm.down.common),]


fc.new.kf.km.up <- new.kf.up[(new.kf.up$NA. %in% new.kf.km.up.common),]
fc.new.kf.km.down <- new.kf.down[(new.kf.down$NA. %in% new.kf.km.down.common),]


fc.new.fm.km.up <- new.fm.up[(new.fm.up$NA. %in% new.fm.km.up.common),]
fc.new.fm.km.down <- new.fm.down[(new.fm.down$NA. %in% new.fm.km.down.common),]




kf.fm <- intersect(new.kf$NA., new.fm$NA.)

kffm.kf <- new.kf[(new.kf$NA. %in% kf.fm),]
kffm.fm <- new.fm[(new.fm$NA. %in% kf.fm),]



common.genes <- data.frame(kffm.kf$NA., kffm.kf$logFC, kffm.kf$logCPM, kffm.kf$PValue,
                                        kffm.fm$logFC, kffm.fm$logCPM, kffm.fm$PValue,
                                        kffm.kf$symbol, kffm.kf$name, kffm.kf$GOID, 
                                        kffm.kf$term)
names(common.genes) <- c("flybase", 
                         "KF-logfc", "KF-logCPM", "KF-PValue",
                         "FM-logfc", "FM-logCPM", "FM-PValue", 
                         "symbol", "name", "goid", "term")


updown <- common.genes[which(common.genes$`KF-logfc` > 0 & common.genes$`FM-logfc` < 0),]
downup <- common.genes[which(common.genes$`KF-logfc` < 0 & common.genes$`FM-logfc` > 0),]
downdown <- common.genes[which(common.genes$`KF-logfc` < 0 & common.genes$`FM-logfc` < 0),]
upup <- common.genes[which(common.genes$`KF-logfc` > 0 & common.genes$`FM-logfc` > 0),]




write.xlsx(new.kf, "new_PVAL.xlsx", sheetName = "new_kf", append = TRUE)
write.xlsx(new.fm, "new_PVAL.xlsx", sheetName = "new_fm", append = TRUE)
write.xlsx(new.km, "new_PVAL.xlsx", sheetName = "new_km", append = TRUE)

