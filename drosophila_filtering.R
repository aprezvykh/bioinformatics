library(xlsx)
old.kf <- read.xlsx("~/counts/dmel_memory_old/results/K-F/Results edgeR.xlsx", sheetIndex = 3)
old.fm <- read.xlsx("~/counts/dmel_memory_old/results/F-mem/Results edgeR.xlsx", sheetIndex = 3)
old.km <- read.xlsx("~/counts/dmel_memory_old/results/K-mem/Results edgeR.xlsx", sheetIndex = 3)


new.kf <- read.xlsx("~/counts/dmel_memory_new/results/K-F/Results edgeR.xlsx", sheetIndex = 3)
new.fm <- read.xlsx("~/counts/dmel_memory_new/results/F-mem/Results edgeR.xlsx", sheetIndex = 3)
new.km <- read.xlsx("~/counts/dmel_memory_new/results/K-mem/Results edgeR.xlsx", sheetIndex = 3)

old.kf <- subset(old.kf, old.kf$PValue < 0.01)
old.fm <- subset(old.fm, old.fm$PValue < 0.01) 
old.km <- subset(old.km, old.km$PValue < 0.01) 


new.kf <- subset(new.kf, new.kf$PValue < 0.01)
new.fm <- subset(new.fm, new.fm$PValue < 0.01) 
new.km <- subset(new.km, new.km$PValue < 0.01) 

intersect(rownames(et_annot), new.kf$NA.)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

#subsetting data // 
old.kf.up <- subset(old.kf, old.kf$logFC > 0)
old.kf.down <- subset(old.kf, old.kf$logFC < 0)
old.fm.up <- subset(old.fm, old.fm$logFC > 0)
old.fm.down <- subset(old.fm, old.fm$logFC < 0)
old.km.up <- subset(old.km, old.km$logFC > 0)
old.km.down <- subset(old.km, old.km$logFC < 0)

new.kf.up <- subset(new.kf, new.kf$logFC > 0)
new.kf.down <- subset(new.kf, new.kf$logFC < 0)
new.fm.up <- subset(new.fm, new.fm$logFC > 0)
new.fm.down <- subset(new.fm, new.fm$logFC < 0)
new.km.up <- subset(new.km, new.km$logFC > 0)
new.km.down <- subset(new.km, new.km$logFC < 0)

#new.old common
new.old.kf.up.common <- intersect(old.kf.up$NA., new.kf.up$NA.)
new.old.kf.down.common <- intersect(old.kf.down$NA., new.kf.down$NA.)

new.old.fm.up.common <- intersect(old.fm.up$NA., new.fm.up$NA.)
new.old.fm.down.common <- intersect(old.fm.down$NA., new.fm.down$NA.)

new.old.km.up.common <- intersect(old.km.up$NA., new.km.up$NA.)
new.old.km.down.common <- intersect(old.km.down$NA., new.km.down$NA.)

fc.new.old.kf.up.common <- new.kf[(new.kf.up$NA. %in% new.old.kf.up.common),]
fc.new.old.kf.down.common <- new.kf[(new.kf.down$NA. %in% new.old.kf.down.common),]



#new
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


#old
old.kf.fm.up.common <- intersect(old.kf.up$NA., old.fm.up$NA.)
old.kf.fm.down.common <- intersect(old.kf.down$NA., old.fm.down$NA.)

old.kf.km.up.common <- intersect(old.kf.up$NA., old.km.up$NA.)
old.kf.km.down.common <- intersect(old.kf.down$NA., old.km.down$NA.)

old.fm.km.up.common <- intersect(old.fm.up$NA., old.km.up$NA.)
old.fm.km.down.common <- intersect(old.fm.down$NA., old.km.down$NA.)



fc.old.kf.fm.up <- old.kf.up[(old.kf.up$NA. %in% old.kf.fm.up.common),]
fc.old.kf.fm.down <- old.kf.down[(old.kf.down$NA. %in% old.kf.fm.down.common),]


fc.old.kf.km.up <- old.kf.up[(old.kf.up$NA. %in% old.kf.km.up.common),]
fc.old.kf.km.down <- old.kf.down[(old.kf.down$NA. %in% old.kf.km.down.common),]


fc.old.fm.km.up <- old.fm.up[(old.fm.up$NA. %in% old.fm.km.up.common),]
fc.old.fm.km.down <- old.fm.down[(old.fm.down$NA. %in% old.fm.km.down.common),]



##########write
write.xlsx(old.kf, "all_drosophila_oldandnew.xlsx", sheetName = "old_kf", append = TRUE)
write.xlsx(old.fm, "all_drosophila_oldandnew.xlsx", sheetName = "old_fm", append = TRUE)
write.xlsx(old.km, "all_drosophila_oldandnew.xlsx", sheetName = "old_km", append = TRUE)
write.xlsx(new.kf, "all_drosophila_oldandnew.xlsx", sheetName = "new_kf", append = TRUE)
write.xlsx(new.fm, "all_drosophila_oldandnew.xlsx", sheetName = "new_fm", append = TRUE)
write.xlsx(new.km, "all_drosophila_oldandnew.xlsx", sheetName = "new_km", append = TRUE)



