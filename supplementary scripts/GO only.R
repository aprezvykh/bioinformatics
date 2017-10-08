library(topGO)
et_annot <- read.xlsx(file = "motoneurons marker.xlsx", sheetName = "moto", append = TRUE)

data(go.sets.mm)
data(go.subs.mm)

foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE,set.size = c(10, 100), rank.test = TRUE)
lapply(gomfres, head)
gomfres <- as.data.frame(gomfres)
gomfres <- gomfres[complete.cases(gomfres), ]
write.xlsx(gomfres, file = "GO.xlsx", sheetName = "GO_MF", append = TRUE)

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- as.data.frame(gobpres)
gobpres <- gobpres[complete.cases(gobpres), ]
write.xlsx(gobpres, file = "GO.xlsx", sheetName = "GO_BP", append = TRUE)

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- as.data.frame(goccres)
goccres <- goccres[complete.cases(goccres), ]
write.xlsx(goccres, file = "GO.xlsx", sheetName = "GO_CC", append = TRUE)

library(DOSE)
