library("gage")
library("gageData")
df <- read.table(file = "C:/Users/rezvykh/Dropbox/ALS-mice project/counts_trimmed_geo/deg_bm_100_fc_1.5/csv/2_late_tg_vs_ctrl_tg.csv", sep = ",", header= TRUE)
df <- df[,colSums(is.na(df))<nrow(df)]
foldchanges = df$log2FoldChange
names(foldchanges) = df$entrez
### GO ###
data(go.sets.mm)
data(go.subs.mm)
gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE)
lapply(gomfres, head)
write.csv(gomfres, file = "6_MF.csv")


gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
write.csv(gobpres, file = "6_BP.csv")

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
write.csv(goccres, file = "6_CC.csv")


### (foldchanges)
require(clusterProfiler)
require(reactome.db)

df <- as.character(df$entrez)
print(df)
summary(df)
class(df)
x = enrichPathway(gene=df, organism = "mouse", minGSSize=10, readable = TRUE )
write.csv(x, file = "GSEA.csv")
head(as.data.frame(x))
dotplot(x, showCategory=12)
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7 )
cnetplot(x, categorySize="pvalue", foldChange=foldchanges)

cnetplot()

