library("rentrez")
library("DOSE")
library("clusterProfiler")
library("ReactomePA")
experimental <- as.data.frame(read.csv("~/bioinformatics/counts/ALS Mice/experimental.csv"))
microglia <- as.data.frame(read.csv("~/bioinformatics/counts/ALS Mice/microglia.csv"))
motoneurons <- as.data.frame(read.csv("~/bioinformatics/counts/ALS Mice/motoneurons.csv"))
rownames(experimental) <- experimental$X
rownames(microglia) <- microglia$X
int <- intersect(experimental$X, microglia$X)
sdiff <- function(a,b){
  setdiff(union(a,b), intersect(a,b))}
mt <- sdiff(int, experimental$X)
moto <- experimental[mt,]
moto <- moto[complete.cases(moto), ]
int2 <- intersect(motoneurons$X, moto$X)
final <- moto[int2,]



dfb <- as.character(final$entrez)
y <- enrichKEGG(gene=dfb, organism = "mouse", minGSSize=1)
y
barplot(y, showCategory=10,  font.size = 9)

