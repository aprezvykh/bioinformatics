library("rentrez")
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


entrez_dbs()
r_search <- entrez_search(db="medgen", term="Sox9")
ids <- r_search$ids
for (f in ids){
    summary <- entrez_summary(db="medgen", id=paste(f))
    summary$definition
}