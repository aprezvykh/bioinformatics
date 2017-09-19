### Script finds top genes from two data sets, that are overexpressed, but heve minimum differential expressed
### avg is BaseMean/gene counts. It can be made with DESeq2 from DEFlow script in this git

moto.avg <- 609.685
exp.avg <- 637.596
  
experimental <- as.data.frame(read.csv("~/motoneurons compare/experimental.csv"))
motoneurons <- as.data.frame(read.csv("~/motoneurons compare/motoneurons.csv"))

common <- intersect(experimental$X, motoneurons$X)
common <- as.list(common)

dfdiff <- function(mer_exp, mer_moto) {
mer_exp <- data.frame()
mer_moto <- data.frame()
for (i in common){
  m <- match(i, experimental$X)
  print(m)
  mer_exp <- rbind(experimental[m,], mer_exp)
}

for (i in common){
  m <- match(i, motoneurons$X)
  print(m)
  mer_moto <- rbind(motoneurons[m,], mer_moto)
}

mer_bm <- data.frame(mer_exp$X, mer_exp$symbol, mer_exp$name, mer_exp$entrez, mer_exp$baseMean, mer_moto$baseMean)

mer_bm$mer_exp.baseMean <- as.numeric(as.character(mer_bm$mer_exp.baseMean))/exp.avg
mer_bm$mer_moto.baseMean <- as.numeric(as.character(mer_bm$mer_moto.baseMean))/moto.avg
mer_bm$Delta <- mer_bm$mer_exp.baseMean/mer_bm$mer_moto.baseMean

sub <- as.data.frame(subset(mer_bm, Delta > 0.9))
sub_low <- as.data.frame(subset(sub, Delta < 1.1))

return(sub_low)
}

dataset <- dfdiff(mer_exp, mer_moto) 
