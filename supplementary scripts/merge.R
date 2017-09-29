### Script finds top genes from two data sets, that are overexpressed, but heve minimum differential expressed
### avg is BaseMean/gene counts. It can be made with DESeq2 from DEFlow script in this git
install.packages("compare")
library(gtools)
library(xlsx)
library(compare)
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


kf <- read.csv("~/Fly memory project/merge/K_vs_F.csv")
kf24 <- read.csv("~/Fly memory project/merge/K_vs_F24.csv")
mlow <- read.xlsx("~/Fly memory project/merge/m_downreg.xlsx", sheetIndex = 1)
mup <- read.xlsx("~/Fly memory project/merge/m_upreg.xlsx", sheetIndex = 1)
m <- smartbind(mlow, mup)

commonf24 <- intersect(kf24$X, m$Gene)
commonf<- intersect(kf$X, m$Gene)
dff <- data.frame()
dff24 <- data.frame()
for (f in commonf24){
  u <- match(f, kf24$X)
  newrow <-data.frame(kf24[u,])
  dff24 <- rbind(dff24, newrow)
  
}
for (f in commonf){
  u <- match(f, kf$X)
  newrow <-data.frame(kf[u,])
  dff <- rbind(dff, newrow)
  
}
dff <- data.frame(dff$X, dff$log2FoldChange, dff$symbol)
dff24 <- data.frame(dff24$X, dff24$log2FoldChange, dff24$symbol)

mdff24 <- data.frame()
mdff <- data.frame()
for (f in commonf24){
  u <- match(f, m$Gene)
  newrow <-data.frame(m[u,])
  mdff24 <- rbind(mdff24, newrow)
  
}
for (f in commonf){
  u <- match(f, m$Gene)
  newrow <-data.frame(m[u,])
  mdff <- rbind(mdff, newrow)
  
}


mdff <- data.frame(mdff$Gene, mdff$Fold.Difference, mdff$Name)
mdff24 <- data.frame(mdff24$Gene, mdff24$Fold.Difference, mdff24$Name)

write.xlsx(dff, file = "dff.xlsx", sheetName = "dff", append = TRUE)
write.xlsx(dff24, file = "dff24.xlsx", sheetName = "dff", append = TRUE)
write.xlsx(mdff, file = "mdff.xlsx", sheetName = "dff", append = TRUE)
write.xlsx(mdff24, file = "mdff24.xlsx", sheetName = "dff", append = TRUE)