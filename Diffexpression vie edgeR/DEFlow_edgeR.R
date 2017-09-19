install.packages("statmod")
library("edgeR")
library("org.Dm.eg.db")
library("xlsx")

### Experimental design
setwd("~/Fly memory project/experimental/F_vs_F24(memory)/")
files <- c("fly_F1.counts", "fly_F2.counts", "fly_F24.counts", "fly_F24_2.counts")
class <- c("cont", "cont", "case", "case")
label <- c("cont1", "cont2", "case1", "case2")
s_table <- cbind(files, class, label)
y <- readDGE(files = files, group = class, labels = label)
head(y$counts)
