install.packages("statmod")
library("edgeR")
library("org.Dm.eg.db")

### Experimental design
setwd("~/Fly memory project/experimental/F_vs_F24(memory)/")
files <- c("fly_F1.counts", "fly_F2.counts", "fly_F24.counts", "fly_F24_2.counts")
class <- c("cont", "cont", "case", "case")
label <- c("cont1", "cont2", "case1", "case2")
s_table <- cbind(files, class, label)
data <- readDGE(files = files, group = class, labels = label)



y <- DGEList (counts = data, group = data$samples$group)

y$tags$Symbol <- mapIds(org.Dm.eg.db, rownames(y),
                         keytype = "FLYBASE", column = "SYMBOL")
y <- y[!is.na(y$tags$Symbol), ]
keep <- rowSums(cpm(y) > 0.5) >= 2
table(keep)

y <- calcNormFactors(y)

design <- model.matrix(~0+myclass)
colnames(design) <- levels(myclass)
y <- estimateDisp(y, design, robust=TRUE)

fit <- glmQLFit(y, design, robust=TRUE)

tr <- glmTreat(fit, lfc=log2(1.5))
topTags(tr)
go <- goana(tr, species="Dm")
topGO(go, n=15)
