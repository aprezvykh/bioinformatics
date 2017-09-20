install.packages("statmod")
library("edgeR")
library("org.Dm.eg.db")
library("xlsx")
library("pheatmap")

### Experimental design
setwd("~/Fly memory project/experimental/K_vs_F24/")
files <- c("fly_K1.counts", "fly_K2.counts", "fly_F24.counts", "fly_F24_2.counts")
class <- c("cont", "cont", "case", "case")
label <- c("cont1", "cont2", "case1", "case2")
ExpTable <- cbind(files, class, label)
y <- readDGE(files = files, group = class, labels = label)


y <- DGEList(y, group=class)
keep <- rowSums(cpm(y) > 100) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

CountsTable <- as.data.frame(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
top <- topTags(et)
et_annot <- as.data.frame(et$table)

et_annot$symbol <- mapIds(org.Dm.eg.db, 
                        keys=row.names(et_annot), 
                        column="SYMBOL", 
                        keytype="FLYBASE",
                        multiVals="first")

et_annot$name <- mapIds(org.Dm.eg.db, 
                          keys=row.names(et_annot), 
                          column="GENENAME", 
                          keytype="FLYBASE",
                          multiVals="first")
et_annot$logFC <- et_annot$logFC*(-1)

et_annot <- as.data.frame(subset(et_annot, PValue < 0.05))

CountsTable$symbol <- mapIds(org.Dm.eg.db, 
                          keys=row.names(CountsTable), 
                          column="SYMBOL", 
                          keytype="FLYBASE",
                          multiVals="first")

CountsTable$name <- mapIds(org.Dm.eg.db, 
                        keys=row.names(CountsTable), 
                        column="GENENAME", 
                        keytype="FLYBASE",
                        multiVals="first")

select <- et_annot[order(et_annot$logFC),]
select <- et_annot[1:20,]
b <- rownames(select)
z <- cpm(y, log=TRUE, prior.count = 1)
pheatmap(scale(z, center = TRUE)[b,])
dev.off()

