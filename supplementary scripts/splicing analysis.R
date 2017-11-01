library("gplots")
library(edgeR)
library("org.Mm.eg.db")
library(stringr)
library(GO.db)
col.pan <- colorpanel(100, "blue", "white", "red")

isoforms <- read.delim("~/counts/ALS Mice/experimental/splicing/Run_2017-10-23_14-14-03.isoforms.quantification.tsv")
names(isoforms) <- c("tg3_1", "tg3_2", "tg3_3", "tg3_4", "tg3_5", "cnt1_1", "cnt1_2", "cnt1_3", "tg1_1", "tg1_2", "tg1_3", "tg1_4", "tg2_1", "tg2_2", "tg1_5", "cnt1_4", "cnt3_1", "cnt3_2", "cnt3_3", "cnt1_5", "tg2_3", "tg2_4", "nth_1", "nth_2", "nth_3", "cnt3_5")
genes <- read.delim("~/counts/ALS Mice/experimental/splicing/Run_2017-10-23_14-14-03.genes.quantification.tsv")


sampleCondition <- c("case", "case", "case", "case", "case", 
                     "control", "control", "control", "control", "control")

sampleNames <- colnames(sub)
sampleTable <- data.frame(sampleName = sampleNames, sampleCondition= sampleCondition)

a <- DGEList(counts=sub, group = sampleTable$sampleCondition) 
CountsTable <- as.data.frame(y$counts)
cpm <- cpm(a) 
keep <- rowSums(cpm > 0.5) >= ncol(sampleTable) 
a <- a[keep, , keep.lib.sizes=FALSE] 
a <- calcNormFactors(a, method = "TMM") 
design <- model.matrix(~sampleTable$sampleCondition) 
a <- estimateDisp(a,design) 
fit <- glmQLFit(a,design, robust = TRUE) 
qlf <- glmQLFTest(fit,coef=ncol(fit$design))
et_annot <- as.data.frame(qlf$table) 
et_annot_non_filtered <- as.data.frame(qlf$table)
top <- as.data.frame(topTags(qlf))
logCPM <- as.data.frame(cpm(a, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
et <- qlf

et_annot <- subset(et_annot, PValue < 0.05)
rownames(et_annot) <- sub("(.*?)_.*", "\\1", rownames(et_annot))

et_annot$gene <- mapIds(org.Mm.eg.db, 
                         keys=row.names(et_annot), 
                         column="ENSEMBL", 
                         keytype="ENSEMBLTRANS",
                         multiVals="first")

et_annot$GOID <-   mapIds(org.Mm.eg.db, 
                          keys=et_annot$gene, 
                          column="GO", 
                          keytype="ENSEMBL",
                          multiVals="first")

et_annot$term <- mapIds(GO.db, 
                        keys=et_annot$GOID, 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")

et_annot$term <- as.character(et_annot$term)
et_annot <- et_annot[complete.cases(et_annot),]
u <- as.data.frame(unique(et_annot$gene))
df <- data.frame()
diff.splice <- data.frame()

for (f in u$`unique(et_annot$gene)`){
    z <- grep(paste(f), et_annot$gene)
    df <- et_annot[z,]
    df$disp <- var(df$logFC, na.rm = TRUE)
    diff.splice <- rbind(df, diff.splice)
}

diff.splice <- diff.splice[complete.cases(diff.splice),]
diff.splice <- subset(diff.splice, disp > 1)

diff.splice$symbol <- mapIds(org.Mm.eg.db, 
                        keys=diff.splice$gene, 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")

diff.splice$name <- mapIds(org.Mm.eg.db, 
                             keys=diff.splice$gene, 
                             column="GENENAME", 
                             keytype="ENSEMBL",
                             multiVals="first")


write.xlsx(diff.splice, file = "DIFFSPLICETG1-TG3.xlsx")

#HEATMAP

cpm <- cpm(a, log = TRUE)
h <- grep("igg", rownames(cpm), ignore.case = TRUE)
hm <- scale(cpm[h,])
heatmap.2(hm, col = col.pan)
