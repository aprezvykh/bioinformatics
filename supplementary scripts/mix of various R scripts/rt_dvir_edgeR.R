library(edgeR)
library(xlsx)
dir <- "~/counts/D.virilis.RT/"
setwd(dir)
g <- grep("txt", list.files(dir), value = TRUE)
df <- data.frame()
for (f in g){
  z <- read.delim(file = paste(f), header = FALSE, sep = "")
  p <- as.data.frame(z$V3)
  names(p) <- paste(f)
  p <- t(p)
  df <- rbind(p, df)
}

df <- as.data.frame(t(df))

rownames(df) <- z$V1

write.xlsx(cpm, file = "101_cpm.xlsx", sheetName = "trans")


sampleCondition <- c("N", "N", "M", "M")
sampleTable <- data.frame(sampleName = colnames(df), condition = sampleCondition)
a <- DGEList(counts=df, group = sampleCondition)
a <- calcNormFactors(a, method = "TMM") 
design <- model.matrix(~sampleTable$condition) 
a <- estimateDisp(a,design) 
fit <- glmQLFit(a,design = design, robust = TRUE) 
qlf <- glmQLFTest(fit,coef=ncol(fit$design))
et_annot <- as.data.frame(topTags(qlf, n = 222, adjust.method = "BH"))

et_annot <- subset(et_annot, et_annot$FDR < 0.05)
write.xlsx(et_annot, file = "101_RT.xlsx", sheetName = "trans")
cpm <- as.data.frame(cpm(a))
