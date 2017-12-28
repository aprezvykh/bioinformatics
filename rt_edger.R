library(edgeR)
library(xlsx)

dir <- c("~/counts/ALS Mice/experimental/results/suppl/mge/")
a <- grep("new", list.files(dir), value = TRUE)
setwd(dir)
length <- read.delim(file = "seqlength.txt", header = FALSE, sep = ",")
gr_control <- c("control_3")
gr_case <- c("tg_3")

files_control <- grep(paste(gr_control),list.files(dir),value=TRUE)
files_case <- grep(paste(gr_case),list.files(dir),value=TRUE)
sampleFiles <- c(files_control, files_case)
cond_control <- rep(paste(gr_control), length(files_control))
cond_case <- rep(paste(gr_case), length(files_case))
sampleCondition <- c(cond_control, cond_case)
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)        
a <- DGEList(counts=y, group = sampleTable$condition) 



CountsTable <- as.data.frame(y$counts)
cpm <- cpm(y) 
cpm <- as.data.frame(cpm(y)) 
cpm$rowsum <- rowSums(cpm) 
keep <- rowSums(cpm > 0.5) >= ncol(sampleTable) 
logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
logCPM <- logCPM[keep,]
a <- a[keep,] 
a <- calcNormFactors(a, method = "TMM") 
design <- model.matrix(~sampleTable$condition) 
a <- estimateDisp(a,design) 
fit <- glmQLFit(a,design, robust = TRUE) 
qlf <- glmQLFTest(fit,coef=ncol(fit$design))
et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
top <- as.data.frame(topTags(qlf, n = 20))



et_annot <- as.data.frame(subset(et_annot, logCPM > 0))
et_annot <- as.data.frame(subset(et_annot, PValue < 0.05))
et_annot <- as.data.frame(subset(et_annot, FDR < 0.05))
et_annot <- as.data.frame(subset(et_annot, logFC > 0.5 | logFC < -0.5))
et_annot <- et_annot[complete.cases(et_annot), ]
l <- length[(length$V1 %in% rownames(et_annot)),]
et_annot$length <- l$V2

write.xlsx(et_annot, file = "RTs edgeR.xlsx", append = TRUE, sheetName = paste(gr_control, "_", gr_case))


