library(AnnotationDbi)
library(org.Mm.eg.db)
library(Rcpp)
library(edgeR)
library(xlsx)

### PASTE 1 IF YOU WANT TO ANALYZE ALL SAMPLES. PASTE 0 IF YOU WANT TO
pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5

### Statistical analysis
directory <- '~/GitHub/counts/ALS Mice/experimental/'
setwd(directory)
sampleFiles <- grep('control_1',list.files(directory),value=TRUE)
sampleCondition <- c('1', '1', '1', '1', '1')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
keep <- rowSums(cpm(y) > cpm_cutoff) > 1
y <- y[keep, ] 
cpm <- as.data.frame(cpm(y))
cpm$avg <- rowSums(cpm)
### ANNOTATE

cpm$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(cpm), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
cpm$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(cpm), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

write.csv(cpm, file = "control_1.csv")

cpm <- cpm[order(cpm$avg, decreasing = TRUE),]
cpm <- cpm[seq(1:1000),]
exp <- read.xlsx("~/GitHub/counts/ALS Mice/experimental/results/fc1/tg_1-tg_3/Results edgeR.xlsx", sheetIndex = 2)

int <- intersect(rownames(cpm), exp$NA.)
df <- data.frame()
data <- data.frame()

for (f in int){
  a <- grep(paste(f), exp$NA.)
  df <- exp[a,]
  data <- rbind(df, data)
  
}

g <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg3/glia/deg_glia_13.csv")
m <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg3/moto/deg_moto_13.csv")
o <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg3/other/deg_out_13.csv")

gl_c <- intersect(data$NA., g$NA.)
mo_c <- intersect(data$NA., m$NA.)
ou_c <- intersect(data$NA., o$NA.)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


tg1 <- read.csv("~/GitHub/counts/ALS Mice/CD4/tg1_expression_profile.csv")
tg1 <- tg1[order(tg1$sum, decreasing = TRUE),]
tg1 <- tg1[seq(1:1000),]
intersect(tg1$X, data$NA.)
outersect(a, data$NA.)

