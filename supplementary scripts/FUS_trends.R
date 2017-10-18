library(edgeR)
library(xlsx)
directory <- '~/GitHub/counts/ALS Mice/experimental/FUS_trends/'
setwd(directory)
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('1', '2', '3', '4', '5', '6', '7', '8')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)

cpm <- as.data.frame(cpm(y))
fits <- data.frame()
for (f in seq(1:nrow(cpm))){
        a <- cpm[f,]
        a <- t(a)
        a <- as.data.frame(a)
        gene <- colnames(a)
        a$num <- paste(seq(1:8))
        names(a) <- c("gene", "num")
        a$gene <- as.numeric(a$gene)
        a$num <- as.numeric(a$num)
        fit <- lm(a)
        df <- data.frame(fit$coefficients)
        d <- data.frame(gene, df[1,1], df[2,1])
        fits <- rbind(d, fits)
}
names(fits) <- c("Gene", "Intercept", "num")
backup <- fits
 
fits <- as.data.frame(subset(fits, Intercept > 0 | num > 0))
fits$angle <- NULL

exp <- read.xlsx("~/GitHub/counts/ALS Mice/experimental/results/tg_2-tg_3/Results edgeR.xlsx", sheetIndex = 2)

top_fits <- as.data.frame(subset(fits, num > 10))
intersect(top_fits$Gene, exp$NA.)
abline(1, 8)
