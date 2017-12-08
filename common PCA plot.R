library(DESeq2)
library(beepr)
library(devtools)
library(ggbiplot)
directory <- "~/counts/for PCAplot/"
setwd(directory)
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3',
                     'Moto-SOD', 'Moto-SOD',
                     'Control-TDP', 'Control-TDP', 'Control-TDP', 'Control-TDP', 
                     'Tg-TDP', 'Tg-TDP', 'Tg-TDP', 'Tg-TDP', 
                     'Tg-1-FUS', 'Tg-1-FUS', 'Tg-1-FUS', 'Tg-1-FUS', 'Tg-1-FUS',
                     'Tg-1-SOD', 'Tg-1-SOD', 'Tg-1-SOD', 
                     'Tg-2-FUS', 'Tg-2-FUS', 'Tg-2-FUS', 'Tg-2-FUS', 
                     'Tg-2-SOD', 'Tg-2-SOD', 'Tg-2-SOD', 
                     'Tg-3-FUS', 'Tg-3-FUS', 'Tg-3-FUS', 'Tg-3-FUS', 'Tg-3-FUS',
                     'Tg-3-SOD', 'Tg-3-SOD', 'Tg-3-SOD')    


sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
col <- as.vector(sampleTable$sampleName)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res <- results(dds, tidy = FALSE )
rld<- rlogTransformation(dds, blind=TRUE)
pdf(file = "PCA.pdf", width = 12, height = 17, family = "Helvetica")
print(plotPCA(rld, intgroup=c('condition')))
dev.off()
beep()



y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
p <- prcomp(y$counts, center = TRUE, scale. = TRUE)
biplot(p)
