library(xlsx)
directory <- ("~/GitHub/counts/ALS Mice/experimental/results/mge/")
setwd(directory)
SampleFiles <- grep("mouse", list.files(directory), value = TRUE)

sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('control_early', 'control_early', 'control_early', 'control_early', 'control_early', 
                     'control_late', 'control_late', 'control_late', 'control_late', 'control_late', 
                     'tg_early', 'tg_early', 'tg_early', 'tg_early', 'tg_early', 
                     'tg_mid', 'tg_mid', 'tg_mid', 'tg_mid', 
                     'tg_late', 'tg_late', 'tg_late', 'tg_late', 'tg_late')

sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
final <- read.table("~/GitHub/counts/ALS Mice/experimental/results/mge/seqlength.txt", sep = ",")
tab <- data.frame()

for (f in SampleFiles){
  a <- read.delim(f)
  final <- cbind(a$X0, final)
}
ncol(final)
names(final) <- c('control_early', 'control_early', 'control_early', 'control_early', 'control_early', 
                  'control_late', 'control_late', 'control_late', 'control_late', 'control_late', 
                  'tg_early', 'tg_early', 'tg_early', 'tg_early', 'tg_early', 
                  'tg_mid', 'tg_mid', 'tg_mid', 'tg_mid', 
                  'tg_late', 'tg_late', 'tg_late', 'tg_late', 'tg_late', 'name', 'length')

rownames(final) <- final$name
final$name <- NULL
newfinal <- final
cpm <- read.xlsx(file = "CPMtrans.xlsx", sheetIndex = 2)
rownames(cpm) <- cpm$NA.
cpm$NA. <- NULL

#cpm$avg <- rowSums(cpm)/ncol(cpm)

keep <- rowSums(cpm > 0.5) >= ncol(sampleTable)
cpm <- cpm[keep,]

cont <- grep("tg_early", colnames(cpm_all))
cont <- cpm_all[,cont]
cont$avg <- rowSums(cont)/nrow(cont)
case <- grep("tg_late", colnames(cpm_all))
case <- cpm_all[,case]
case$avg <- rowSums(case)/nrow(case)

df <- data.frame(rownames(cont), cont$avg, case$avg)


df$logfc <- log2(df$case.avg/df$cont.avg)
df <- subset(df, df$logfc > 0.5 | df$logfc < -0.5)


