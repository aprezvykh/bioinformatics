library("pheatmap")
library(devtools)
library(org.Mm.eg.db)
library(ggbiplot)
library(RColorBrewer)
library(gplots)
library("GenomicRanges")
library(SAJR)
library("genefilter")

ncol.pan <- colorpanel(100, "blue", "white", "red")
setwd("~/Documents/intron_retention/")
data <- readRDS("data.rds")
data.f <- readRDS("data_filtered.rds")
all.alts <- readRDS("alts.rds")

mod = list(f=factor(c('tg1', 'tg1', 'tg1', 'tg1', 'tg1', 
                      'tg2', 'tg2', 'tg2', 'tg2', 
                      'tg3', 'tg3', 'tg3', 'tg3', 'tg3', 
                      'wt1', 'wt1', 'wt1', 'wt1', 'wt1', 
                      'wt3', 'wt3', 'wt3', 'wt3', 'wt3')))

n <- c('tg1', 'tg1', 'tg1', 'tg1', 'tg1', 
       'tg2', 'tg2', 'tg2', 'tg2', 
       'tg3', 'tg3', 'tg3', 'tg3', 'tg3', 
       'wt1', 'wt1', 'wt1', 'wt1', 'wt1', 
       'wt3', 'wt3', 'wt3', 'wt3', 'wt3')



data.f.glm = fitSAGLM(data.f,terms(x~mod$f), mod, .parallel = TRUE)


data.f.pv = calcSAPvalue(data.f.glm)
data.f.pv[,2] = p.adjust(data.f.pv[,2],method='BH')
data.sign = data.f[data.f.pv[,2] <=0.05,]
sign.res <- data.sign[order(abs(data.sign$ir[,1]-data.sign$ir[,2]),decreasing = TRUE),]

##basic statistics things

data_for_analysis <- data
data <- NULL
data_for_analysis$ir[is.na(data_for_analysis$ir)] <- 0

p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$ir)>0,]$ir), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Inclusion rate, PCA")

p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$i)>0,]$i), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Introns, PCA")

p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$e)>0,]$e), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Exons, PCA")


pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$ir)>0,]$ir, method = "pearson"),color = col.pan,main = "Inclusion ratio, Spearman corr.")
pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$i)>0,]$i, method = "pearson"),color = col.pan,main = "Introns, Spearman corr.")
pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$e)>0,]$e, method = "pearson"),color = col.pan,main = "Exons, Spearman corr.")

df <- NULL
df <- data.frame(as.character(colMeans(data_for_analysis$e)))
df <- as.numeric(as.character(df$mean))
df$sds <- as.character(apply(data_for_analysis$ir, 2, sd))
df$gr <- n
names(df) <- c("mean", "SD", "group")
df$num <- seq(1:length(n))
ggplot(data=df) + geom_bar(aes(x = df$num, y = df$mean, fill = df$group), stat = "identity") 


