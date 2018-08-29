install.packages("gridExtra")
library("pheatmap")
library(devtools)
library(org.Mm.eg.db)
library(ggbiplot)
library(RColorBrewer)
library(gplots)
library("GenomicRanges")
library(SAJR)
library("genefilter")
library(reshape2)
library("gridExtra")

col.pan <- colorpanel(100, "blue", "white", "red")
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

plot_boxplots <- function(data, group_vec){
    df <- data.frame(as.character(colMeans(data$ir)),stringsAsFactors = F)
    names(df) <- c("mean")
    df$mean <- as.numeric(as.character(df$mean))
    df$sds <- as.character(apply(data$ir, 2, sd))
    df$sds <- as.numeric(as.character(df$sds))
    df$gr <- group_vec
    names(df) <- c("mean", "SD", "group")
    for_m <- data.frame(df$mean,df$group)
    names(for_m) <- c("mean","group")
    g1 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) +
                               ggtitle("Inclusion ratio by groups") + 
                               theme_bw()

    
    df <- data.frame(as.character(colMeans(data$i)),stringsAsFactors = F)
    names(df) <- c("mean")
    df$mean <- as.numeric(as.character(df$mean))
    df$sds <- as.character(apply(data$i, 2, sd))
    df$sds <- as.numeric(as.character(df$sds))
    df$gr <- group_vec
    names(df) <- c("mean", "SD", "group")
    for_m <- data.frame(df$mean,df$group)
    names(for_m) <- c("mean","group")
    g2 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) + 
                              ggtitle("Introns") + theme_bw()
    
    
    df <- data.frame(as.character(colMeans(data$e)),stringsAsFactors = F)
    names(df) <- c("mean")
    df$mean <- as.numeric(as.character(df$mean))
    df$sds <- as.character(apply(data$e, 2, sd))
    df$sds <- as.numeric(as.character(df$sds))
    df$gr <- group_vec
    names(df) <- c("mean", "SD", "group")
    for_m <- data.frame(df$mean,df$group)
    names(for_m) <- c("mean","group")
    g3 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) + 
                               ggtitle("Exons") + theme_bw()
    grid.arrange(g1,g2,g3)
}
plot_boxplots(data_for_analysis, n)
all_vars_pca <- function(data,group_vec){
    big.df <- data.frame()
    df <- data$ir
    big.df <- rbind(df, big.df)
    df <- data$i
    big.df <- rbind(df, big.df)
    df <- data$e
    big.df <- rbind(df, big.df)
    p <- prcomp(t(big.df), center = TRUE)
    ggbiplot(p, var.axes = F, groups = group_vec, ellipse = T) + theme_bw() + ggtitle("All vars PCA")
}
all_vars_pca(data_for_analysis,n)

df <- data.frame(colMeans(data_for_analysis$i), colMeans(data_for_analysis$e))
df$gr <- n
names(df) <- c("i", "e", "groups")
ggplot(data=df) + geom_point(aes(x = i, y = e, fill = groups, color = groups)) + theme_bw()

