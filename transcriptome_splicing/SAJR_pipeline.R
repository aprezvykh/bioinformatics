library(pheatmap)
library(devtools)
library(org.Mm.eg.db)
library(ggbiplot)
library(RColorBrewer)
library(gplots)
library(GenomicRanges)
library(SAJR)
library(genefilter)
library(reshape2)
library(gridExtra)
library(plotrix)

col.pan <- colorpanel(100, "blue", "white", "red")
setwd("~/SAJR/")
data <- readRDS("data.rds")
data.f <- readRDS("data_filtered.rds")
all.alts <- readRDS("alts.rds")
gdata <- readRDS("gdata.rds")

###COMMON CORALLATION AND STUFF
data_for_analysis <- data
data <- NULL
data_for_analysis$ir[is.na(data_for_analysis$ir)] <- 0


p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$ir)>0,]$ir), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Inclusion rate, PCA")

p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$i)>0,]$i), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Inclusion reads, PCA")

p <- prcomp(t(data_for_analysis[rowSums(data_for_analysis$e)>0,]$e), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Exclusion reads, PCA")


pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$ir)>0,]$ir, method = "pearson"),color = col.pan,main = "Inclusion ratio, Spearman corr.")
pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$i)>0,]$i, method = "pearson"),color = col.pan,main = "Inclusion reads, Spearman corr.")
pheatmap(cor(data_for_analysis[rowSums(data_for_analysis$e)>0,]$e, method = "pearson"),color = col.pan,main = "Exclusion reads, Spearman corr.")

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
    ggtitle("Inclusion reads, mean by group") + theme_bw()
  
  
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
    ggtitle("Exclusion reads, mean by group") + theme_bw()
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

###Counts

p <- prcomp(t(gdata$cnts), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Counts, PCA")
pheatmap(cor(gdata$cnts[rowSds(gdata$cnts)>1,], method = "pearson"), color = col.pan,main = "Inclusion ratio, Spearman corr.")

rpkm <- calcRPKM(data,gdata,F)
nrow(rpkm$rpkm)
rpkm <- rpkm[rowSums(rpkm$rpkm)>0,]$rpkm
rpkm[is.na(rpkm)] <- 0
pheatmap(cor(rpkm, method = "pearson"), color = col.pan,main = "Inclusion ratio, Spearman corr.")

###STATISTICAL TESTS
samples <- c(seq(1:24))

transgenity <- c("tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt")

age <- c("60d", "60d", "60d", "60d", "60d", 
         "90d", "90d", "90d", "90d", 
         "120d", "120d", "120d", "120d", "120d", 
         "60d", "60d", "60d", "60d", "60d", 
         "120d", "120d", "120d", "120d", "120d")

groups <- c('tg1', 'tg1', 'tg1', 'tg1', 'tg1', 
            'tg2', 'tg2', 'tg2', 'tg2', 
            'tg3', 'tg3', 'tg3', 'tg3', 'tg3', 
            'wt1', 'wt1', 'wt1', 'wt1', 'wt1', 
            'wt3', 'wt3', 'wt3', 'wt3', 'wt3')


meta <- list(samples=samples,transgenity=transgenity, age = age, groups = groups)

data.f.glm = fitSAGLM(data.f,terms(x~transgenity),meta)
data.f.pv = calcSAPvalue(data.f.glm)
data.f$seg[,c("pvalue")]=data.f.pv[,"age"]
data.f$seg[,c("padj")] = p.adjust(data.f$seg[,c("pvalue")],method='BH')
data.f$seg$gene_symbol < mapIds(org.Mm.eg.db, 
                            keys=as.character(data.f$seg$gene_id), 
                            column="SYMBOL", 
                            keytype="ENSEMBL",
                            multiVals="first")

data.f$seg$gene_name <- mapIds(org.Mm.eg.db, 
                            keys=as.character(data.f$seg$gene_id), 
                            column="GENENAME", 
                            keytype="ENSEMBL",
                            multiVals="first")

cond1 <- meta$samples[which(meta$transgenity == "wt")]
cond2 <- meta$samples[which(meta$transgenity == "tg")]

ir <- data.frame(data.f$ir)
ir[is.na(ir)] <- 0
names(ir) <- meta$samples
ir$dpsi <- rowMeans(ir[names(ir) %in% cond2], na.rm = T) - rowMeans(ir[names(ir) %in% cond1], na.rm = T)
ir$pvalue <- data.f$seg$pvalue
ir$padj <- data.f$seg$padj
ir <- ir[abs(ir$dpsi)>0.1,]
ir <- ir[order(ir$padj,decreasing = F),]

