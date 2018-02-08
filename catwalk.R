setwd("~/")
library(xlsx)
library(ggbiplot)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
a <- read.xlsx(file = "trails_10 days.xlsx", sheetIndex = 8)
wt <- a[which(a$Group == "wt"),]
fus <- a[which(a$Group == "FUS"),]
control1 <- wt[which(wt$age.days. < 60),]
control1$exp <- c("Control-1")
control3 <- wt[which(wt$age.days. > 60),]
control3$exp <- c("Control-3")
tg1 <- fus[which(fus$age.days. < 60),]
tg1$exp <- c("Tg-1")
tg2 <- fus[which(fus$age.days. > 60 & fus$age.days. < 90),]
tg2$exp <- c("Tg-2")
tg3 <- fus[which(fus$age.days. > 90),]
tg3$exp <- c("Tg-3")

df <- bind_rows(control1, control3, tg1, tg2, tg3)


ggplot(df) + geom_boxplot(aes(x = df$exp, y = df$LH_BodySpeed_.cm.s._Mean, fill = df$exp))


pca <- df[,8:170]
pca <- as.data.frame(t(pca))

pca$NA..1 <- NULL
pca$NA. <- NULL
pca$X169 <- NULL
pca <- pca[complete.cases(pca),]
p <- prcomp(pca, center = TRUE, scale. = TRUE)
lab.pca <- as.data.frame(t(df))

ggbiplot(p, var.axes = FALSE, ellipse = TRUE, choices = c(1,2), labels = rownames(pca))



pdf(file = "PCA(3-4 comp).pdf", width = 10, height = 10)
ggbiplot(p, groups = df$exp, var.axes = FALSE, ellipse = TRUE, choices = c(3,4))
dev.off()

pdf(file = "PCA(5-6 comp).pdf", width = 10, height = 10)
ggbiplot(p, groups = df$exp, var.axes = FALSE, ellipse = TRUE, choices = c(5,6))
dev.off()

c <- df[,8:150]
c <- c[complete.cases(c),]
s <- cor(c, method = "spearman")

pdf(file = "Corellation matrix Spearman.pdf", width = 10, height = 10)
heatmap.2(s, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=0.4, cexCol=0.4, density.info="none",
          margin=c(12,12), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")
dev.off()


g <- grep("RF", colnames(df))
z <- df[,g]



z <- df[which(df$exp =="Control-1"),]
control1.means <- colMeans(z[,8:170])
z <- df[which(df$exp =="Control-3"),]
control3.means <- colMeans(z[,8:170])
z <- df[which(df$exp =="Tg-1"),]
tg1.means <- colMeans(z[,8:170])
z <- df[which(df$exp =="Tg-2"),]
tg2.means <- colMeans(z[,8:170])
z <- df[which(df$exp =="Tg-3"),]
tg3.means <- colMeans(z[,8:170])


means <- bind_rows(control1.means, control3.means, tg1.means, tg2.means, tg3.means)
rownames(means) <- c("Control-1", "Control-3", "Tg-1", "Tg-2", "Tg-3")

barplot(means$RF_Swing_.s._Mean, names.arg = rownames(means))

