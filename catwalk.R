setwd("~/")
library(xlsx)
library(ggbiplot)
library(ggplot2)
library(dplyr)
library(ggsignif)
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



ggplot(df) + geom_boxplot(aes(x = df$exp, y = df$,
                              fill = df$exp)) + 
                              theme_bw()

df$
pca <- df[,8:15]
pca$NA..1 <- NULL
pca$NA. <- NULL
pca$X169 <- NULL
pca <- pca[complete.cases(pca),]
p <- prcomp(pca, center = TRUE, scale. = TRUE)

ggbiplot(p, groups = df$exp, var.axes = FALSE, ellipse = TRUE, choices = c(1,2))

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
