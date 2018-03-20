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

g <- grep("Mean", colnames(df), ignore.case = FALSE)

m <- df[,g]

pca <- m[,1:4]
pca <- pca[complete.cases(pca),]
p <- prcomp(pca, center = TRUE, scale. = TRUE)
lab.pca <- as.data.frame(t(df))

ggbiplot(p, var.axes = FALSE, ellipse = TRUE, choices = c(1,2), groups = df$exp)


d <- dist(m[,1:4])
fit <- cmdscale(d,eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Metric	MDS",	type="n")
text(x, y, labels = row.names(m), cex=.7)


c <- cor(m, method = "spearman")

pdf(file = "Corellation matrix Spearman111-MICE.pdf", width = 10, height = 10)
heatmap.2(c, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=0.4, cexCol=0.4, density.info="none",
          margin=c(12,12), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")
dev.off()



ñ.df <- as.data.frame(c)

