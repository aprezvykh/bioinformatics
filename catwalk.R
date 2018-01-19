setwd("~/")
library(xlsx)
library(ggbiplot)
library(ggplot2)
library(dplyr)
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

pca <- df[,8:180]
pca$NA..1 <- NULL
pca$NA. <- NULL
pca$X169 <- NULL
pca <- pca[complete.cases(pca)]
p <- prcomp(pca, center = TRUE, scale. = TRUE)


ggbiplot(p, groups = df$exp, var.axes = FALSE, ellipse = TRUE)
