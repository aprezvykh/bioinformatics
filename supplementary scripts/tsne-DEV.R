install.packages("Rtsne")
library("Rtsne")
library(org.Mm.eg.db)
library(GO.db)
library(ggplot2)
library(cluster) 
cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv", header = TRUE)
rownames(cpm) <- make.names(cpm$X)
cpm$X <- NULL

t <- Rtsne(cpm, verbose = TRUE, check_duplicates = FALSE)

df <- as.data.frame(t$Y)
fit <- kmeans(df, 8)
clusplot(df, fit$cluster, color=TRUE, shade = TRUE, lines = 2)







