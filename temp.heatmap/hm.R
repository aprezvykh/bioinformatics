library(gplots, quietly = TRUE)
library(pheatmap, quietly = TRUE)
x <- read.table("~/adm/temp.heatmap/temp.txt")

col.pan <- colorpanel(80, "blue", "white", "red")
dat <- x$V3
dat <- gsub("°C", "", dat)
dat <- gsub("+", "", dat)
dat <- as.numeric(dat)
mat <- matrix(dat, nrow = 8, ncol = 8)
colnames(mat) <- c("CPU #1", "CPU #2", "CPU #3", "CPU #4", 
                   "CPU #5", "CPU #6", "CPU #7", "CPU #8")
rownames(mat) <- as.character(seq(1:8))
d <- system("date", intern = TRUE)
png(paste(d, ".png", sep = ""))
pheatmap(mat, cluster_cols = F, cluster_rows = F, col = col.pan, main = "CPU temperature", breaks=c(20:80))
dev.off()
