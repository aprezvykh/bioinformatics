library(ggbiplot)
library(gplots)
library(ggplot2)
library(xlsx)
sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                     'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                     'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                     'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')

tg <- c('Control', 'Control', 'Control', 'Control', 'Control', 
        'Control', 'Control', 'Control', 'Control', 'Control', 
        'Tg', 'Tg', 'Tg', 'Tg', 'Tg', 
        'Tg', 'Tg', 'Tg', 'Tg', 
        'Tg', 'Tg', 'Tg', 'Tg', 'Tg')


col.pan <- colorpanel(100, "blue", "white", "red")
a <- read.xlsx("~/counts/ALS Mice/experimental/results/suppl/mge/CPMtrans.xlsx", sheetIndex = 1)

rownames(a) <- make.names(a$NA., unique = TRUE)
a$NA. <- NULL
length <- a$length
names(length) <- rownames(a)
a$length <- NULL
orig <- colnames(a)
colnames(a) <- sampleCondition

b <- data.matrix(a)
b <- b[,16:24]
b <- as.data.frame(b)

b$rs <- rowSums(b)
b <- b[which(b$rs > 24),]
b <- b[order(b$rs, decreasing = TRUE),]
b$rs <- NULL
b$length <- NULL
b <- b[1:70,]
hm <- t(scale(t(b)))


heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=0.6, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "", na.rm = TRUE)

df <- data.frame()
b$tg2_mean <- apply(b[,1:4], 1, mean)
b$tg2_sd <- apply(b[,1:4], 1, sd)                                                               
b$tg2_var <- apply(b[,1:4], 1, var)   

b$tg3_mean <- apply(b[,5:9], 1, mean)
b$tg3_sd <- apply(b[,5:9], 1, sd)                                                               
b$tg3_var <- apply(b[,5:9], 1, var)   

b$com_mean <- apply(b[,1:9], 1, mean)
b$com_sd <- apply(b[,1:9], 1, sd)
b$com_var <- apply(b[,1:9], 1, var)
b$com_sum <- apply(b[,1:9], 1, sum)
nrow(b)
b <- b[which(b$com_mean > 1),]
nrow(b)
b <- b[which(b$com_sum > 50),]
nrow(b)
b$lfc <- log2(b$tg3_mean/b$tg2_mean)
b$f <- b$tg3_sd/b$tg2_sd
b <- b[order(b$f, decreasing = TRUE),]

b <- b[is.finite(b$f),]

plot(b$lfc, b$f)


apply(b[10,1:9], 1, barplot)




