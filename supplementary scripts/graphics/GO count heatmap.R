library(gplots)
library(pheatmap)
library(xlsx)
control3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Control-3/go-up.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-1/go-up.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-2/go-up.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-3/go-up.csv")
sel <- read.xlsx("~/counts/ALS Mice/experimental/results/all/GO_selected (1).xlsx", sheetIndex = 1)
#control3 <- subset(control3, P.DE < 0.05)
#tg1 <- subset(tg1, P.DE < 0.05)
#tg2 <- subset(tg2, P.DE < 0.05)
#tg3 <- subset(tg3, P.DE < 0.05)
length(intersect(sel$GO.ID, tg2$X))/nrow(sel)

i <- intersect(control3$Term, tg1$Term)
j <- intersect(tg2$Term, tg3$Term)
k <- intersect(i,j)

c_control3 <- control3[(control3$Term %in% k),]
c_tg1 <- tg1[(tg1$Term %in% k),]
c_tg2 <- tg2[(tg2$Term %in% k),]
c_tg3 <- tg3[(tg3$Term %in% k),]

c_control3 <- c_control3[(order(c_control3$Term)),]
c_tg1 <- c_tg1[order(c_tg1$Term),]
c_tg2 <- c_tg2[order(c_tg2$Term),]
c_tg3 <- c_tg3[order(c_tg3$Term),]



df <- data.frame(c_control3$Term, c_control3$Ont, c_control3$perc, c_tg1$perc, c_tg2$perc, c_tg3$perc, c_control3$N)
names(df) <- c("term", "ont", "Control-3", "Tg-1", "Tg-2", "Tg-3", "N")
df <- df[which(df$ont == "BP"),]
df <- subset(df, N < 200)
df <- subset(df, N > 10)
df <- df[order(df$`Tg-3`, decreasing = TRUE),]
df <- df[seq(1:100),]
df <- df[complete.cases(df),]
df$N <- NULL
rownames(df) <- df$term

df$term <- NULL
df$ont <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "GO HEATMAP COUNT.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=FALSE, scale="none",
          trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "GO terms", border_color = NA 
          )
dev.off()



###KEGG HEATMAP


control3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Control-3/kegg_up.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-1/kegg_up.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-2/kegg_up.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-3/kegg_up.csv")

#control3 <- subset(control3, P.DE < 0.05)
#tg1 <- subset(tg1, P.DE < 0.05)
#tg2 <- subset(tg2, P.DE < 0.05)
#tg3 <- subset(tg3, P.DE < 0.05)

i <- intersect(control3$Pathway, tg1$Pathway)
j <- intersect(tg2$Pathway, tg3$Pathway)
k <- intersect(i,j)

c_control3 <- control3[(control3$Pathway %in% k),]
c_tg1 <- tg1[(tg1$Pathway %in% k),]
c_tg2 <- tg2[(tg2$Pathway %in% k),]
c_tg3 <- tg3[(tg3$Pathway %in% k),]

c_control3 <- c_control3[(order(c_control3$Pathway)),]
c_tg1 <- c_tg1[order(c_tg1$Pathway),]
c_tg2 <- c_tg2[order(c_tg2$Pathway),]
c_tg3 <- c_tg3[order(c_tg3$Pathway),]



df <- data.frame(c_control3$Pathway, c_control3$perc, c_tg1$perc, c_tg2$perc, c_tg3$perc, c_control3$N)
names(df) <- c("pathway", "Control-3", "Tg-1", "Tg-2", "Tg-3", "N")
df <- subset(df, N < 100)
df
df <- df[order(df$`Tg-3`, decreasing = TRUE),]
df <- df[seq(1:100),]
df$N <- NULL
rownames(df) <- df$pathway
df$pathway <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "KEGG HEATMAP COUNT.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
         trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
         margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "GO terms", border_color = NA)
dev.off()


### KEGG_PLOT
color <- c("purple", "green", "yellow", "red")
p <- df[which(df$`Tg-2` > 0 & df$`Tg-3` > 0),]
p <- p[5,]
plot.name <- p$pathway
plot.n <- p$N
rownames(p) <- plot.name
p$pathway <- NULL
p$N <- NULL
p <- as.data.frame(t(p))
p$cond <- rownames(p)
names(p) <- c("path", "Condition")
pdf(file = paste(plot.name, ".pdf", sep = ""), width = 10, height = 10, family = "Helvetica")
ggplot(p, aes(x = Condition, y = path, fill = Condition)) + 
       geom_bar(stat = "identity") + 
       scale_x_discrete(name = "Experimental groups") + 
       scale_y_continuous(name = "% of involved genes in pathway") +
       theme_bw() + ggtitle(paste(plot.name, "\n", "Total genes: ", plot.n)) + 
       theme(plot.title = element_text(hjust = 0.5))
dev.off()







###GO HEATMAP DOWNREG
control3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Control-3/go-down.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-1/go-down.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-2/go-down.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-3/go-down.csv")
sel <- read.xlsx("~/counts/ALS Mice/experimental/results/all/GO_selected (1).xlsx", sheetIndex = 2)
#control3 <- subset(control3, P.DE < 0.05)
#tg1 <- subset(tg1, P.DE < 0.05)
tg2 <- subset(tg2, P.DE < 0.05)
tg3 <- subset(tg3, P.DE < 0.05)

i <- intersect(control3$Term, tg1$Term)
j <- intersect(tg2$Term, tg3$Term)
k <- intersect(i,j)

c_control3 <- control3[(control3$Term %in% k),]
c_tg1 <- tg1[(tg1$Term %in% k),]
c_tg2 <- tg2[(tg2$Term %in% k),]
c_tg3 <- tg3[(tg3$Term %in% k),]

c_control3 <- c_control3[(order(c_control3$Term)),]
c_tg1 <- c_tg1[order(c_tg1$Term),]
c_tg2 <- c_tg2[order(c_tg2$Term),]
c_tg3 <- c_tg3[order(c_tg3$Term),]



df <- data.frame(c_control3$Term, c_control3$Ont, c_control3$perc, c_tg1$perc, c_tg2$perc, c_tg3$perc, c_control3$N)
names(df) <- c("term", "ont", "Control-3", "Tg-1", "Tg-2", "Tg-3", "N")
df <- df[which(df$ont == "BP"),]
df <- subset(df, N < 500)
df <- subset(df, N > 10)
df <- df[order(df$`Tg-3`, decreasing = TRUE),]
df <- df[seq(1:100),]
df$N <- NULL
df <- df[complete.cases(df),]
rownames(df) <- df$term
df$term <- NULL
df$ont <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "GO HEATMAP COUNT DOWNREG.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
         trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
         margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "GO terms downreg", border_color = NA)
dev.off()


###KEGG DOWNREG

control3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Control-3/kegg_down.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-1/kegg_down.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-2/kegg_down.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-3/kegg_down.csv")

#control3 <- subset(control3, P.DE < 0.05)
#tg1 <- subset(tg1, P.DE < 0.05)
#tg2 <- subset(tg2, P.DE < 0.05)
#tg3 <- subset(tg3, P.DE < 0.05)

i <- intersect(control3$Pathway, tg1$Pathway)
j <- intersect(tg2$Pathway, tg3$Pathway)
k <- intersect(i,j)

c_control3 <- control3[(control3$Pathway %in% k),]
c_tg1 <- tg1[(tg1$Pathway %in% k),]
c_tg2 <- tg2[(tg2$Pathway %in% k),]
c_tg3 <- tg3[(tg3$Pathway %in% k),]

c_control3 <- c_control3[(order(c_control3$Pathway)),]
c_tg1 <- c_tg1[order(c_tg1$Pathway),]
c_tg2 <- c_tg2[order(c_tg2$Pathway),]
c_tg3 <- c_tg3[order(c_tg3$Pathway),]



df <- data.frame(c_control3$Pathway, c_control3$perc, c_tg1$perc, c_tg2$perc, c_tg3$perc, c_control3$N)
names(df) <- c("pathway", "Control-3", "Tg-1", "Tg-2", "Tg-3", "N")
df <- df[order(df$`Tg-3`, decreasing = TRUE),]
df <- df[seq(1:100),]
df$N <- NULL
rownames(df) <- df$pathway
df$pathway <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "KEGG HEATMAP COUNT DOWNREG.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
         trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
         margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "GO terms", border_color = NA)
dev.off()


#############################################
control3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Control-3/kegg_up.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-1/kegg_up.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-2/kegg_up.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/enrichments/Control-1-Tg-3/kegg_up.csv")
sel <- read.xlsx("~/counts/ALS Mice/experimental/results/all/kegg selected.xlsx", sheetIndex = 1)


sel1 <- control3[(control3$X %in% sel$path),]
sel2 <- tg1[(tg1$X %in% sel$path),]
sel3 <- tg2[(tg2$X %in% sel$path),]
sel4 <- tg3[(tg3$X %in% sel$path),]
df <- data.frame(sel1$Pathway, sel2$perc, sel3$perc, sel4$perc)
names(df) <- c("pathway", "Tg-1", "Tg-2", "Tg-3")
rownames(df) <- df$pathway
df$pathway <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "kegg upregulated.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
         trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
         margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "Significantly upregulated KEGG terms in microglia", border_color = NA,
         cluster_rows=T, cluster_cols=F)
dev.off()

sel
