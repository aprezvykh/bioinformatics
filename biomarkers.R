library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(dplyr)
col.pan <- colorpanel(100, "blue", "white", "red")
##устанавливаем домашнюю директорию
setwd("~/Documents/test_task/")
#cчитываем данные, так как у нас тут Nan вкачестве проебанных ячеек, так и пишем
df_main <- read.delim("Hugo.tsv", na.strings = "NaN")
biomarkers <- read.delim("biomarkers_info.tsv")
df_big <- read.delim("GSE91061.tsv", na.strings = "NaN")
pat <- read.delim("Hugo_annotation.tsv", header = FALSE)
#красим пациентов на основании ответа
pat$col <- ifelse(pat$V2 == "R", "red", "blue")
pat <- pat[order(pat$V1),]
#чистим датафреймы от говна
df_main <- df_main[complete.cases(df_main),]
#пишем названия маркеров в отдельную переменную, саму колонку стираем
#ассоциируем их с цветаами, красный - повышение, синий - падение
names.bio.main <- data.frame(df_main$Biomarker)
colnames(names.bio.main) <- c("names")
names.bio.main$color <- ifelse(names.bio.main$names == "up", "red", "blue")
df_main$Biomarker <- NULL

#переводим все в числа
df_main <-as.data.frame(sapply(df_main, as.numeric))
rownames(df_main) <- names.bio.main$names
biomarkers$col <- ifelse(biomarkers$Direction == "up", "red", "blue")

df_main <- data.frame(t(df_main))
df_main <- df_main[order(rownames(df_main)),]

#простейшие матрицы корелляций, по датасету. цвет - ответ на лечение
annot <- data.frame(pat$V1, pat$col)

df_main.cor <- cor(t(df_main), method = "spearman")
heatmap.2(df_main.cor, col = col.pan, Rowv=TRUE, Colv= TRUE, scale="none",
          trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6),
          main = "Corellation", RowSideColors = pat$col,
          ColSideColors = pat$col)
legend("topright", legend = unique(pat$V2), col = pat$col, lty= 1, lwd = 5, cex=.7)


##сделаем анализ главных компонент по всем образцам, c окраской по response-no response
pr.main <- prcomp(df_main, center = TRUE)
plot(pr.main$x, 
     col = biomarkers$col,
     main = "Principle Component Analysis",
     pch = ".")
text(pr.main$x, labels = rownames(df_main), col = pat$col)


###дендрограмма
dd <- dist(scale(df_main), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc)



