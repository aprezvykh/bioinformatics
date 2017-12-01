pvalue_cutoff <- 0.05
logfchigh_cutoff <- 1
logfclow_cutoff <- -1
cpm_cutoff <- 0.5
directory <- '~/counts/ALS Mice/experimental/'
setwd(directory)
gr_control <- c("Control-1")
gr_case <- c("Tg-3")

files_control <- grep(paste(gr_control),list.files(directory),value=TRUE)
files_case <- grep(paste(gr_case),list.files(directory),value=TRUE)
sampleFiles <- c(files_control, files_case)
cond_control <- rep(paste(gr_control), length(files_control))
cond_case <- rep(paste(gr_case), length(files_case))
sampleCondition <- c(cond_control, cond_case)
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)

setwd("results")
stattest <- paste(gr_control, gr_case, sep = "-")
results.dir <- paste(directory, "results", sep = "")
setwd(results.dir)
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )


if (analyze_all_samples == FALSE){
  dir.create(stattest)
  setwd(stattest)
} else if (analyze_all_samples == TRUE){
  dir.create("all")
  setwd("all")
}


a <- DGEList(counts=y, group = sampleTable$condition) 
CountsTable <- as.data.frame(y$counts)
cpm <- cpm(y) 
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
cpm <- as.data.frame(cpm(y)) 
cpm$rowsum <- rowSums(cpm) 
keep <- rowSums(cpm > cpm_cutoff) >= ncol(sampleTable) 
logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
logCPM <- logCPM[!(row.names(logCPM) %in% row.names.remove), ]
logCPM <- logCPM[keep,]
a <- a[keep, , keep.lib.sizes=FALSE] 
a <- calcNormFactors(a, method = "TMM") 
design <- model.matrix(~sampleTable$condition) 
a <- estimateDisp(a,design) 
fit <- glmQLFit(a,design, robust = TRUE) 
qlf <- glmQLFTest(fit,coef=ncol(fit$design))
et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
et_annot_non_filtered <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
top <- as.data.frame(topTags(qlf, n = 20))
et <- qlf


et_annot$entrez <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot), 
                          column="ENTREZID", 
                          keytype="ENSEMBL",
                          multiVals="first")

et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]


counts_control <- CountsTable[,grep(gr_control, names(CountsTable))]
counts_case <- CountsTable[,grep(gr_case, names(CountsTable))]
counts_control$rowsum.control <- rowSums(counts_control)
counts_case$rowsum.case <- rowSums(counts_case)
diff <- data.frame(counts_control$rowsum.control, counts_case$rowsum.case)
rownames(diff) <- rownames(CountsTable) 

fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
c <- grep(paste(lfgene), rownames(diff))
dfc <- as.data.frame(diff[c,])
lfgenefc <- fc[4,1]
stat <- dfc$counts_control.rowsum.control > dfc$counts_case.rowsum.case

if (stat == TRUE & lfgenefc < 0){
  print("No Correction Needed!")
} else {
  et_annot$logFC <- et_annot$logFC*(-1)
  print("Correction protocol executed, logFC have been inverted!")
}

et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))


goana_up <- goana(de = et_annot_high$entrez, species = "Mm")
go_up <- topGO(goana_up, n=50000)
go_up$perc = (go_up$DE/go_up$N)*100
write.csv(go_up, "go-up.csv")

goana_down <- goana(de = et_annot_low$entrez, species = "Mm")
go_down <- topGO(goana_down, n=50000)
go_down$perc = (go_down$DE/go_down$N)*100
write.csv(go_up, "go-down.csv")


keg_up <- kegga(de = et_annot_high$entrez, species="Mm")
tk_up <- topKEGG(keg_up, n=1000)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
write.csv(tk_up, "kegg_up.csv")

keg_down <- kegga(de = et_annot_low$entrez, species="Mm")
tk_down <- topKEGG(keg_down, n=1000)
tk_down$perc <- (tk_down$DE/tk_down$N)*100
write.csv(tk_down, "kegg_down.csv")

df_high <- et_annot_high$entrez
x <- enrichPathway(gene=df_high, organism = "mouse", minGSSize=5, readable = TRUE, pvalueCutoff = 1)
write.csv(x, "r_up.csv")

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "mouse", minGSSize=gs_size, readable = TRUE, pvalueCutoff = 1)
write.csv(x, "r_down.csv")



###GO HEATMAP


control3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Control-3/go-up.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-1/go-up.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-2/go-up.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-3/go-up.csv")

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
rownames(df) <- df$term
df$term <- NULL
df$ont <- NULL
df <- as.matrix(df)
#df <- scale(df)
setwd("~/counts/ALS Mice/experimental/results/all/")
col.pan <- colorpanel(100, "white", "red")
pdf(file = "GO HEATMAP COUNT.pdf", width = 12, height = 17, family = "Helvetica")
pheatmap(df, col = col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,10), lhei=c(2,10), lwid=c(4,4), main = "GO terms", border_color = NA)
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


control3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Control-3/go-down.csv")
tg1 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-1/go-down.csv")
tg2 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-2/go-down.csv")
tg3 <- read.csv("~/counts/ALS Mice/experimental/results/Control-1-Tg-3/go-down.csv")

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
