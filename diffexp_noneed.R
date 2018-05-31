### DESEQ2 PART (ADDITIONAL)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
res <- results(dds, tidy = FALSE )
rld<- rlogTransformation(dds, blind=TRUE)
setwd(results.dir)

if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


res <- as.data.frame(res)
png(file = "PCAPlot.png", width = 1024, height = 768)
pca <- plotPCA(res, intgroup=c('condition'))
print(plotPCA(rld, intgroup=c('condition')))
dev.off()

pdf(file = "MAplot(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotMA(dds,ylim=c(-10,10),main='DESeq2')
dev.off()

pdf(file = "Dispestimate(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotDispEsts(dds)
dev.off()

pdf(file = "Sparsity(DESeq2).pdf", width = 12, height = 17, family = "Helvetica")
plotSparsity(dds)
dev.off()

if (logging == TRUE){
  sink(type="message")
  close(zz)
}


### COMPARE DESEQ2 AND EDGER

res_df <- as.data.frame(res)
res_nf <- as.data.frame(res)
res_df <- res_df[complete.cases(res_df), ]
res_df <- as.data.frame(subset(res_df, res_df$baseMean > baseMean_cutoff))
res_df <- as.data.frame(subset(res_df, res_df$padj < pvalue_cutoff))
res_df <- as.data.frame(subset(res_df, res_df$log2FoldChange > logfchigh_cutoff | res_df$log2FoldChange < logfclow_cutoff))

res_df_high <- as.data.frame(subset(res_df, res_df$log2FoldChange > 0))
res_df_low <- as.data.frame(subset(res_df, res_df$log2FoldChange < 0))

com_deseq_edger_high <- as.data.frame(intersect(rownames(res_df_high), rownames(et_annot_high)))
com_deseq_edger_low <- as.data.frame(intersect(rownames(res_df_low), rownames(et_annot_low)))

proc_high<- (nrow(com_deseq_edger_high)/nrow(et_annot_high))*100
proc_low <- (nrow(com_deseq_edger_low)/nrow(et_annot_low))*100

sum <- NULL
sum <- data.frame(paste(proc_high), paste(proc_low))
names(sum) <- c("Upreg common genes", "Downreg common genes")

sign_proc_edger <- (nrow(et_annot)/nrow(et_annot_non_filtered))*100
sign_proc_deseq <- (nrow(res_df)/nrow(res_nf))*100

aproc <- data.frame(paste(sign_proc_edger), paste(sign_proc_deseq))

names(aproc) <- c("edgeR, significance %", "DESeq2, significance %")

write.xlsx(sum, file = "DESeq and edgeR comparsion.xlsx", sheetName = "Common genes of upreg and downreg", append = TRUE)
write.xlsx(aproc, file = "DESeq and edgeR comparsion.xlsx", sheetName = "Percent of significance, w.filtering", append = TRUE)

### COMPARE EDGER AND DESEQ2
com <- as.data.frame(intersect(rownames(res_nf), rownames(et_annot_non_filtered)))
names(com) <- c("sas")
edger_com <- data.frame()
deseq_com <- data.frame()

edger_com <- et_annot_non_filtered[(rownames(et_annot_non_filtered))%in% com$sas ,]
deseq_com <- res_nf[(rownames(res_nf))%in% com$sas ,]

common <- data.frame()
common <- data.frame(edger_com$logFC, deseq_com$log2FoldChange)
rownames(common) <- rownames(edger_com)
common <- transform(common, SD=apply(common,1, sd, na.rm = TRUE))

common <- common[order(common$SD, decreasing = TRUE),]

png(file = "Deviance between edgeR and deseq2.png")
plot(common$SD, type= "l", main= "Deviance between edgeR and deseq2(should be approximately zero)")
dev.off()

deviant <- common[seq(1:10),]

deviant$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(deviant), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")
deviant$Name <- mapIds(org.Ce.eg.db, 
                       keys=row.names(deviant), 
                       column="GENENAME", 
                       keytype="WORMBASE",
                       multiVals="first")

write.xlsx(deviant, file = "Top 10 deviant genes between deseq2 and edgeR.xlsx")


### SOME NUMBERS SHOULD BE FIXED dependent do dataset
dir.create("GO boxplots")
setwd("GO boxplots")

goana_straight <- goana(de = et_annot$entrez, species = "Ce")
go_all_1000 <- topGO(goana_straight, n=1000)
go_for_test <- topGO(goana_straight, n=1000)
go_all_1000 <- go_all_1000[!(go_all_1000$Ont %in% remove), ] 
go_all_1000$perc = (go_all_1000$DE/go_all_1000$N)*100
go_all_1000 <- subset(go_all_1000, go_all_1000$P.DE < 0.05)
go_all_1000 <- subset(go_all_1000, go_all_1000$perc > 20)
#go_all_1000 <- subset(go_all_1000, go_all_1000$N < 100)
go_all_1000 <- subset(go_all_1000, go_all_1000$DE <= 10)
#go_all_1000 <- subset(go_all_1000, go_all_1000$DE > 5)
go_all_1000 <- go_all_1000[order(go_all_1000$P.DE, decreasing = TRUE),]
go_all_1000 <- go_all_1000[seq(1:20),]
go_all_1000


for (i in 1:nrow(go_all_1000)){
  
  str <- go_all_1000[i,]
  f <- rownames(str)
  plot.name <- str$Term
  plot.ont <- str$Ont
  plot.perc <- round(str$perc, digits = 2)
  plot.pval <- round(str$P.DE, digits = 6)
  z <- allGO[[f]]
  z <- as.data.frame(z)
  thm <- cpm[(cpm$entrez %in% z$z),]
  rownames(thm) <- thm$Symbol
  thm$Symbol <- NULL
  thm$Name <- NULL
  thm$GOID <- NULL
  thm$term <- NULL
  thm$entrez <- NULL
  colnames(thm) <- sampleTable$condition
  thm <- log2(thm)
  thm <-as.data.frame(t(thm))
  thm$Ñondition <- rownames(thm)
  thm <- melt(thm, id.vars = "Ñondition")
  
  g12 <- ggplot(thm, aes(x = variable, y = value)) +
    geom_boxplot(data = thm, aes(fill = Ñondition), alpha = 0.5) + 
    scale_x_discrete(name = "Experimental Groups") + 
    scale_y_continuous(name = "logCPM") + 
    theme_bw() +
    geom_signif(comparisons = list(c("tg_2", "tg_3")), map_signif_level = TRUE) +
    #facet_wrap(~variable, scales="free") + 
    theme(axis.title.x=element_blank())
  g12 <- g12 + ggtitle(paste(f, ":", plot.name, "\n", "Ontology:", plot.ont, "\n", "Percent of significant genes: ", plot.perc, "%", "\n", "Term p-value:", plot.pval)) + theme(plot.title = element_text(hjust = 0.5, size = 10)) 
  ggsave(filename = paste(i, "_", plot.name, ".pdf", sep = ""), plot = g12, width = 20, height = 20, units = "cm")
  
}  


beep(sound = "coin")

################
setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

dir.create("GO K-means")
setwd("GO K-means")

go_all_1000 <- topGO(goana_straight, n=10000, ontology = "BP")
go_all_1000$g <- rownames(go_all_1000)
cpm <- as.data.frame(cpm(y))
names(cpm) <- sampleTable$condition
cpm <- cpm[keep,]
cpm$Symbol <- mapIds(org.Ce.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="WORMBASE",
                     multiVals="first")

cpm$entrez <- mapIds(org.Ce.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="WORMBASE",
                     multiVals="first")

sub <- go_all_1000[(go_all_1000$Term %in% go_up_30),]

for (f in rownames(sub)){
  print(f)
  z <- allGO[paste(f)]
  n <- go_all_1000[which(go_all_1000$g == f),]
  plot.name <- n$Term
  x <- cpm[which(cpm$entrez %in% z[[f]]),]
  e <- et_annot[(rownames(et_annot) %in% rownames(x)),]
  mean.1 <- apply(x[,grep(paste(gr_control), names(x))], 1, mean)
  mean.2 <- apply(x[,grep(paste(gr_case), names(x))], 1, mean)
  df <- data.frame()
  df <- data.frame(mean.1, mean.2)
  df$mean.1 <- log10(df$mean.1)
  df$mean.2 <- log10(df$mean.2)
  rownames(df) <- make.names(x$Symbol, unique = TRUE)
  cluster <- kmeans(df, 3, nstart = 20)
  df$gene <- rownames(df)
  df$Cluster <- as.factor(cluster$cluster)
  df$pval <- e$FDR
  df$lfc <- e$logFC
  df$Log2FoldChange <- ifelse(df$lfc > 0, "Positive", "Negative")
  
  g11 <- ggplot(df) + geom_point(aes(x = mean.1, y = mean.2, color = Cluster,  shape = Log2FoldChange)) + 
    theme_bw() + 
    geom_smooth(aes(x = mean.1, y = mean.2),method = "lm", se = FALSE, color = "red") + 
    geom_smooth(aes(x = mean.1, y = mean.2), se = FALSE, color = "green") + 
    geom_text(aes(x = mean.1, y = mean.2, label = ifelse(abs(lfc) > 2 & pval < 0.01, as.character(gene), ''), color = Cluster), hjust = 0, vjust = 0, size = 3) + 
    scale_x_continuous(name = paste("Genes in ", gr_control, " log10(CPM)", sep = "")) +
    scale_y_continuous(name = paste("Genes in ", gr_case, " log10(CPM)", sep = "")) + 
    ggtitle(paste("K-means clustering", "\n", f, ":", plot.name)) + theme(plot.title = element_text(hjust = 0.5, size = 10)) + 
    geom_abline(aes(intercept = 0, slope = 1), color = "blue") + 
    scale_colour_manual(name="Lines", values=c("green", "blue", "red"))
  ggsave(filename = paste(plot.name, ".png", sep = ""), plot = g11, width = 20, height = 20, units = "cm")
  
  
}


### GRAPHICAL SUMMARY


if (summary == TRUE){
  setwd(results.dir)
  if (analyze_all_samples == TRUE){
    setwd("all")
  } else {
    setwd(stattest)
  }
  
  filename <- "graphical summary.xlsx"
  wb <- createWorkbook(type="xlsx")
  sheet <- createSheet(wb, sheetName = "Model distributions")
  sheet1 <- createSheet(wb, sheetName = "DE visualization")
  sheet2 <- createSheet(wb, sheetName = "distributions")
  sheet3 <- createSheet(wb, sheetName = "enrichments")
  
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g)})
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g11)})
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g2)}, startCol = 18, startRow = 1)
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g3)}, startRow = 1, startCol = 10)
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g4)}, startRow = 25, startCol = 10)
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g5)}, startRow = 1, startCol = 18)
  xlsx.addPlot(wb, sheet, plotFunction = function(){print(g6)}, startRow = 25, startCol = 18)
  
  xlsx.addPlot(wb, sheet1, plotFunction = function(){print(volc)})
  xlsx.addPlot(wb, sheet1, plotFunction = function(){print(pca)})
  xlsx.addPlot(wb, sheet1, plotFunction = function(){heatmap.2(x, col=col.pan, Rowv=TRUE, scale="none",
                                                               trace="none", dendrogram="both", cexRow=1, cexCol=1, density.info="none",
                                                               margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")}, startCol = 10, startRow = 1)
  xlsx.addPlot(wb, sheet1, plotFunction = function(){plotSmear(qlf, de.tags = rownames(tt$table), lowess = TRUE, smooth.scatter = TRUE)}, startCol = 10, startRow = 25)
  xlsx.addPlot(wb, sheet1, plotFunction = function(){heatmap.2(logCPMpval, col=col.pan, Rowv=TRUE, scale="none",
                                                               trace="none", dendrogram="both", cexRow=0.5, cexCol=1, density.info="none",
                                                               margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top FDR genes, p < 0.05")}, startCol = 18, startRow = 1)
  xlsx.addPlot(wb, sheet1, plotFunction = function(){heatmap.2(top100cpm, col=col.pan, Rowv=TRUE, scale="none",
                                                               trace="none", dendrogram="both", cexRow=0.5, cexCol=1, density.info="none",
                                                               margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top FDR genes, p < 0.05")}, startCol = 18, startRow = 25)
  
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d1)})
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d2)}, startRow = 1, startCol = 9)
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d3)}, startRow = 25, startCol = 2)
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d4)}, startRow = 25, startCol = 9)
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d5)}, startRow = 50, startCol = 2 )
  xlsx.addPlot(wb, sheet2, plotFunction = function(){print(d6)}, startRow = 50, startCol = 9)
  xlsx.addPlot(wb, sheet3, plotFunction = function(){print(g_u)})
  xlsx.addPlot(wb, sheet3, plotFunction = function(){print(g_d)}, startCol = 10, startRow = 1)
  xlsx.addPlot(wb, sheet3, plotFunction = function(){print(r.bp.up)}, startCol = 1, startRow = 25)
  xlsx.addPlot(wb, sheet3, plotFunction = function(){print(r.bp.down)}, startCol = 10, startRow = 25)
  
  saveWorkbook(wb, "graphical summary.xlsx")
  xlsx.openFile("graphical summary.xlsx")
}






r <- grep("WBGene00000733", rownames(cpm), ignore.case = TRUE)
thm <- cpm[r,]
rownames(thm) <- thm$Symbol
plot.name <- as.character(rownames(thm))
plot.description <- as.character(thm$Name)
plot.term <- as.character(thm$term)
thm$Symbol <- NULL
thm$Name <- NULL
thm$GOID <- NULL
thm$term <- NULL
thm$entrez <- NULL
colnames(thm) <- sampleCondition
thm <-as.data.frame(t(thm))
thm$Condition <- rownames(thm)
thm$Group <- thm$Condition
names(thm) <- c("gene", "Condition", "Group")

g <- ggplot(thm, aes(x = Condition, y = gene)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5, coef = 1) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  geom_line(aes(x = Condition, y = mean(gene))) + 
  scale_x_discrete(name = "Experimental Groups") + 
  scale_y_continuous(name = "Log10(Counts per million)") + 
  theme_bw() + 
  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
  #geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
  ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
  theme(plot.title = element_text(hjust = 0.5))



g

ggsave(filename = paste(plot.name, "pdf", sep = "."), plot = g, height = 20, width = 20, units = "cm")