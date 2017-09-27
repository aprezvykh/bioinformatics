### NEW HEATMAP
y$genes$Symbol <- mapIds(org.Dm.eg.db, 
                         keys=row.names(y), 
                         column="SYMBOL", 
                         keytype="FLYBASE",
                         multiVals="first")

y$genes$Name <- mapIds(org.Dm.eg.db, 
                       keys=row.names(y), 
                       column="GENENAME", 
                       keytype="FLYBASE",
                       multiVals="first")
fit <- glmQLFit(y, robust=FALSE)
logCPM <- cpm(y, prior.count=2, log=TRUE)

rownames(logCPM) <- y$genes$Symbol 
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:100],]
logCPM <- t(scale(t(logCPM)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Main Heatmap_drosophila.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

dev.off()


### Gene set heatmap. Paste in let your word interested in!
let <- c("Hsp70")
logCPM <- NULL
logCPM <- cpm(y, prior.count=2, log=TRUE)
rownames(logCPM) <- y$genes$Symbol
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
sub <- logCPM[grep(paste(let), rownames(logCPM))]
sub <- t(scale(t(sub)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "paste(let).pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(sub, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()

### MULTIPLE BOXPLOTS OF INTERESTING GENES
logdf <- as.data.frame(sub)
for (i in seq(1:nrow(logdf))){
  png(file = paste(rownames(logdf[i,]), "_CPM.png", sep=""))
  boxplot.default(logdf[i,],outline = TRUE,  main = paste(rownames(logdf[i,])))
  dev.off()
}


### Multiple MDPlots

for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}

pdf(file = "MDplot_common.pdf", width = 12, height = 17, family = "Helvetica")
plotMD(tr, values=c(1,-1), col=c("red","blue"),
       legend="topright")
dev.off()

logdf <- as.data.frame(logCPM)
for (i in seq(1:nrow(logdf))){
  png(file = paste(rownames(logdf[i,]), "_CPM.png", sep=""))
  boxplot.default(logdf[i,],outline = TRUE,  main = paste(rownames(logdf[i,])))
  dev.off()
}
