als <- read.xlsx('~/ALS-associated genes.xlsx', sheetIndex = 1, header = FALSE)
als <- as.data.frame(als$X1)
names(als) <- c("gene")

logCPM <- NULL
logCPM <- cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes)
nColCount <- ncol(logCPM)
logCPM <- as.data.frame(logCPM)
logCPM$Name <- mapIds(org.Mm.eg.db, 
                      keys=row.names(logCPM), 
                      column="GENENAME", 
                      keytype="ENSEMBL",
                      multiVals="first")

rownames(logCPM) <- make.names(y$genes$Symbol, unique=TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM)[nColCount+1] <- c("Name")

disset <- data.frame()
for (f in als$gene){
  sub <- logCPM[grepl(paste(f), rownames(logCPM), ignore.case = TRUE, fixed = FALSE),]
  disset <- rbind(sub, disset)
}
disset$Name <- NULL
disset[,25] <- NULL
disset <- t(scale(t(disset)))
pdf(file = "ALS_heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(disset, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()