###
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


if (analyze_all_samples == FALSE){
  dir.create(stattest)
  setwd(stattest)
} else if (analyze_all_samples == TRUE){
  dir.create("all")
  setwd("all")
}

if (qlm_test == TRUE){ 
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
  et <- exactTest(a)
} else if (qlm_test == FALSE){ 
  design <- model.matrix(~sampleTable$condition) 
  normalized_lib_sizes <- calcNormFactors(y, method = "TMM") 
  CountsTable <- as.data.frame(y$counts) 
  raw_counts <- as.data.frame(y$counts) 
  y <- estimateCommonDisp(y) 
  y <- estimateTagwiseDisp(y) 
  y <- estimateDisp(y, design = design) 
  nf <- exactTest(y) 
  no_filtered <- as.data.frame(nf$table) 
  keep <- rowSums(cpm(y) > cpm_cutoff) >= ncol(sampleTable) 
  cpm <- cpm(y) 
  cpm <- as.data.frame(cpm(y)) 
  cpm$rowsum <- rowSums(cpm) 
  y <- y[keep, , keep.lib.sizes=FALSE] 
  logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
  et <- exactTest(y) 
  top <- as.data.frame(topTags(et)) 
  et_annot <- as.data.frame(et$table) 
  et_annot_non_filtered <- as.data.frame(et$table) 
}
y$genes$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
y$genes$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(et_annot_non_filtered), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

et_annot$symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

et_annot$name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(et_annot), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")


et_annot$entrez <- mapIds(org.Mm.eg.db, 
                          keys=row.names(et_annot), 
                          column="ENTREZID", 
                          keytype="ENSEMBL",
                          multiVals="first")

et_annot$GOID <-     mapIds(org.Mm.eg.db, 
                            keys=row.names(et_annot), 
                            column="GO", 
                            keytype="ENSEMBL",
                            multiVals="first")

et_annot$term <- mapIds(GO.db, 
                        keys=et_annot$GOID, 
                        column="TERM", 
                        keytype="GOID",
                        multiVals="first")

et_annot$term <- as.character(et_annot$term)

et_annot_non_filtered$Symbol <- mapIds(org.Mm.eg.db, 
                                       keys=row.names(et_annot_non_filtered), 
                                       column="SYMBOL", 
                                       keytype="ENSEMBL",
                                       multiVals="first")
et_annot_non_filtered$Symbol <- mapIds(org.Mm.eg.db, 
                                       keys=row.names(et_annot_non_filtered), 
                                       column="ENTREZID", 
                                       keytype="ENSEMBL",
                                       multiVals="first")


top$Symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(top), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")
top$Name<- mapIds(org.Mm.eg.db, 
                  keys=row.names(top), 
                  column="GENENAME", 
                  keytype="ENSEMBL",
                  multiVals="first")


cpm <- as.data.frame(cpm(y))

cpm$Symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$Name <- mapIds(org.Mm.eg.db, 
                   keys=row.names(cpm), 
                   column="GENENAME", 
                   keytype="ENSEMBL",
                   multiVals="first")

cpm$entrez <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$GOID <-     mapIds(org.Mm.eg.db, 
                       keys=row.names(cpm), 
                       column="GO", 
                       keytype="ENSEMBL",
                       multiVals="first")

cpm$term <- mapIds(GO.db, 
                   keys=cpm$GOID, 
                   column="TERM", 
                   keytype="GOID",
                   multiVals="first")
cpm$term <- as.character(cpm$term)
et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))



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

goana_up <- goana(de = et_annot_high$entrez, species = "Mm")
go_up_500 <- topGO(goana_up, n=40000)
go_up_500$perc = (go_up_500$DE/go_up_500$N)*100
go_up_500$perc <- round(go_up_500$perc, digits = 4)

goana_down <- goana(de = et_annot_high$entrez, species = "Mm")
go_down_500 <- topGO(goana_down, n=40000)
go_down_500$perc = (go_down_500$DE/go_down_500$N)*100
go_down_500$perc <- round(go_down_500$perc, digits = 4)

write.csv(go_up_500, file = "go-up.csv")
write.csv(go_down_500, file = "go-down.csv")


keg_up <- kegga(de = et_annot_high$entrez, species="Mm")
tk_up <- topKEGG(keg_up, n=500)
#tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.csv(tk_up, file = "kegg_up.csv")

keg_down <- kegga(de = et_annot_low$entrez, species="Mm")
tk_down <- topKEGG(keg_down, n=500)
#tk_down <- subset(tk_down, tk_down$P.DE < 0.05)
tk_down$perc <- (tk_down$DE/tk_down$N)*100
tk_down <- tk_down[order(tk_down$perc, decreasing = TRUE),]
write.csv(tk_down, file = "kegg_down.csv")
beep()