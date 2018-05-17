library(edgeR)
library(biomaRt)
library(org.Dm.eg.db)
library(GO.db)
analyze_all_samples <- FALSE
qlm_test <- TRUE
cpm_cutoff <- -1
directory <- '~/counts/BAP.cutted.counts/'
setwd(directory)
gr_control <- c("CONTROL_ADULT")
gr_case <- c("10MM_ADULT")

if (analyze_all_samples == TRUE){
  sampleFiles <- grep('fly',list.files(directory),value=TRUE)
  sampleFiles
  sampleCondition <- c("old_k", "old_k", "old_n", "old_n", 
                       "old_f", "old_f", "old_m", "old_m", 
                       "new_f", "new_f", "new_k", "new_k", 
                       "new_m", "new_m")
  
  sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
  col <- as.vector(sampleTable$sampleName)
  y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
} else if (analyze_all_samples == FALSE){
  files_control <- grep(paste(gr_control),list.files(directory),value=TRUE)
  files_case <- grep(paste(gr_case),list.files(directory),value=TRUE)
  sampleFiles <- c(files_control, files_case)
  cond_control <- rep(paste(gr_control), length(files_control))
  cond_case <- rep(paste(gr_case), length(files_case))
  sampleCondition <- c(cond_control, cond_case)
  sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
  y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
}



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

# PLOTTING HTSEQ QUALITY BARPLOTS
row.names.remove <- c("__ambiguous", "__alignment_not_unique", "__no_feature", "__too_low_aQual", "__not_aligned" )


### DIFFEXPRESSION STATISTICAL ANALYSIS - EXACT NEGATIVE-
### BINOMIAL OR QLM TEST

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
  a <- calcNormFactors(a, method = "TMM", doWeighting = TRUE) 
  design <- model.matrix(~sampleTable$condition)
  a <- estimateDisp(a,design, tagwise = TRUE, trend.method = "loess") 
  fit <- glmQLFit(a,design, robust = TRUE, abundance.trend = TRUE) 
  qlf <- glmQLFTest(fit,coef=ncol(fit$design))
  et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
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
  #et_annot <- as.data.frame(et$table)
  et_annot <- as.data.frame(topTags(et, n = nrow(logCPM), adjust.method = "BH", sort.by = "PValue"))
  
  
}

CountsTable <- as.data.frame(a$counts)
counts_control <- CountsTable[,grep(gr_control, names(CountsTable))]
counts_case <- as.data.frame(CountsTable[,grep(gr_case, names(CountsTable))])
counts_control$rowsum.control <- rowSums(counts_control)
counts_case$rowsum.case <- rowSums(counts_case)
diff <- data.frame(counts_control$rowsum.control, counts_case$rowsum.case)
rownames(diff) <- rownames(CountsTable) 

fc <- et_annot[order(et_annot$logFC),]
lfgene <- rownames(fc[4,])
lfgene
c <- grep(paste(lfgene), rownames(diff))
dfc <- as.data.frame(diff[c,])
lfgenefc <- fc[4,1]
stat <- dfc$counts_control.rowsum.control > dfc$counts_case.rowsum.case
lfgene
if (stat == TRUE & lfgenefc < 0){
  print("No Correction Needed!")
} else {
  et_annot$logFC <- et_annot$logFC*(-1)
  print("Correction protocol executed, logFC have been inverted!")
}

taxon = 'Drosophila Melanogaster'
taxon = tolower(taxon)
tmp = unlist(strsplit(x = taxon, split = ' '))
dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
needed.attributes = c("ensembl_gene_id","external_gene_name", "transcript_length")

gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(et_annot), mart=mart)
gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]

###transcript length-logFC
et_annot.trlen <- et_annot
trlen.common <- intersect(rownames(et_annot.trlen), gmt_flt$ensembl_gene_id)
et_annot.trlen.common <- et_annot.trlen[(rownames(et_annot.trlen) %in% trlen.common),]
gmt_flt.common <- gmt_flt[(gmt_flt$ensembl_gene_id %in% trlen.common),]
et_annot.trlen.common <- et_annot.trlen.common[order(rownames(et_annot.trlen.common), decreasing = TRUE),]
gmt_flt.common <- gmt_flt.common[order(gmt_flt.common$ensembl_gene_id, decreasing = TRUE),]
trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM, et_annot.trlen.common$PValue, et_annot.trlen.common$FDR)
names(trlen.distrib) <- c("fc", "tl", "cpm", "pvalue", "FDR")

for.cor <- trlen.distrib[which(trlen.distrib$pvalue < 1),]
for.cor$tl <- log2(for.cor$tl)

fit <- lm(fc~tl, data = for.cor)
coefs <- coef(fit)
b1 <- round(coefs[2],2)
correl <- cor(for.cor$fc, for.cor$tl, method = "spearman")
png("Raw.png")
plot.fc.raw <- plot(fc ~ tl, data = for.cor, pch = ".",
                    xlab = "Log2(Transcript Length)",
                    ylab = "Log2(Fold Change)",
                    main = paste("Raw data:", "\n", "slope = ",b1, ", Transcript Length vs LogFC (all genes, logCPM>-1)", "\n", "Correlation = ", correl, sep = ""))
abline(fit)
dev.off()

tr.len.median <- median(trlen.distrib$tl)
m.low <- trlen.distrib[which(trlen.distrib$tl < tr.len.median),]$cpm 
m.low <- mean(10^m.low)
m.high <- trlen.distrib[which(trlen.distrib$tl > tr.len.median),]$cpm
m.high <- mean(10^m.high)
tr.length <- log2(trlen.distrib$tl)
assym.coef <- log2(m.low/m.high)
assym.coef

a <- DGEList(counts=y, group = sampleTable$condition) 
df <- as.data.frame(a$counts)

bins = c(0, 8, 9.6, 10.0, 10.35, 10.7, 11.05, 11.4, 11.8, 20)
adjusted_pseudo.counts.table.per.bin = vector(mode = 'list')
adjusted_pseudo.counts.table <- data.frame()
for(bin_n in 1:(length(bins) - 1)){
  current.counts.table = df[tr.length >= bins[bin_n] & tr.length < bins[bin_n + 1], ]
  d.bin = DGEList(counts = current.counts.table)
  d.bin = calcNormFactors(d.bin, method="TMM")
  print(d.bin$samples$norm.factors)
  current.adjusted_pseudo.counts.table = t(t(current.counts.table) / d.bin$samples$lib.size * mean(d.bin$samples$lib.size) / d.bin$samples$norm.factors )
  adjusted_pseudo.counts.table.per.bin[[bin_n]] = current.adjusted_pseudo.counts.table
  adjusted_pseudo.counts.table <- rbind(current.adjusted_pseudo.counts.table, adjusted_pseudo.counts.table)
}

adjusted_pseudo.counts.table = round(adjusted_pseudo.counts.table)
cs = colSums(adjusted_pseudo.counts.table)
dummy.row = t(as.data.frame(max(cs) - cs))
rownames(dummy.row) = 'dummy'
adjusted_pseudo.counts.table = rbind(adjusted_pseudo.counts.table, dummy.row)

d = DGEList(counts = adjusted_pseudo.counts.table)
d = calcNormFactors(d, method = 'none')
 

