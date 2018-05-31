library(RDAVIDWebService)
library(r2excel)
library(BatchJobs)
library(BiocParallel)
library(ENCODExplorer)
library(DT)
library(reshape2)
library(htmlwidgets)
library(visNetwork)
library(beepr)
library(AnnotationDbi)
library(Rcpp)
library(gplots)
library(edgeR)
library(org.Mm.eg.db)
library(gplots)
library(pheatmap)
library(xlsx)
library(gage)
library(gageData)
library(topGO)
library(ggplot2)
library(ReactomePA)
library(reshape)
require(clusterProfiler)
require(reactome.db)
library(data.table)
library(pathview)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(Matrix)
library(PANTHER.db)
library(GO.db)
library(Hmisc)
library(DESeq2)
library(checkmate)
library(biomaRt)
library(RcisTarget)
library(devtools)
library(enrichR)
library(zoo)
library(rvest)
library(XML)
library(plyr)
library(AnnotationDbi)
library(ggsignif)
library(grid)
library(pcaExplorer)
library(ggbiplot)
library(cowplot)
library(gapminder)
library(gridExtra)

### TFES
library(RcisTarget.mm9.motifDatabases.20k)
#data("mm9_10kbpAroundTss_motifRanking")
#data("mm9_direct_motifAnnotation")

###BIOCBARALLEL SETTINGS
default <- registered()
register(BatchJobsParam(workers = 10), default = TRUE)
options(MulticoreParam=quote(MulticoreParam(workers=8)))
param <- SnowParam(workers = 2, type = "SOCK")


  heatmaps <- TRUE
  custom_genes_plots <- FALSE
  analyze_all_samples <- TRUE
  disease_association <- FALSE
  kegg_plots <- TRUE
  panther_analysis <- TRUE
  deseq2_part <- TRUE
  qlm_test <- TRUE
  logging <- FALSE
  motiv <- TRUE
  boxplots <- TRUE
  biotype <- FALSE
  ### CONSTANTS BLOCK
  
  pvalue_cutoff <- 0.05
  logfchigh_cutoff <- 0
  logfclow_cutoff <- 0
  cpm_cutoff <- -1000
  gs_size <- 10
  diseases_set <- 50
  number_of_kegg_plots <- 100
  go_terms_set <- 50
  pathways_set <- 30
  genes_in_term <- 3
  filter_thresh <- 5
  baseMean_cutoff <- 1.5
  significant_enriched_motif <- 5
  go_heatmap_count <- 20
  stattest_number <- 1
  

#a <- read.xlsx("tests.xlsx", sheetIndex = 1)
#a <- data.frame(a$control, a$case)
#gr_control <- as.character(a[1,1])
#gr_case <- as.character(a[1,2])

  directory <- '~/counts/ALS Mice/experimental/'
  setwd(directory)
  gr_control <- c("Tg-1")
  gr_case <- c("Tg-2")
  
  ### BUILDING A SPECIFIC DESIGN TABLE
  if (logging == TRUE){
    zz <- file("error.log", open="wt")
    sink(zz, type="message")
  }
  
  col.pan <- colorpanel(100, "blue", "white", "red")
  ###DIRECTORY WHERE SAMPLES ARE LOCATED
  library(dplyr)
  
  if (analyze_all_samples == TRUE){
          sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
          sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                               'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                               'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                               'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                               'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')
    
  
          sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
          col <- as.vector(sampleTable$sampleName)
          names(col) <- c("yellow", "yellow", "yellow", "yellow", "yellow", 
                          "purple", "purple", "purple", "purple", "purple", 
                          "green", "green", "green", "green", "green", 
                          "red", "red", "red", "red", 
                          "blue", "blue", "blue", "blue", "blue")
  
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
    a <- calcNormFactors(a, method = "TMM") 
    design <- model.matrix(~sampleTable$condition) 
    a <- estimateDisp(a,design) 
    fit <- glmQLFit(a,design = design, robust = TRUE) 
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
    y <- calcNormFactors(y, method = "TMM") 
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

#write.csv(cpm, "~/counts/Early.markers/all_cpm.csv")
####Dispresion

cpm$rowsum <- NULL
cpm <- as.data.frame(cpm)
g_tg1 <- grepl("mouse_Control", colnames(cpm), ignore.case = T)
g_tg2 <- grepl("mouse_Tg", colnames(cpm), ignore.case = T)
tg1 <- cpm[,g_tg1]
tg2 <- cpm[,g_tg2]

tg1$sd_tg1 <- apply(tg1,1, sd, na.rm = TRUE)

tg2$sd_tg2 <- apply(tg2,1, sd, na.rm = TRUE)

sas <- data.frame(tg1$sd_tg1, tg2$sd_tg2)
rownames(sas) <- rownames(tg1)
sas <- sas[rownames(sas) %in% rownames(et_annot),]
et_annot$disp1 <- sas$tg1.sd_tg1
et_annot$disp2 <- sas$tg2.sd_tg2



et_annot$disp_div <- et_annot$disp2/et_annot$disp1
et_annot_disp <- et_annot[which(et_annot$disp_div > 0.98 & et_annot$disp_div < 1.02),]

et_annot_disp <- et_annot_disp[which(et_annot_disp$PValue < 0.05),]
et_annot_disp <- et_annot_disp[which(et_annot_disp$FDR < 0.05),]
et_annot <- et_annot_disp

plot(et$table$logCPM, y$tagwise.dispersion, pch = ".")

z <- as.data.frame(et)
pallete = c("#F46D43", "#66C2A5", "#cd8845", "#3288BD", "#a8bf32", "#5E4FA2", "#D53E4F", "#d6d639", "#8ed384", "#9E0142", "#ebba2f")
density.cols = colorRampPalette(pallete)(dim(y$counts)[2])
pdf(file = "All samples density.pdf", height = 10, width = 10, family = "Helvetica")
lty.n = 1
lty.count = 4
plot(density(log(y$counts[,1])),col=density.cols[1], xlim=c(-4,10), main = 'Log2(CPM) density, non-normalized', lty = lty.n)
for (x in 2:dim(y$counts)[2]){
  lty.n = (lty.n +1) %% lty.count 
  lines(density(log(y$counts[,x])),col=density.cols[x], lty = lty.n)
}
legend.x.pos = (par("usr")[2] - par("usr")[1])*0.73 + par("usr")[1]
legend.y.pos = (par("usr")[4] - par("usr")[3])*0.93 + par("usr")[3]
legend(legend.x.pos, legend.y.pos, legend = colnames(y$counts), col = density.cols, cex = 0.8, lty = 1:4)
dev.off()  

getwd()


pdf(file = "Average expression level difference.pdf", height = 10, width = 10, family = "Helvetica")
barplot(colSums(y$counts), legend = colnames(y$counts),
        col=seq(1:length(colnames(y$counts))),
        main = "Average expression level difference")
dev.off()



y$samples


libz <- data.frame(y$counts)
libz.colsums <- data.frame(colSums(libz))
lib <- data.frame(fit$samples$lib.size, sampleTable$condition)
names(lib) <- c("size", "group")
lib
g <- ggplot(lib) + geom_boxplot(aes(x = group, y = size, fill = group)) +
              scale_x_discrete(name = "Experimental Groups") + 
              scale_y_continuous(name = "Reads counts") + 
              theme_bw() + 
              ggtitle("Library Size") + 
              theme(plot.title = element_text(hjust = 0.5))
pdf(file = "Library Size.pdf", height = 10, width = 10, family = "Helvetica")
g
dev.off()


disp <- data.frame(fit$dispersion, et_annot_non_filtered$logFC)
names(disp) <- c("dispersion", "LogFC")
disp$threshold = as.factor(et_annot_non_filtered$FDR < 0.05)

g1 <- ggplot(disp) + geom_point(aes(x = dispersion, y = LogFC, color = threshold), alpha = 1/2) +
     scale_x_continuous(name = "Dispersion") + 
     scale_y_continuous(name = "Log-2 Fold Change") + 
     theme_bw() + 
     ggtitle("Dispersion") + 
     theme(plot.title = element_text(hjust = 0.5)) +
     geom_smooth(aes(x = dispersion, y = LogFC))
pdf(file = "Dispersion-LogFC.pdf", height = 10, width = 10, family = "Helvetica")
g1
dev.off()


log <- data.frame(a$AveLogCPM, et_annot_non_filtered$logFC)
names(log) <- c("logCPM", "logFC")
log$threshold = as.factor(et_annot_non_filtered$FDR < 0.05)
g2 <- ggplot(log) + geom_point(aes(x = logCPM, y = logFC, color = threshold), alpha = 1/2) +
  scale_x_continuous(name = "logCPM") + 
  scale_y_continuous(name = "Log-2 Fold Change") + 
  theme_bw() + 
  ggtitle("logCPM-logFC trend") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_smooth(aes(x = logCPM, y = logFC))

pdf(file = "LogCPM-logFC.pdf", height = 10, width = 10, family = "Helvetica")
g2
dev.off()


pv <- data.frame(-log10(et$table$PValue), -log10(et_annot_non_filtered$FDR))
names(pv) <- c("pvalue", "fdr")
pv$thershold <- as.factor(pv$fdr < 0.05)
g3 <- ggplot(pv, aes(x = pvalue, y = fdr)) + 
    geom_point(aes(x = pvalue, y = fdr, color = thershold), alpha = 1/10) +
    scale_x_continuous(name = "-log10(P-value(exactTest))") + 
    scale_y_continuous(name = "-log10(FDR(glmQLFit))") + 
    theme_bw() + 
    ggtitle("P-value/FDR") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_smooth(aes(x = pvalue, y = fdr), method = "lm", color = "blue")

pdf(file = "PValue-FDR.pdf", height = 10, width = 10, family = "Helvetica")
g3
dev.off()


fc <- data.frame(et$table$logFC, et_annot_non_filtered$logFC)
names(fc) <- c("et", "glm")

g4 <- ggplot(fc, aes(x = glm, y = et)) + 
  geom_point(aes(x = glm, y = et), alpha = 1/10) + 
  scale_x_continuous(name = "Log2FoldChange - GLM") + 
  scale_y_continuous(name = "Log2FoldChange - Fisher Exact Test") + 
  theme_bw() + 
  ggtitle("Log2FoldChange GLM/Fisher Exact Test") + 
  theme(plot.title = element_text(hjust = 0.5))
  
pdf(file = "LogFC GLM-ET.pdf", height = 10, width = 10, family = "Helvetica")
g4
dev.off()

plot()

fc <- data.frame(qlf$table$logCPM, a$tagwise.dispersion, qlf$table$logFC, qlf$table$PValue)
names(fc) <- c("cpm", "tag", "logfc", "pvalue")
fc$threshold <- as.factor(fc$pvalue < 0.05)

g5 <- ggplot(fc, aes(x = cpm, y = tag)) + 
  geom_point(aes(x = cpm, y = tag, color = threshold), alpha = 1) + 
  scale_x_continuous(name = "LogCPM") + 
  scale_y_continuous(name = "Tagwise dispresion") + 
  theme_bw() + 
  ggtitle("LogCPM/Tagwise dispersion") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_smooth(se = FALSE)

#pdf(file = "LogCPM/tagwise dispersion.pdf", height = 10, width = 10, family = "Helvetica")
g5
dev.off()


fc <- data.frame(qlf$table$logCPM, a$trended.dispersion)
names(fc) <- c("cpm", "trend")

g6 <- ggplot(fc, aes(x = cpm, y = trend)) + 
  geom_point(aes(x = cpm, y = trend), alpha = 1/10) + 
  scale_x_continuous(name = "LogCPM") + 
  scale_y_continuous(name = "Trended dispersion") + 
  theme_bw() + 
  ggtitle("LogCPM/dispersion trend") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_smooth(se = FALSE)

pdf(file = "LogCPM/trended dispersion.pdf", height = 10, width = 10, family = "Helvetica")
g6
dev.off()

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
et_annot_non_filtered$entrez <- mapIds(org.Mm.eg.db, 
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


top <- top[complete.cases(top),]
### CPM
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


###TopTags Boxplots

if (boxplots == TRUE){
  dir.create("TopTags Boxplots")
  setwd("TopTags Boxplots")
    for (f in rownames(top)){
      r <- grep(paste(f), rownames(cpm), ignore.case = TRUE)
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
      #pdf(file = paste(plot.name, "pdf", sep = "."), width = 10, height = 10, family = "Helvetica")
      #png(paste(plot.name, "-", src, ".png", sep=""), height = 600, width = 800)
      g <- ggplot(thm, aes(x = Condition, y = gene)) + 
        geom_boxplot(data = thm, aes(fill = Group), alpha = 0.5) + 
        scale_x_discrete(name = "Experimental Groups") + 
        scale_y_continuous(name = "Counts per million") + 
        theme_bw() + 
        geom_signif(comparisons = list(c(paste(gr_control)), paste(gr_case)), map_signif_level = TRUE) 
        
      g <- g + ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
        theme(plot.title = element_text(hjust = 0.5))
                    
        ggsave(filename = paste(plot.name, "pdf", sep = "."), plot = g)

    }
}


setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}
getwd()


### transcript length distribution
taxon = 'Mus Musculus'
taxon = tolower(taxon)
tmp = unlist(strsplit(x = taxon, split = ' '))
dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
mart <- useMart("ensembl", dataset=dataset.name)
needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene", "transcript_length")

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


trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM, et_annot.trlen.common$PValue, ed_annot.trlen.common$FDR)
names(trlen.distrib) <- c("fc", "tl", "cpm", "pvalue", "FDR")

##expressed
trlen.distrib <- trlen.distrib[which(trlen.distrib$cpm > 2),]
trlen.distrib$tl <- log2(trlen.distrib$tl)
fit <- lm(fc~tl, data = trlen.distrib)
rmse <- round(sqrt(mean(resid(fit)^2)), 2)
coefs <- coef(fit)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(fit)$r.squared, 2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                r^2 == .(r2) * "," ~~ RMSE == .(rmse))
png("Transcript length - LogFC distribution, expressed.png", height = 20, width = 20)
plot(fc ~ tl, data = trlen.distrib, pch = ".", xlab = "Log2(Transcript Length)", ylab = "Log2(Fold Change)", main = "Transcript Length vs LogFC (all genes, logCPM>2)")
abline(fit)
dev.off()

###significant
trlen.distrib <- trlen.distrib[which(trlen.distrib$pvalue < pvalue_cutoff && trlen.distrib$FDR < pvalue_cutoff),]
trlen.distrib$tl <- log2(trlen.distrib$tl)
fit <- lm(fc~tl, data = trlen.distrib)
rmse <- round(sqrt(mean(resid(fit)^2)), 2)
coefs <- coef(fit)
b0 <- round(coefs[1], 2)
b1 <- round(coefs[2],2)
r2 <- round(summary(fit)$r.squared, 2)

eqn <- bquote(italic(y) == .(b0) + .(b1)*italic(x) * "," ~~ 
                r^2 == .(r2) * "," ~~ RMSE == .(rmse))
png("Transcript length - LogFC distribution, significant.png", height = 20, width = 20)
plot(fc ~ tl, data = trlen.distrib, pch = ".", xlab = "Log2(Transcript Length)", ylab = "Log2(Fold Change)", main = "Transcript Length vs LogFC (significant)")
abline(fit)
dev.off()


### FILTRATION

et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))
et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]


if(nrow(et_annot < 50)){
  for(i in seq(1:2)){
    beep()
    Sys.sleep(0.1)
    beep()
    Sys.sleep(0.2)
    print("Too low number of differentialy expressed genes found!
           reduce LogFC cutoff or check your data!")  
  }
  break
}


###PRcomp of differentially expressed genes 
аdir.create("prcomp")
setwd("prcomp")
fp <- et_annot
fp$is.upreg <- ifelse(fp$logFC > 0, "Upregulated", "Downregulated")
fp$is.expressed <- ifelse(fp$logCPM > 1, "High-expressed", "Low-expressed")

p <- prcomp(fp[,1:5], center = TRUE, scale. = TRUE)

pdf(file = "PCA-plot dim 1-2-FoldChange.pdf", width = 12, height = 17, family = "Helvetica")
ggbiplot(p, choices = c(1,2), labels = fp$symbol, groups = fp$is.upreg)
dev.off()
pdf(file = "PCA-plot dim 1-2-CPM.pdf", width = 12, height = 17, family = "Helvetica")
ggbiplot(p, choices = c(1,2), labels = fp$symbol, groups = fp$is.expressed)
dev.off()

pdf(file = "PCA-plot dim 3-4-FoldChange.pdf", width = 12, height = 17, family = "Helvetica")
ggbiplot(p, choices = c(3,4), labels = fp$symbol, groups = fp$is.upreg)
dev.off()
pdf(file = "PCA-plot dim 3-4-CPM.pdf", width = 12, height = 17, family = "Helvetica")
ggbiplot(p, choices = c(3,4), labels = fp$symbol, groups = fp$is.expressed)
dev.off()

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}


dir.create("densities")
setwd("densities")
pdf(file = "Fold Change density.pdf", width = 12, height = 17, family = "Helvetica")
plot(density(et_annot$logFC), main = "Fold Change Density")
dev.off()
pdf(file = "logCPM density.pdf", width = 12, height = 17, family = "Helvetica")
plot(density(et_annot$logCPM), main = "logCPM Density")
dev.off()
pdf(file = "FDR density.pdf", width = 12, height = 17, family = "Helvetica")
plot(density(et_annot$FDR), main = "FDR Density")
dev.off()
pdf(file = "PValue density.pdf", width = 12, height = 17, family = "Helvetica")
plot(density(et_annot$PValue), main = "PValue Density")
dev.off()



### BIOTYPE ANNOT FIX IT!
if (biotype == TRUE){
    taxon = 'Mus Musculus'
    taxon = tolower(taxon)
    tmp = unlist(strsplit(x = taxon, split = ' '))
    dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
    mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
    use.official.gene.symbol <- TRUE
    needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene")
    if(use.official.gene.symbol == TRUE){
      gmt_flt = getBM(attributes=needed.attributes,filters="external_gene_name",values=et_annot$symbol, mart=mart)
      gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"external_gene_name"])),]
      rownames(gmt_flt) = gmt_flt[,"external_gene_name"]
    } else {
      gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(et_annot), mart=mart)
      gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
      rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]
    }
    
    # et_annot$biotype <- gmt_flt$gene_biotype
    gmt_flt_freq <- as.data.frame(table(unlist(gmt_flt$gene_biotype)))
    
    pdf(file = "Filtered genes biotype distribution.pdf", width = 12, height = 17, family = "Helvetica")
    bp <- ggplot(gmt_flt_freq, aes(x="", y=Freq, fill=Var1))+
      geom_bar(width = 1, stat = "identity")
    pie <- bp + coord_polar("y", start=0) + theme_bw() +theme(axis.title.x=element_blank(),
                                                                 axis.text.x=element_blank(),
                                                                 axis.ticks.x=element_blank(),
                                                                 axis.title.y=element_blank(),
                                                                 axis.text.y=element_blank(),
                                                                 axis.ticks.y=element_blank())
    
    pie
    dev.off()
}

### TESTING A HYPOTESIS
counts_control <- CountsTable[,grep(gr_control, names(CountsTable))]
counts_case <- as.data.frame(CountsTable[,grep(gr_case, names(CountsTable))])
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


filt <- et_annot[which(abs(et_annot$FDR) < 0.05),]


### Distribution, with moda and median
pdf(file = "Distributions.pdf", width = 10, height = 10, family = "Helvetica")
par(mfrow=c(3,2))
hist(et_annot_non_filtered$logCPM, main = "logCPM distribution, non-filtered", freq = FALSE, col = col.pan, xlab = "logCPM", breaks = 50)
hist(et_annot$logCPM, main = "logCPM distribution, filtered", freq = FALSE, col = col.pan, xlab = "logCPM", breaks = 50)
hist(et_annot_non_filtered$PValue, main = "p-value distribution, non-filtered", freq = FALSE, col = col.pan, xlab = "p-value", breaks = 50)
hist(et_annot$PValue, main = "p-value distribution, filtered", freq = FALSE, col = col.pan,  xlab = "p-value", breaks = 50)
hist(et_annot_non_filtered$logFC, main = "logFC distribution, non-filtered", freq = FALSE, col = col.pan, xlab = "LogFC", breaks = 50)
hist(et_annot$logFC, main = "logFC distribution, filtered", freq = FALSE, col = col.pan, xlab = "LogFC", breaks = 50)
dev.off()

d1 <- ggplot(data=et_annot_non_filtered, aes(et_annot_non_filtered$logCPM)) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "logCPM") + 
      ggtitle("logCPM, non filtered") + 
      theme(plot.title = element_text(hjust = 0.5))

d2 <- ggplot(data=et_annot, aes(et_annot$logCPM)) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "logCPM") + 
      ggtitle("logCPM,filtered") + 
      theme(plot.title = element_text(hjust = 0.5))

d3 <- ggplot(data=et_annot_non_filtered, aes(-log10(et_annot_non_filtered$PValue))) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "-log10(PValue)") + 
      ggtitle("-log10(PValue), non filtered") + 
      theme(plot.title = element_text(hjust = 0.5))

d4 <- ggplot(data=et_annot, aes(-log10(et_annot$PValue))) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "-log10(PValue)") + 
      ggtitle("-log10(PValue), filtered") + 
      theme(plot.title = element_text(hjust = 0.5))

d5 <- ggplot(data=et_annot_non_filtered, aes(et_annot_non_filtered$logFC)) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "logFC") + 
      ggtitle("logFC, non filtered") + 
      theme(plot.title = element_text(hjust = 0.5))

d6 <- ggplot(data=et_annot, aes(et_annot$logFC)) + 
      geom_histogram(binwidth = 0.05) + 
      scale_x_continuous(name = "logFC") + 
      ggtitle("logFC, filtered") + 
      theme(plot.title = element_text(hjust = 0.5))


setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else if (analyze_all_samples == FALSE){
  setwd(stattest)  
}

### Simple summary
all <- nrow(y$counts)
allpadj <- sum(et_annot$PValue < pvalue_cutoff, na.rm=TRUE)
avg_cpm <- mean(et_annot$logCPM)
up <- sum(et_annot$logFC > logfchigh_cutoff, na.rm=TRUE)
down <- sum(et_annot$logFC < logfclow_cutoff, na.rm=TRUE)
header <- c('all genes', 'mean of logCPM', 'padj<0,05', 'genes with > high', 'genes with < low')
meaning <- c(print(all), print(avg_cpm), print(allpadj), print(up), print(down))
df <- data.frame(header, meaning)

cpm_subset <- cpm[(rownames(cpm) %in% rownames(et_annot)),]
cpm_subset <- cpm_subset[,1:length(sampleCondition)]
m <- data.frame()
m <- cpm_subset
m <- cbind(et_annot, m)
# WRITE RESULTS
write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(top, file = "Results edgeR.xlsx", sheetName = "Top Tags (with FDR)", append = TRUE)
write.xlsx(m, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)
m <- NULL
##TEXT SUMMARY
et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))
if(nrow(et_annot_high > nrow(et_annot_low))){
  updown <- c(" Преобладает апререгуляция. ")
} else if(nrow(et_annot_high < nrow(et_annot_low))){
  updown <- c(" Преобладает даунререгуляция. ")
} else {
  updown <- c(" Равное количество апрегулированных и даунрегулированных генов. ")
}


if(nrow(et_annot) < 200){
  cond <- c(" Малое количество диффэспрессных генов, это скажется на обогащении. Рекомендуется уменьшить отсечку по logFC.")
} else {
  cond <- c(" Нормальное количество диффэкспрессных генов.")
}

nrow(subset(go_for_test, go_for_test$P.DE < 0.05))

p <-paste("Всего генов: ", nrow(y), "," , "из них экспрессных: ", nrow(logCPM), "(", (round(nrow(logCPM)/nrow(y)*100)), "%)", ", ",
      "достоверно диффэкспрессных: ", nrow(et_annot), " (", (round(nrow(et_annot)/nrow(logCPM)*100)), "% от всех экспрессных)", 
      " с отсечкой LogFC по модулю ", abs(logfchigh_cutoff), "," ," отсечкой p-value ", pvalue_cutoff,  ", ", "апрегулированных ",
      nrow(et_annot_high), ", ", "даунрегулированных ", nrow(et_annot_low), ".", cond, updown, sep="")



p

  cat(p, file = "sas.txt")


### ANNOTATE COUNTS
CountsTable$symbol <- mapIds(org.Mm.eg.db, 
                             keys=row.names(CountsTable), 
                             column="SYMBOL", 
                             keytype="ENSEMBL",
                             multiVals="first")

CountsTable$name <- mapIds(org.Mm.eg.db, 
                           keys=row.names(CountsTable), 
                           column="GENENAME", 
                           keytype="ENSEMBL",
                           multiVals="first")


### GO TESTS
data(go.sets.mm)
data(go.subs.mm)

foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE,set.size = c(10, 100), rank.test = TRUE)
lapply(gomfres, head)
gomfres <- as.data.frame(gomfres)
gomfres <- gomfres[complete.cases(gomfres), ]
write.xlsx(gomfres, file = "GO.xlsx", sheetName = "GO_MF", append = TRUE)

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- as.data.frame(gobpres)
gobpres <- gobpres[complete.cases(gobpres), ]
write.xlsx(gobpres, file = "GO.xlsx", sheetName = "GO_BP", append = TRUE)

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- as.data.frame(goccres)
goccres <- goccres[complete.cases(goccres), ]
write.xlsx(goccres, file = "GO.xlsx", sheetName = "GO_CC", append = TRUE)


### GOANA

###KEGG expression profile (without lfc, but with generatio)

kk <- enrichKEGG(gene = df_high, organism = "mmu", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

kk <- enrichKEGG(gene = df_low, organism = "mmu", pvalueCutoff = 0.05)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

###KEGGA
keg_com <- kegga(de = et_annot$entrez, species="Mm")
tk_common <- topKEGG(keg_com, n=100)
tk_common <- subset(tk_common, tk_common$P.DE < 0.05)
tk_common$perc <- (tk_common$DE/tk_common$N)*100
tk_common <- tk_common[order(tk_common$perc, decreasing = TRUE),]
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Mm")
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Mm")
tk_down <- topKEGG(keg_down, n=100)
tk_down <- subset(tk_down, tk_down$P.DE < 0.05)
tk_down$perc <- (tk_down$DE/tk_down$N)*100
tk_down <- tk_down[order(tk_down$perc, decreasing = TRUE),]
write.xlsx(tk_down, file = "kegga.xlsx", sheetName = "Downreg", append = TRUE)
getwd()
rownames(tk_common) <- substring(rownames(tk_common), 6)


write.xlsx(tk_up, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated KEGG pathways", append = TRUE)
write.xlsx(tk_down, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated KEGG pathways", append = TRUE)


# TOP 100 PVALUE GENES
logCPM <- as.data.frame(logCPM)
nrow(logCPM)
rownames(logCPM) <- make.names(y$genes$Symbol, unique = TRUE)
colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
o <- order(et$table$PValue)
logCPMpval <- logCPM[o[1:100],]
logCPMpval <- t(scale(t(logCPMpval)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 Pvalue Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPMpval, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top FDR genes, p < 0.05")
dev.off()

# TOP 100 LOGFC GENES
o <- order(et$table$logFC)
logCPMfc <- logCPM[o[1:100],]
logCPMfc <- logCPMfc[complete.cases(logCPMfc),]
logCPMfc <- t(scale(t(logCPMfc)))
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(file = "Top 100 logFC Heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(logCPMfc, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Top logCPM genes, p < 0.05", na.rm = TRUE)
dev.off()


### High Expressed Genes
topcpm <- as.data.frame(cpm(y))
topcpm$rowsum <- rowSums(topcpm)
topcpm <- topcpm[order(topcpm$rowsum, decreasing = TRUE),]
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[!(row.names(topcpm) %in% row.names.remove), ] 
topcpm <- topcpm[1:100,]

topcpm$rowsum <- NULL
topcpm <- as.data.frame(topcpm)
# colnames(topcpm) <- paste(y$samples$group, 1:2, sep="-")
topcpm$Symbol <- mapIds(org.Mm.eg.db, 
                        keys=row.names(topcpm), 
                        column="SYMBOL", 
                        keytype="ENSEMBL",
                        multiVals="first")
rownames(topcpm) <- make.names(topcpm$Symbol, unique=TRUE)
topcpm$Symbol <- NULL
top100cpm <- t(scale(t(topcpm)))
pdf(file = "All heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(top100cpm, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "",
          labCol = FALSE)
            
dev.off()
### heatmap with all samples
topcpm <- as.data.frame(cpm(y))
topcpm$rowsum <- rowSums(topcpm)
topcpm <- topcpm[!(row.names(topcpm) %in% row.names.remove), ] 
topcpm <- topcpm[order(topcpm$rowsum, decreasing = TRUE),]
topcpm$rowsum <- NULL
topcpm <- topcpm[complete.cases(topcpm), ]
topcpm <- topcpm[1:15000,]
topcpm <- t(scale(t(topcpm)))
pdf(file = "Top 10000 expressed genes.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(topcpm, col = col.pan, Rowv=TRUE, scale="column",
            trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Highest expressed genes",
            labRow = FALSE)
dev.off()

getwd()

### Multiple MDPlots
dir.create("mdplots")
setwd("mdplots")
for (f in 1:ncol(y)){
  png(file = paste(f, ".png", sep=""))
  plotMD(y, column=f)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}
pdf(file = "MDplot_common.pdf", width = 12, height = 17, family = "Helvetica")
plotMD(y, values=c(1,-1), col=c("red","blue"),
       legend="topright")
dev.off()

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}



# KEGG PLOTS
setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggreswr <- as.data.frame(keggres)
write.xlsx(keggreswr, file = "KEGG pathview.xlsx", sheetName = "KEGG", append = TRUE)

if (kegg_plots == TRUE){
keggrespathways = data.frame(id=rownames(keggres$greater),
  keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=number_of_kegg_plots) %>% 
  .$id %>% 
  as.character()
keggresids = substr(keggrespathways, start=1, stop=8)

plot_pathway = function(pid){
         pathview(gene.data=foldchanges, 
         pathway.id=pid, 
         species="mmu", 
         new.signature=FALSE)
}


plot_pathway("mmu03050")

detach("package:dplyr", unload=TRUE)
dir.create("kegg")
setwd("kegg")

tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


getwd()
dir.create("another kegg plots")
setwd("another kegg plots")
plot_pathway("mmu01100")

plot_pathway("mmu01100")

for (f in rownames(tk_common)){
  plot_pathway(f)
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


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


png(file = "PCAPlot.png", width = 1024, height = 768)
pca <- plotPCA(rld, intgroup=c('condition'))
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

deviant$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(deviant), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")
deviant$Name <- mapIds(org.Mm.eg.db, 
                         keys=row.names(deviant), 
                         column="GENENAME", 
                         keytype="ENSEMBL",
                         multiVals="first")

write.xlsx(deviant, file = "Top 10 deviant genes between deseq2 and edgeR.xlsx")
### MOTIV ENRICHMENT


if (motiv == TRUE){
      geneList1 <- et_annot$symbol
      geneLists <- list(geneListName=geneList1)
      motifRankings <- mm9_10kbpAroundTss_motifRanking
      
      motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                      motifAnnot_direct=mm9_direct_motifAnnotation)
      
      motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
      resultsSubset <- motifEnrichmentTable_wGenes_wLogo
      result <- datatable(resultsSubset[,-c("enrichedGenes", "TF_inferred"), with=FALSE], 
                          escape = FALSE, # To show the logo
                          filter="top", options=list(pageLength=5))
      
      
      saveWidget(result, file="TF enrichment.html")
      
      motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
      motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=mm9_direct_motifAnnotation)
      
      signifMotifNames <- motifEnrichmentTable$motif[1:3]
      incidenceMatrix <- getSignificantGenes(geneLists, 
                                             motifRankings,
                                             signifRankingNames=signifMotifNames,
                                             plotCurve=TRUE, maxRank=5000, 
                                             genesFormat="incidMatrix",
                                             method="aprox")$incidMatrix
      
      edges <- melt(incidenceMatrix)
      edges <- edges[which(edges[,3]==1),1:2]
      colnames(edges) <- c("from","to")
      
      
      motifs <- unique(as.character(edges[,1]))
      genes <- unique(as.character(edges[,2]))
      nodes <- data.frame(id=c(motifs, genes),   
                          label=c(motifs, genes),    
                          title=c(motifs, genes), # tooltip 
                          shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                          color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
      
      web <- visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                                     nodesIdSelection = TRUE)
      
      saveWidget(web, file="TF enrichment web.html")
      write.xlsx(motifEnrichmentTable_wGenes, file = "Transcription factor binding motif enrichment.xlsx")
}



### GO BOXPLOTS!
### SOME NUMBERS SHOULD BE FIXED dependent do dataset
dir.create("GO boxplots")
setwd("GO boxplots")

goana_straight <- goana(de = et_annot$entrez, species = "Mm")
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
  thm$Сondition <- rownames(thm)
  thm <- melt(thm, id.vars = "Сondition")
  
  g12 <- ggplot(thm, aes(x = variable, y = value)) +
    geom_boxplot(data = thm, aes(fill = Сondition), alpha = 0.5) + 
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
cpm$Symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$entrez <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="ENSEMBL",
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



### goplot


cpm <- as.data.frame(cpm(y, log = TRUE))

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






sub <- cpm[grepl("Hspa", cpm$Symbol),]


z <- grep("ENSMUSG00000091971", rownames(cpm))

dir.create("HSP")
setwd("HSP")
for (f in rownames(sub)){
  r <- grep(paste(f), rownames(cpm), ignore.case = TRUE)
  thm <- cpm[f,]
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
    geom_boxplot(data = thm, aes(fill = Group), alpha = 0.5) + 
    scale_x_discrete(name = "Experimental Groups") + 
    scale_y_continuous(name = "Counts per million") + 
    geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
    geom_signif(comparisons = list(c("Tg-1", "Tg-2")), map_signif_level = TRUE) +
    theme_bw()
  g <- g + ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
    theme(plot.title = element_text(hjust = 0.5))
  g
  ggsave(filename = paste(plot.name, "png", sep = "."), plot = g, height = 25, width = 25, units = "cm")
}







cpm <- as.data.frame(cpm(y))
write.csv(cpm, "~/counts/Early.markers/all_cpm.csv")

df <- cpm[which(rownames(cpm) == "ENSMUSG00000033777"),]
df <- as.data.frame(t(df))
df$names <- rownames(df)
df$group <- c("FUS-Control", "FUS-Control", "FUS-Control", "FUS-Control", "FUS-Control", 
              "FUS-Control", "FUS-Control", "FUS-Control", "FUS-Control", "FUS-Control", 
              "TDP-Control", "TDP-Control", "TDP-Control", "TDP-Control",
              "TDP-Tg", "TDP-Tg", "TDP-Tg", "TDP-Tg", 
              "FUS-Tg", "FUS-Tg", "FUS-Tg", "FUS-Tg", "FUS-Tg", 
              "SOD-Control", "SOD-Control", "SOD-Control", 
              "FUS-Tg", "FUS-Tg", "FUS-Tg", "FUS-Tg",
              "SOD-Tg", "SOD-Tg", "SOD-Tg")

names(df) <- c("gene", "names", "group")

ggplot(df) + geom_bar(aes(x = names, y = gene, fill = group), stat = "identity") +
             theme(axis.text.x = element_text(angle = 90, hjust = 1))
              
