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
    library(org.Ce.eg.db)
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
    library(knitr)
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
    library(PoissonSeq)
    library(ggplot2)
    library(yaml)
    library(tibble)
    
    
    ###BIOCBARALLEL SETTINGS
  
    default <- registered()
    register(BatchJobsParam(workers = 10), default = TRUE)
    options(MulticoreParam=quote(MulticoreParam(workers=8)))
    param <- SnowParam(workers = 2, type = "SOCK")
    
    
    heatmaps <- TRUE
    custom_genes_plots <- FALSE
    analyze_all_samples <- FALSE
    disease_association <- FALSE
    kegg_plots <- TRUE
    panther_analysis <- TRUE
    deseq2_part <- TRUE
    qlm_test <- FALSE
    logging <- FALSE
    boxplots <- TRUE
    biotype <- FALSE
    distrib <- FALSE
    summary <- FALSE
      ### CONSTANTS BLOCK
      
    
    pvalue_cutoff <- 0.05
    logfchigh_cutoff <- 0
    logfclow_cutoff <- -0
    cpm_cutoff <- -1
    gs_size <- 1
    diseases_set <- 50
    number_of_kegg_plots <- 100
    go_terms_set <- 50
    pathways_set <- 30
    genes_in_term <- 3
    filter_thresh <- 5
    baseMean_cutoff <- 1.5
    go_heatmap_count <- 20
    stattest_number <- 1
    fdr_cutoff <- 0.05
    
  directory <- '~/counts/AIKAR.permutation.test.late/'
  setwd(directory)
  gr_control <- c("control_late")
  gr_case <- c("aikar_late")
      
  ### BUILDING A SPECIFIC DESIGN TABLE
  if (logging == TRUE){
     zz <- file("error.log", open="wt")
        sink(zz, type="message")
      }
  col.pan <- colorpanel(100, "blue", "white", "red")
  ###DIRECTORY WHERE SAMPLES ARE LOCATED
  library(dplyr)
  if (analyze_all_samples == TRUE){
              sampleFiles <- grep('worm',list.files(directory),value=TRUE)
              sampleFiles
              sampleCondition <- c('control_early',
                                   'aikar_late',
                                   'control_late',
                                   'aikar_late',
                                   'aikar_early',
                                   'control_early',
                                   'aikar_early',
                                   'control_early',
                                   'aikar_early',
                                   'control_late',
                                   'aikar_late',
                                   'control_late'
                                  )
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
        #cpm$rowsum <- rowSums(cpm) 
        keep <- rowSums(cpm > cpm_cutoff) >= ncol(sampleTable) 
        logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
        logCPM <- logCPM[!(row.names(logCPM) %in% row.names.remove), ]
        logCPM <- logCPM[keep,]
        nrow(a$counts)
        a <- a[keep, , keep.lib.sizes=FALSE] 
        a <- calcNormFactors(a, method = "TMM") 
        design <- model.matrix(~sampleTable$condition) 
        a <- estimateDisp(a,design) 
        fit <- glmQLFit(a,design, robust = TRUE) 
        qlf <- glmQLFTest(fit,coef=ncol(fit$design))
        et_annot <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BY"))
        
        et_annot_non_filtered <- as.data.frame(topTags(qlf, n = nrow(logCPM), adjust.method = "BH"))
        top <- as.data.frame(topTags(qlf, n = 20, adjust.method = "BH", sort.by = "PValue"))
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
        #y <- y[keep, ,keep.lib.sizes=FALSE] 
        logCPM <- as.data.frame(cpm(y, log = TRUE, lib.size = colSums(counts) * normalized_lib_sizes))
        et <- exactTest(y) 
        top <- as.data.frame(topTags(et)) 
        et_annot <- as.data.frame(topTags(et, n = nrow(logCPM), adjust.method = "BH", sort.by = "PValue"))
        et_annot_non_filtered <- as.data.frame(et$table) 
  }
  
  
      correl <- logCPM
      correl <- as.data.frame(correl)
      correl <- correl[complete.cases(correl),]
      x <- cor(correl, method = "pearson")
      pdf(file = "Corellation matrix Spearman.pdf", width = 10, height = 10)
      heatmap.2(x, col=col.pan, Rowv=TRUE, scale="none",
                trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
                margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")
      
      dev.off()
    
    
      
      
      pch <- c(0,1,2,15,16,17)
    colors <- rep(c("darkgreen", "red", "blue"), 2)
    # pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
    pdf(file = "MDSPlot.pdf", width = 10, height = 10)
    plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition], labels = sampleTable$condition)
    legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
    dev.off()
    

    
    
    
    
    
### DENSITY DISTRIBUTION
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


if (distrib == TRUE){

if(qlm_test == TRUE){

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
} else {
  lib <- data.frame(y$samples$lib.size, sampleTable$condition)
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
}


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


log <- data.frame(y$AveLogCPM, et_annot_non_filtered$logFC)
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


}



y$genes$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")
y$genes$Name <- mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot_non_filtered), 
                         column="GENENAME", 
                         keytype="WORMBASE",
                         multiVals="first")

et_annot$symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")

et_annot$name <- mapIds(org.Ce.eg.db, 
                        keys=row.names(et_annot), 
                        column="GENENAME", 
                        keytype="WORMBASE",
                        multiVals="first")


et_annot$entrez <- mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot), 
                         column="ENTREZID", 
                         keytype="WORMBASE",
                         multiVals="first")

et_annot$GOID <-     mapIds(org.Ce.eg.db, 
                         keys=row.names(et_annot), 
                         column="GO", 
                         keytype="WORMBASE",
                         multiVals="first")

et_annot$term <- mapIds(GO.db, 
                       keys=et_annot$GOID, 
                       column="TERM", 
                       keytype="GOID",
                       multiVals="first")

et_annot$uniprot <-mapIds(org.Ce.eg.db, 
                          keys=row.names(et_annot), 
                          column="UNIPROT", 
                          keytype="WORMBASE",
                          multiVals="first")

et_annot$term <- as.character(et_annot$term)

et_annot_non_filtered$Symbol <- mapIds(org.Ce.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="SYMBOL", 
                          keytype="WORMBASE",
                          multiVals="first")
et_annot_non_filtered$entrez <- mapIds(org.Ce.eg.db, 
                          keys=row.names(et_annot_non_filtered), 
                          column="ENTREZID", 
                          keytype="WORMBASE",
                          multiVals="first")


top$Symbol <- mapIds(org.Ce.eg.db, 
                         keys=row.names(top), 
                         column="SYMBOL", 
                         keytype="WORMBASE",
                         multiVals="first")
top$Name<- mapIds(org.Ce.eg.db, 
                       keys=row.names(top), 
                       column="GENENAME", 
                       keytype="WORMBASE",
                       multiVals="first")


top <- top[complete.cases(top),]
### CPM
cpm <- as.data.frame(cpm(y))

cpm$Symbol <- mapIds(org.Ce.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="WORMBASE",
                     multiVals="first")

cpm$Name <- mapIds(org.Ce.eg.db, 
                   keys=row.names(cpm), 
                   column="GENENAME", 
                   keytype="WORMBASE",
                   multiVals="first")

cpm$entrez <- mapIds(org.Ce.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="WORMBASE",
                     multiVals="first")

cpm$GOID <-     mapIds(org.Ce.eg.db, 
                       keys=row.names(cpm), 
                       column="GO", 
                       keytype="WORMBASE",
                       multiVals="first")

cpm$term <- mapIds(GO.db, 
                   keys=cpm$GOID, 
                   column="TERM", 
                   keytype="GOID",
                   multiVals="first")
cpm$term <- as.character(cpm$term)

write.csv(cpm, "cpm.csv")

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
taxon = 'Caenorhabditis Elegans'
taxon = tolower(taxon)
tmp = unlist(strsplit(x = taxon, split = ' '))
dataset.name = tolower(sprintf('%s%s_gene_ensembl', substr(tmp[1],1,1), tmp[2]))
mart <- useMart("ensembl", dataset=dataset.name) #, host="www.ensembl.org"
needed.attributes = c("ensembl_gene_id","external_gene_name", "description","gene_biotype","entrezgene", "transcript_length")

gmt_flt = getBM(attributes=needed.attributes,filters="ensembl_gene_id",values=rownames(et_annot), mart=mart)
gmt_flt = gmt_flt[!(duplicated(gmt_flt[,"ensembl_gene_id"])),]
rownames(gmt_flt) = gmt_flt[,"ensembl_gene_id"]

###transcript length-logFC

trlen.common <- intersect(rownames(et_annot.trlen), gmt_flt$ensembl_gene_id)
et_annot.trlen.common <- et_annot.trlen[(rownames(et_annot.trlen) %in% trlen.common),]
gmt_flt.common <- gmt_flt[(gmt_flt$ensembl_gene_id %in% trlen.common),]


et_annot.trlen.common <- et_annot.trlen.common[order(rownames(et_annot.trlen.common), decreasing = TRUE),]
gmt_flt.common <- gmt_flt.common[order(gmt_flt.common$ensembl_gene_id, decreasing = TRUE),]


trlen.distrib <- data.frame(et_annot.trlen.common$logFC, gmt_flt.common$transcript_length, et_annot.trlen.common$logCPM)
names(trlen.distrib) <- c("fc", "tl", "cpm")
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

png("Transcript length - LogFC distribution.png", height = 20, width = 20)
plot(fc ~ tl, data = trlen.distrib, pch = ".", xlab = "Log2(Transcript Length)", ylab = "Log2(Fold Change)", main = "Transcript Length vs LogFC")
abline(fit)
dev.off()




### FILTRATION

et_annot <- as.data.frame(subset(et_annot, logCPM > cpm_cutoff))
et_annot <- as.data.frame(subset(et_annot, PValue < pvalue_cutoff))

et_annot <- as.data.frame(subset(et_annot, FDR < pvalue_cutoff))


et_annot <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff | logFC < logfclow_cutoff))
et_annot <- et_annot[complete.cases(et_annot), ]




### BIOTYPE ANNOT FIX IT!
if (biotype == TRUE){
    taxon = 'Drosophila melanogaster'
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

# WRITE RESULTS
write.xlsx(df, file = "Results edgeR.xlsx", sheetName = "Simple Summary", append = TRUE)
write.xlsx(top, file = "Results edgeR.xlsx", sheetName = "Top Tags (with FDR)", append = TRUE)
write.xlsx(et_annot, file = "Results edgeR.xlsx", sheetName = "Filtered Genes, logCPM, logfc", append = TRUE)

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
CountsTable$symbol <- mapIds(org.Ce.eg.db, 
                             keys=row.names(CountsTable), 
                             column="SYMBOL", 
                             keytype="WORMBASE",
                             multiVals="first")

CountsTable$name <- mapIds(org.Ce.eg.db, 
                           keys=row.names(CountsTable), 
                           column="GENENAME", 
                           keytype="WORMBASE",
                           multiVals="first")


### GO TESTS
foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez


### GOANA
et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))

goana_up <- goana(de = et_annot_high$entrez, species = "Ce")
go_up_30 <- topGO(goana_up, n=30)
go_up_30$perc = (go_up_30$DE/go_up_30$N)*100
go_up_100 <- topGO(goana_up, n=100)
go_up_100$perc = (go_up_100$DE/go_up_100$N)*100
go_up_500 <- topGO(goana_up, n=500)
go_up_500$perc = (go_up_500$DE/go_up_500$N)*100
go_up_500$perc <- round(go_up_500$perc, digits = 4)

goana_down <- goana(de = et_annot_low$entrez, species = "Ce")
go_down_30 <- topGO(goana_down, n=30)
go_down_30$perc = (go_down_30$DE/go_down_30$N)*100
go_down_100 <- topGO(goana_down, n=100)
go_down_100$perc = (go_down_100$DE/go_down_100$N)*100
go_down_500 <- topGO(goana_down, n=500)
go_down_500$perc = (go_down_500$DE/go_down_500$N)*100
go_down_500$perc <- round(go_down_500$perc, digits = 4)


write.xlsx(go_up_30, file = "Goana GO tests, upreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_up_100, file = "Goana GO tests, upreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_up_500, file = "Goana GO tests, upreg.xlsx", sheetName = "top500", append = TRUE)

write.xlsx(go_down_30, file = "Goana GO tests, downreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_down_100, file = "Goana GO tests, downreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_down_500, file = "Goana GO tests, downreg.xlsx", sheetName = "top500", append = TRUE)


write.xlsx(go_up_100, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated GO terms", append = TRUE)
write.xlsx(go_down_100, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated GO terms", append = TRUE)

###ANNOTATE LOGCPM FOR GO HEATMAPS
for_go_heatmap <- logCPM

for_go_heatmap$entrez <- mapIds(org.Ce.eg.db, 
                                keys=row.names(for_go_heatmap), 
                                column="ENTREZID", 
                                keytype="WORMBASE",
                                multiVals="first")

for_go_heatmap$symbol <- mapIds(org.Ce.eg.db, 
                                keys=row.names(for_go_heatmap), 
                                column="SYMBOL", 
                                keytype="WORMBASE",
                                multiVals="first")



### GO HEATMAPS OF VARIOUS GENES

all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez
GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)

allGO <- genesInTerm(GOdata)

##UP
remove <- c("MF", "CC")
dir.create("GO heatmap plots upregulated")
setwd("GO heatmap plots upregulated")

go_up_subset <- subset(go_up_500, go_up_500$DE > 40)
go_up_subset <- subset(go_up_subset, go_up_subset$DE < 60)
go_up_subset <- subset(go_up_subset, go_up_subset$N < 400)
go_up_subset <- go_up_subset[!(go_up_subset$Ont %in% remove), ] 
up_terms <- data.frame(go_up_subset$Term)
rownames(up_terms) <- rownames(go_up_subset)


for (f in seq(1:nrow(go_up_subset))){
      z <- allGO[paste(rownames(up_terms))[f]]
      z <- as.data.frame(z)
      goid <- names(z)
      goname <- go_up_subset[f,1]
      goperc <- go_up_subset[f,6]
      gopvalue <- go_up_subset[f,5]
      statperc <- paste("Percent of involved genes", goperc, sep = ":")
      name <-paste(goid, goname, sep = ":")
      name <- as.character(name)
      name <- gsub('[.]', "", name)
      name <- paste(name, statperc, sep = "\n")
      names(z) <- c("sus")
      mat <- for_go_heatmap[(for_go_heatmap$entrez %in% z$sus),]
      rownames(mat) <- mat$symbol
      mat$symbol <- NULL
      mat$entrez <- NULL
      mat <- t(scale(t(mat)))
      pdf(file = paste(goid,"upreg_GO.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
      heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none",
                trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
                margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = paste(name))
      dev.off()
      
}

##DOWN

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

dir.create("GO heatmap plots downregulated")
setwd("GO heatmap plots downregulated")
go_down_subset <- subset(go_down_500, go_down_500$DE > 40)
go_down_subset <- subset(go_down_subset, go_down_subset$DE < 60)
go_down_subset <- subset(go_down_subset, go_down_subset$N < 400)
go_down_subset <- go_down_subset[!(go_down_subset$Ont %in% remove), ] 
down_terms <- data.frame(go_down_subset$Term)
rownames(down_terms) <- rownames(go_down_subset)

for (f in seq(1:nrow(go_down_subset))){
  z <- allGO[paste(rownames(down_terms))[f]]
  z <- as.data.frame(z)
  goid <- names(z)
  goname <- go_down_subset[f,1]
  goperc <- go_down_subset[f,6]
  gopvalue <- go_down_subset[f,5]
  statperc <- paste("Percent of involved genes", goperc, sep = ":")
  name <-paste(goid, goname, sep = ":")
  name <- as.character(name)
  name <- gsub('[.]', "", name)
  name <- paste(name, statperc, sep = "\n")
  names(z) <- c("sus")
  mat <- for_go_heatmap[(for_go_heatmap$entrez %in% z$sus),]
  rownames(mat) <- mat$symbol
  mat$symbol <- NULL
  mat$entrez <- NULL
  mat <- t(scale(t(mat)))
  pdf(file = paste(goid,"downreg_GO.pdf",sep=""), width = 12, height = 17, family = "Helvetica")
  heatmap.2(mat, col=col.pan, Rowv=TRUE, scale="none",
            trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = paste(name))
  dev.off()
  
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

s <- subset(go_up_500, go_up_500$perc > 20 & go_up_500$P.DE < 0.001)
s <- s[which(s$Ont == "BP"),]
s <- s[which(nchar(s$Term) < 50),]
s <- s[order(s$DE, decreasing = FALSE),]
s <- s[seq(1:30),]
s <- s[complete.cases(s),]
png(file = "Go terms upreg.png", width = 1024, height = 768)
#pdf(file = "Go terms upreg.pdf", width = 12, height = 17, family = "Helvetica")
g_u <- ggplot(s, aes(x = reorder(Term, DE), y = DE, fill = P.DE)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(breaks = s$Term, name = "Significant Upregulated Terms") + 
  coord_flip() +
  scale_y_continuous(name = "Number of differentialy expressed genes") +
  theme_bw() + 
  theme(axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"))
g_u
dev.off()

s <- subset(go_down_500, go_down_500$perc > 20 & go_down_500$P.DE < 0.001)
s <- s[which(s$Ont == "BP"),]
s <- s[which(nchar(s$Term) < 50),]
s <- s[order(s$DE, decreasing = FALSE),]
s <- s[seq(1:30),]
s <- s[complete.cases(s),]
png(file = "Go terms downreg.png", width = 1024, height = 768)
#pdf(file = "Go terms downreg.pdf", width = 12, height = 17, family = "Helvetica")
g_d <- ggplot(s, aes(x = reorder(Term, DE), y = DE, fill = P.DE)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(breaks = s$Term, name = "Significant Downregulated Terms") + 
  coord_flip() +
  scale_y_continuous(name = "Number of differentialy expressed genes") +
  theme_bw() + 
  theme(axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"))
g_d
dev.off()

getwd()

### GO with foldchanges
go_fc <- goana(de = et_annot$entrez, species = "Ce")
go_fc_res <- topGO(go_fc, n=1000, ontology = "BP")
df <- data.frame()
for (f in rownames(go_fc_res)){
  x <- allGO[[paste(f)]]
  z <- et_annot[(et_annot$entrez %in% x),]
  mean.logfc<- mean(z$logFC)
  mean.cpm <- mean(z$logCPM)
  b <- go_fc_res[which(rownames(go_fc_res) == f),]
  b$logFC <- mean.logfc
  b$cpm <- mean.cpm
  df <- rbind(b, df)
}
df <- df[complete.cases(df),]
write.xlsx(df, file = "goana tests with foldchanges.xlsx", sheetName = "Top 1000 GO terms")

### FISHER GO TEST WITH ELIMINATION
### WARNING! THIS TEST MIGHT BE USELESS!
all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez

GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 30", append = TRUE)
write.xlsx(top100.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(top500.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 500", append = TRUE)


dir.create("GO graphs")
setwd("GO graphs")

pdf(file = "top 5 GO graph.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


#MF

GOdata = new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 30", append = TRUE)
write.xlsx(top100.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(top500.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 500", append = TRUE)


setwd("GO graphs")

pdf(file = "top 5 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}



##CC


GOdata = new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 30", append = TRUE)
write.xlsx(top100.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 100", append = TRUE)
write.xlsx(top500.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 500", append = TRUE)


setwd("GO graphs")

pdf(file = "top 5 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()


setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}




## CORELLATION MARIX
correl <- logCPM
correl <- as.data.frame(correl)
correl <- correl[complete.cases(correl),]
x <- cor(correl, method = "pearson")
pdf(file = "Corellation matrix Spearman.pdf", width = 10, height = 10)
heatmap.2(x, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1.4, cexCol=1.4, density.info="none",
          margin=c(18,18), lhei=c(2,10), lwid=c(2,6), main = "Spearman corellation")

dev.off()
#legend("bottomright",
#       legend = c("Control-1", "Control-3", "Tg-1", "Tg-2", "Tg-3"),
#       col = c("yellow", "purple", "green", "red", "blue"),
#       lty= 1,
#       lwd = 10)

#heatmap.2(x, margins = c(10,10))      
### MDS PLOT
getwd()

pch <- c(0,1,2,15,16,17)
colors <- rep(c("darkgreen", "red", "blue"), 2)
# pdf(file = "PCAPlot.pdf", width = 12, height = 17, family = "Helvetica")
pdf(file = "MDSPlot.pdf", width = 10, height = 10)
plotMDS(y, col=colors[sampleTable$condition], pch = pch[sampleTable$condition], labels = sampleTable$condition)
legend("topleft", legend=levels(sampleTable$condition), pch=pch, col=colors, ncol=2)
dev.off()



### VOLCANO PLOT
allgenes <- nrow(et_annot_non_filtered)
max <- max(-log10(et_annot_non_filtered$PValue))
et_annot_non_filtered$threshold = as.factor(abs(et_annot_non_filtered$logFC) > logfchigh_cutoff & et_annot_non_filtered$PValue < 0.05/allgenes)
pdf(file = "Volcano plot.pdf", width = 10, height = 10)
volc = ggplot(data=et_annot_non_filtered, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=1, size=2) +
  labs(legend.position = "none") +
  xlim(c(-6, 6)) + ylim(c(1.30103, max)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw()
volc
dev.off()

#SMEAR PLOT
if (qlm_test == TRUE){
  tt <- topTags(qlf, n=nrow(y), p.value=0.05)
} else {
  tt <- topTags(et, n=nrow(y), p.value=0.05)
  
}
pdf(file = "Smear plot.pdf", width = 10, height = 10)
if (qlm_test == TRUE){
    plotSmear(qlf, de.tags = rownames(tt$table), lowess = TRUE, smooth.scatter = TRUE)
} else {
    plotSmear(et, de.tags = rownames(tt$table), lowess = TRUE, smooth.scatter = TRUE)
  
}
dev.off()

pdf(file = "BCV plot.pdf", width = 10, height = 10)
plotBCV(y = a, xlab = "Average LogCPM", ylab = "Biological coefficient of variation")
dev.off()


pdf(file = "Quasi-likelhood dispersion.pdf", width = 10, height = 10)
plotQLDisp(fit)
dev.off()


pdf(file = "Barcode plot-logFC.pdf", width = 10, height = 10)
barcodeplot(statistics = et_annot$PValue, gene.weights = et_annot$logFC, weights.label = "Log2 Fold Change", worm = TRUE, xlab = "P-value")
dev.off()


pdf(file = "Barcode plot-logCPM.pdf", width = 10, height = 10)
barcodeplot(statistics = et_annot$logCPM, gene.weights = et_annot$logFC, weights.label = "Log2 Fold Change", worm = TRUE, xlab = "logCPM")
dev.off()



### REACTOME PART
dfa <- as.character(et_annot$entrez)
x_com <- enrichPathway(gene = dfa, pAdjustMethod = "BH", universe = et_annot$entrez, organism = "celegans", minGSSize=gs_size, readable = TRUE)
write.xlsx(x_com, "Reactome.xlsx", sheetName = "All reactome", append = TRUE)
if (nrow(x) == 0){
  for(i in seq(1:2)){
    beep()
    Sys.sleep(0.1)
    beep()
    Sys.sleep(0.2)
    print("No terms found!")  
  }
}
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot.pdf", width = 12, height = 17, family = "Helvetica")
r.bp.com <- barplot(x_com, showCategory=50,  font.size = 9)
r.bp.com
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
r.em.com <- enrichMap(x_com, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
r.cn.com <- cnetplot(x_com, foldChange = foldchanges, categorySize="pvalue", showCategory = 10)
dev.off()

#HIGH

df_high <- et_annot_high$entrez
x_up <- enrichPathway(gene=df_high, organism = "celegans", minGSSize=gs_size, readable = TRUE, pAdjustMethod = "BH")
write.xlsx(x_up, "Reactome.xlsx", sheetName = "High", append = TRUE)
write.xlsx(x_up, file = "Results edgeR.xlsx", sheetName = "top 100 Upregulated Reactome pathways", append = TRUE)
head(as.data.frame(x_up))
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
r.bp.up <- barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x_down <- enrichPathway(gene=df_low, organism = "celegans", minGSSize=gs_size, readable = TRUE )
write.xlsx(x_down, "Reactome.xlsx", sheetName = "Low", append = TRUE)
write.xlsx(x_down, file = "Results edgeR.xlsx", sheetName = "top 100 Downregulated Reactome pathways", append = TRUE)
head(as.data.frame(x_down))

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
r.bp.down <- barplot(x_down, showCategory=30,  font.size = 9)
dev.off()


###KEGG expression profile (without lfc, but with generatio)
df_kk_high <- et_annot_high$uniprot
kk <- enrichKEGG(gene = df_kk_high, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

df_kk_low <- et_annot_low$uniprot
kk <- enrichKEGG(gene = df_kk_low, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=30,  font.size = 9)
dev.off()

###KEGGA
library(limma)
keg_com <- kegga(de = et_annot$entrez, species="Ce", convert = TRUE)
tk_common <- topKEGG(keg_com, n=100)
tk_common <- subset(tk_common, tk_common$P.DE < 0.05)
tk_common$perc <- (tk_common$DE/tk_common$N)*100
tk_common <- tk_common[order(tk_common$perc, decreasing = TRUE),]
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Ce", convert = TRUE)
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Ce", convert = TRUE)
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
getwd()


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
topcpm$Symbol <- mapIds(org.Ce.eg.db, 
                        keys=row.names(topcpm), 
                        column="SYMBOL", 
                        keytype="WORMBASE",
                        multiVals="first")
rownames(topcpm) <- make.names(topcpm$Symbol, unique=TRUE)
topcpm$Symbol <- NULL
top100cpm <- t(scale(t(topcpm)))
pdf(file = "All heatmap.pdf", width = 12, height = 17, family = "Helvetica")
heatmap.2(top100cpm, col=col.pan, Rowv=TRUE, scale="column",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "")
            
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


library(curl)
plot_pathway = function(pid){
         pathview(gene.data=foldchanges, 
         pathway.id=pid, 
         species="cel", 
         new.signature=FALSE)
}


dir.create("another kegg plots")
setwd("another kegg plots")
getwd()

tk_common <- tk_common[which(rownames(tk_common) != "cel04215"),]
tk_common

for (f in rownames(tk_common)){
  plot_pathway(f)
}

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}


p <- prcomp(t(y$counts), center = TRUE, scale. = TRUE)

ggbiplot(p, choices = c(3,4))


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






r <- grep("WBGene00006575", rownames(cpm), ignore.case = TRUE)
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


r <- grep("WBGene00000913", rownames(cpm), ignore.case = TRUE)
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
  geom_signif(comparisons = list(c("Tg-1", "Tg-2")), map_signif_level = TRUE) + 
  geom_signif(comparisons = list(c("Control-1", "Tg-1")), map_signif_level = TRUE) + 
  
  
  #geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
  ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
  theme(plot.title = element_text(hjust = 0.5))


g

