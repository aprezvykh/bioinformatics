library("ballgown") 
library("dplyr") 
library("limma") 
library("org.Dm.eg.db") 
library(KEGG.db) 
dir <- "~/ballgown.Dmel.memory/" 
setwd(dir) 

samples <- grep("trimmed", list.files(dir), value = TRUE) 
sampleCondition <- as.data.frame(samples) 
sampleCondition$condition <- c("F", "F", "mem", "mem", 
                               "control", "control", "N", "N") 
sampleCondition$cont <- c("0", "0", "0", "0", "1", "1", "0", "0") 
design <- model.matrix(~sampleCondition$condition) 

data <- ballgown(dataDir = dir, samplePattern = "trimmed", pData = sampleCondition) 

results_transcripts = stattest(data, feature="transcript", getFC=TRUE, meas="FPKM", covariate = "cont") 


results_genes = stattest(data, feature="gene", getFC=TRUE, meas="FPKM", covariate = "cont") 

results_transcripts = data.frame(geneNames=ballgown::geneNames(data), 
                                 geneIDs=ballgown::geneIDs(data), results_transcripts) 

indices <- match(results_genes$id, texpr(data, 'all')$gene_id) 
gene_names_for_result <- texpr(data, 'all')$gene_name[indices] 
results_genes <- data.frame(geneNames=gene_names_for_result, results_genes) 

fpkm <- texpr(data, meas = "FPKM") 
exons <- eexpr(data, meas = "rcount") 
transcript_id_by_exon = indexes(data)$e2t$t_id[match(unique(indexes(data)$e2t$e_id), indexes(data)$e2t$e_id)] 
transcript_id_by_gene = indexes(data)$t2g$t_id #transcript/gene mapping 
geneID = indexes(data)$t2g$g_id[match(transcript_id_by_exon, transcript_id_by_gene)] 

intersect(results_genes$id, ts$GeneID) 


design <- model.matrix(~sampleCondition$condition) 
fit <- lmFit(eexpr(data), design = design, method = "ls", correlation = TRUE) 
ex <- diffSplice(fit, geneid = geneID) 
ts <- topSplice(ex, number = 40000) 
ts <- subset(ts, ts$FDR < 0.05) 
ts <- subset(ts, ts$P.Value < 0.05) 

plotTranscripts(data, gene = "MSTRG.8272", samples = sampleCondition$samples) 
plotSplice(ex, geneid = "MSTRG.11277") 

plotLatentTranscripts(gene='MSTRG.11277', gown=data, k=3, method='kmeans', returncluster=FALSE) 
dev.off()