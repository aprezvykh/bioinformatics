library(DEXSeq)
dir <- "~/transcriptomes/reads/intron_retention/merge/dexseq"
flat_gtf <- "~/transcriptomes/reads/intron_retention/merge/dexseq/flatten.gtf"
setwd(dir)
full <- TRUE
if(full){
    files <- grep("counts", list.files(dir),value = T)
    
    groups <- c("tdp-wt","tdp-wt","tdp-wt","tdp-wt",
                "tdp-tg", "tdp-tg", "tdp-tg", "tdp-tg", 
                "sod-wt", "sod-wt", "sod-tg", "sod-tg",
                "tg1", "tg1", "tg1", "tg1", "tg1", 
                "tg2", "tg2", "tg2", "tg2", 
                "tg3", "tg3", "tg3", "tg3", "tg3", 
                "wt1", "wt1", "wt1", "wt1", "wt1", 
                "wt3", "wt3", "wt3", "wt3", "wt3")
    
    tissue <- c("cortex", "cortex", "cortex", "cortex", 
                "cortex", "cortex", "cortex", "cortex", 
                "motoneurons", "motoneurons", 
                "motoneurons", "motoneurons", 
                "spinal", "spinal", "spinal", "spinal", "spinal", 
                "spinal", "spinal", "spinal", "spinal", 
                "spinal", "spinal", "spinal", "spinal", "spinal", 
                "spinal", "spinal", "spinal", "spinal", "spinal", 
                "spinal", "spinal", "spinal", "spinal", "spinal")
    
    model <- c("tdp", "tdp", "tdp", "tdp", 
               "tdp", "tdp", "tdp", "tdp", 
               "sod", "sod", "sod", "sod", 
               "fus", "fus", "fus", "fus", "fus", 
               "fus", "fus", "fus", "fus", 
               "fus", "fus", "fus", "fus", "fus", 
               "fus", "fus", "fus", "fus", "fus",
               "fus", "fus", "fus", "fus", "fus")
    
    transgenity <- c("wt", "wt", "wt", "wt", 
                     "tg", "tg", "tg", "tg", 
                     "wt", "wt", "tg", "tg", 
                     "tg", "tg", "tg", "tg", "tg", 
                     "tg", "tg", "tg", "tg", 
                     "tg", "tg", "tg", "tg", "tg", 
                     "wt", "wt", "wt", "wt", "wt", 
                     "wt", "wt", "wt", "wt", "wt")
    
    meta <- data.frame(files=files, groups = groups, tissue = tissue,
                       model = model, transgenity = transgenity)

}

setwd(dir)
#files <- grep("counts", list.files(dir), value = T)

#genes_subset <- read.csv("~/transcriptomes/reads/intron_retention/fus/cash/all_significant_events.csv")
#ref <- read.csv("~/transcriptomes/reads/intron_retention/fus/cash/all_significant_diff.csv")

#genes_subset <- genes_subset[genes_subset$st == "tg2tg3" & genes_subset$SplicingType == "IR",]

#genes_subset$ens <- mapIds(org.Mm.eg.db, 
#                    keys=as.character(genes_subset$AccID), 
#                    column="ENSEMBL", 
#                    keytype="SYMBOL",
#                    multiVals="first")

#g_sub <- as.character(genes_subset$ens)

mm <- model.matrix(~meta$groups)

dxd = DEXSeqDataSetFromHTSeq(
  files,
  sampleData=meta,
  design= ~ groups + exon + groups:exon,
  flattenedfile=flat_gtf)

#dxd = dxd[geneIDs(dxd) %in% g_sub,]

ncol(counts(dxd))
dxd = estimateSizeFactors(dxd)
dxd = estimateDispersions(dxd)

saveRDS(dxd, "~/transcriptomes/SPICA/data/dexseq/intronic_nonest_disp.rds")



formulaFullModel = ~ files + cond + cond:exon + files:exon
formulaReducedModel = ~ files + exon

dxd = testForDEU(dxd,
                 fullModel = formulaFullModel,
                 reducedModel = formulaReducedModel)

dxd = estimateExonFoldChanges(dxd,fitExpToVar="cond")
dxr = DEXSeqResults(dxd)

dxr$gene_name <-mapIds(org.Mm.eg.db, 
                                            keys=as.character(dxr$groupID), 
                                            column="SYMBOL", 
                                            keytype="ENSEMBL",
                                            multiVals="first")

View(genes_subset)

plotDEXSeq(dxr, "ENSMUSG00000025137", displayTranscripts=TRUE, legend=TRUE,
            cex.axis=1.3, cex=1.1, lwd=1,fitExpToVar = "cond")

dev.off()
df <- data.frame(dxr$exonBaseMean,
                 dxr$dispersion,
                 dxr$stat,
                 dxr$pvalue,
                 dxr$padj,
                 dxr$case,
                 dxr$control,
                 dxr$gene_name,stringsAsFactors = F,
                 dxr$log2fold_control_case)

names(df) <- c("bm", "disp","stat","pval","padj","case","control","gene_name","fc")
dev.off()

ggplot(data=df,aes(x = control, y = case)) + geom_point(col = ifelse(df$padj<0.05, "red", "black")) + 
  geom_text(aes(label = ifelse(df$padj<0.05, gene_name, "")),hjust=1, vjust=1,check_overlap = F)

mat <- dxr$countData
colnames(mat) <- cond
m <- mat[order(rowSds(mat), decreasing = T),][1:100,]
pheatmap(t(scale(t(m))))
