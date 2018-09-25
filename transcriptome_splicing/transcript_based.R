library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(Rsubread)
library(annotables)

dir <- "~/transcriptomes/reads/intron_retention/all_counts/fc_transcripts/"
setwd(dir)
g <- grep("parsed.counts", list.files(dir),value = T)

samples <- g
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

meta <- data.frame(samples = samples, groups = groups, tissue = tissue,
             model = model, transgenity = transgenity)

#g <- g[grep("fus", meta$model)]

#meta <- meta[grep("fus", meta$model),]

meta$samples <- as.character(meta$samples)
meta$groups <- as.character(meta$groups)
meta$tissue <- as.character(meta$tissue)
meta$model <- as.character(meta$model)
meta$transgenity <- as.character(meta$transgenity)

design <- model.matrix(~meta$groups)
design
dge <- readDGE(g,path = dir,columns = c(1,7), group = meta$groups)
dge <- calcNormFactors(dge, method = "TMM")
#dge <- DGEList(dge,lib.size = colSums(dge$counts))
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
cpm <- cpm(dge,normalized.lib.sizes = T,log = F)

cpm$enstrans <- rownames(cpm)
cpm <- as.data.frame(cpm)
cpm <- cpm %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

saveRDS(cpm,"~/transcriptomes/SPICA/data/ensmust.rds")


#fit <- glmFit(dge, design)
#glm <- glmLRT(fit, design)
#tt <- data.frame(topTags(glm,adjust.method = "BH", n = 100))
n <- nrow(dge$counts)
et <- exactTest(dge, pair = c("wt1","wt3"))
tt <- data.frame(topTags(et,n = n))

tt$enstrans <- rownames(tt)

tt <- tt %>% 
  dplyr::arrange(FDR) %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

tt <- tt[tt$FDR<0.05,]
tt$symbol <- mapIds(org.Mm.eg.db, 
                     keys=as.character(tt$ensgene), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

tt$name <- mapIds(org.Mm.eg.db, 
                   keys=as.character(tt$ensgene), 
                   column="GENENAME", 
                   keytype="ENSEMBL",
                   multiVals="first")

tt$stattest <- "wt1wt3"
tab <- rbind(tt, tab)

write.xlsx(tab, "/mnt/raid/illumina/AlexR/diffexp_isoforms.xlsx", sheetName = "iso")

frq <- data.frame(table(t(tt$ensgene)))
frq <- frq[frq$Freq>1,]
mult.trans.genes <- as.character(frq$Var1)

for(f in mult.trans.genes){
  print(f)
  print(tt[tt$ensgene == f,]$logFC)
}

cpm <- as.data.frame(cpm(dge,log = F))
cpm$enstrans <- rownames(cpm)
View(tt)

cpm <- cpm %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

sub <- cpm[which(cpm$ensgene == "ENSMUSG00000037509"),]
rownames(sub) <- sub$enstrans
sub$enstrans <- NULL
sub$ensgene <- NULL

pheatmap(sub)
png("~/transcriptomes/reads/intron_retention/arhgef_trans_expr.png")
pheatmap(sub[,which(meta$model == "fus")])
dev.off()

#int.gene <- c("Arhgef4")
#ENSMUST00000159747 - up
#ENSMUST00000047664 - down