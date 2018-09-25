col.pan <- colorpanel(100, "blue", "white", "red")
dir <- "~/transcriptomes/reads/intron_retention/all_counts/htseq/"
setwd(dir)
g <- grep("counts", list.files(dir), value = T)
files <- g
cond <- c("tdp-wt","tdp-wt","tdp-wt","tdp-wt",
          "tdp-tg", "tdp-tg", "tdp-tg", "tdp-tg", 
          "sod-wt", "sod-wt", "sod-tg", "sod-tg",
          "tg1", "tg1", "tg1", "tg1", "tg1", 
          "tg2", "tg2", "tg2", "tg2", 
          "tg3", "tg3", "tg3", "tg3", "tg3", 
          "wt1", "wt1", "wt1", "wt1", "wt1", 
          "wt3", "wt3", "wt3", "wt3", "wt3")

dge <- readDGE(files = files,path = "~/transcriptomes/reads/intron_retention/all_counts/htseq/", group = cond)
norm.factors <- calcNormFactors.default(dge,method = "TMM")
dgelist <- DGEList(dge, lib.size = colSums(dge$counts), norm.factors = norm.factors,samples = g,group = cond)
cpm <- cpm(dgelist)

cpm_sub <- cpm[order(rowSds(cpm), decreasing = T),][1:100,]
colnames(cpm_sub) <- cond
fh <- t(scale(t(cpm_sub)))
pheatmap(fh, col = col.pan,show_rownames = F)

cpm_sub <- cpm[order(rowSums(cpm), decreasing = T),][1:1000,]
colnames(cpm_sub) <- cond
fh <- t(scale(t(cpm_sub)))
pheatmap(fh, col = col.pan,show_rownames = F)


View(data$seg)
