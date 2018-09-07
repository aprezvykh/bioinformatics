library("ballgown") 
library("dplyr") 
library("limma") 
library("org.Dm.eg.db") 
library(KEGG.db) 
dir <- "~/transcriptomes/reads/intron_retention/fus/ballgown/" 
setwd(dir) 

topSplice <- function(fit, coef=ncol(fit), test="simes", number=10L, FDR=1, sort.by="p"){
  #	Check fit is as produced by diffSplice
  if(is.null(fit$gene.genes$NExons)) stop("fit should be fit object produced by diffSplice")
  
  #	Can only specify one coefficient
  coef <- coef[1]
  
  test <- match.arg(test,choices=c("simes","F","f","t"))
  if(test=="f") test <- "F"
  
  sort.by <- match.arg(sort.by,choices=c("p","none","logFC","NExons"))
  if(sort.by=="logFC" & test!="t") stop("Sorting by logFC only available with Simes test")
  if(sort.by=="NExons" & test=="t") stop("Sorting by NExons only available with gene-level tests")
  
  #	Assemble data.frame of results for this coef
  switch(test,
         t = {
           out <- fit$genes
           out$logFC <- as.matrix(fit$coefficients)[,coef]
           out$t <- as.matrix(fit$t)[,coef]
           out$P.Value <- as.matrix(fit$p.value)[,coef]
         },
         
         F = {
           out <- fit$gene.genes
           out$F <- as.matrix(fit$gene.F)[,coef]
           out$P.Value <- as.matrix(fit$gene.F.p.value)[,coef]
         },
         
         simes = {
           out <- fit$gene.genes
           out$P.Value <- as.matrix(fit$gene.simes.p.value)[,coef]
         }
  )
  out$FDR <- p.adjust(out$P.Value, method="BH")
  
  #	Reduce to significant genes
  if(FDR<1) out <- out[out$FDR <= FDR,]
  
  #	Is the number of rows requested more than number available?
  number <- min(number, nrow(out))
  if(number <= 1L) return(out)
  
  #	Sort rows
  o <- switch(sort.by,
              p = order(out$P.Value, decreasing=FALSE),
              logFC = order(abs(out$logFC), decreasing=TRUE),
              NExons = order(out$NExons, -out$P.Value, decreasing=TRUE),
              none = 1:nrow(out)
  )
  o <- o[1:number]
  out[o,]
}
plotSplice <- function(fit, coef=ncol(fit), geneid=NULL, genecolname=NULL, rank=1L, FDR=0.05){
  if(is.null(genecolname)) 
    genecolname <- fit$genecolname
  else
    genecolname <- as.character(genecolname)
  
  if(is.null(geneid)) {
    #		Find gene from specified rank 
    if(rank==1L)
      i <- which.min(fit$gene.F.p.value[,coef])
    else
      i <- order(fit$gene.F.p.value[,coef])[rank]
    geneid <- paste(fit$gene.genes[i,genecolname], collapse=".")
  } else {
    #		Find gene from specified name
    geneid <- as.character(geneid)
    i <- which(fit$gene.genes[,genecolname]==geneid)[1]
    if(!length(i)) stop(paste("geneid",geneid,"not found"))
  }
  
  #	Row numbers containing exons
  j <- fit$gene.firstexon[i]:fit$gene.lastexon[i]
  
  exoncolname <- fit$exoncolname
  
  #	Get strand if possible		
  strcol <- grepl("strand", colnames(fit$gene.genes), ignore.case=TRUE)
  if(any(strcol)) geneid <- paste0(geneid, " (", as.character(fit$gene.genes[i, strcol])[1], ")")
  
  if(is.null(exoncolname)) {
    plot(fit$coefficients[j, coef], xlab="Exon", ylab="logFC (this exon vs rest)", main=geneid, type="b")
  } else {
    exon.id <- fit$genes[j, exoncolname]
    xlab <- paste("Exon", exoncolname, sep=" ")
    plot(fit$coefficients[j, coef], xlab="", ylab="logFC (this exon vs rest)", main=geneid, type="b", xaxt="n")
    axis(1, at=1:length(j), labels=exon.id, las=2, cex.axis=0.5)
    mtext(xlab, side=1, padj=5.2)
  }
  
  #	Mark the topSpliced exons
  top <- topSplice(fit, coef=coef, number=Inf, test="t", sort.by="none")
  m <- which(top[, genecolname] %in% fit$gene.genes[i, genecolname])
  fdr <- top$FDR[m]
  sig <- fdr < FDR
  if(any(sig)){
    fdr.sig <- fdr[sig]
    if(length(unique(fdr.sig))==1)
      cex <- 1.5
    else {
      abs.fdr <- abs(log10(fdr.sig))
      from <- range(abs.fdr)
      to <- c(1,2)
      cex <- (abs.fdr - from[1])/diff(from) * diff(to) + to[1]
    }
    points((1:length(j))[sig], fit$coefficients[j[sig], coef], col="red", pch=16, cex=cex)
  }
  
  abline(h=0,lty=2)
  invisible()
}
grep_and_plot_fpkm <- function(x,gv){
  u <- unique(results_transcripts[grep(x, results_transcripts$geneNames),]$geneIDs)
  plotMeans(u, data, groupvar=gv, meas='FPKM', colorby='transcript',labelTranscripts = TRUE)
}

grep_and_plot_kmeans <- function(x,nk){
  u <- unique(results_transcripts[grep(x, results_transcripts$geneNames),]$geneIDs)
  plotLatentTranscripts(gene=as.character(u), gown=data, k=nk, method='kmeans', returncluster=FALSE) 
  
}

grep_and_plot_transcripts <- function(x,n,bb){
  u <- unique(results_transcripts[grep(x, results_transcripts$geneNames),]$geneIDs)
  plotTranscripts(u, data, meas = "FPKM",blackBorders = bb,samples = sampleNames(data)[n],labelTranscripts = T)
}

plot_fpkm <- function(id,meas){
  df <- data.frame(fpkm[id,])
  df$samples <- meas
  names(df) <- c("fpkm", "group")
  ggplot(data=df) + geom_boxplot(aes(x = group, y = fpkm, fill = group)) + ggtitle(meas)
}


samples <- list.files(dir)

all <- rep("all", 24)

transgenity <- c("tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt")

age <- c("60d", "60d", "60d", "60d", "60d", 
         "90d", "90d", "90d", "90d", 
         "120d", "120d", "120d", "120d", "120d", 
         "60d", "60d", "60d", "60d", "60d", 
         "120d", "120d", "120d", "120d", "120d")

groups <- c('tg1', 'tg1', 'tg1', 'tg1', 'tg1', 
            'tg2', 'tg2', 'tg2', 'tg2', 
            'tg3', 'tg3', 'tg3', 'tg3', 'tg3', 
            'wt1', 'wt1', 'wt1', 'wt1', 'wt1', 
            'wt3', 'wt3', 'wt3', 'wt3', 'wt3')


meta <- list(samples=samples,transgenity=transgenity, age = age, groups = groups, all = all)

meta <- as.data.frame(meta)

data <- ballgown(dataDir = dir, samplePattern = "s", pData = meta,samples = list.files(dir)) 

results_transcripts = stattest(data, feature="transcript", getFC=TRUE, meas="FPKM", covariate = "transgenity") 
results_genes = stattest(data, feature="gene", getFC=TRUE, meas="FPKM", covariate = "transgenity") 

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


design <- model.matrix(~meta$groups) 
fit <- lmFit(eexpr(data), design = design, method = "ls", correlation = TRUE) 
ex <- diffSplice(fit, geneid = geneID) 
ts <- topSplice(ex,coef = ncol(fit$design),FDR = 0.05,number = 10000)
head(ts)

plotSplice(fit,FDR = 0.05,coef = ncol(fit$design),geneid = geneID,genecolname = 1)




