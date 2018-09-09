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
  u <- unique(results_transcripts[grep(paste("\\b",x,"\\b",sep = ""), results_transcripts$geneNames),]$geneIDs)
  plotMeans(u, data, groupvar=gv, meas='FPKM', colorby='transcript',labelTranscripts = TRUE)
}
grep_and_plot_kmeans <- function(x,nk){
  u <- unique(results_transcripts[grep(paste("\\b",x,"\\b",sep = ""),results_transcripts$geneNames),]$geneIDs)
  plotLatentTranscripts(gene=as.character(u), gown=data, k=nk, method='kmeans', returncluster=FALSE) 
  
}
grep_and_plot_transcripts <- function(x,n,bb){
  u <- unique(results_transcripts[grep(paste("\\b",x,"\\b",sep = ""), results_transcripts$geneNames),]$geneIDs)
  plotTranscripts(u, data, meas = "FPKM",blackBorders = bb,samples = sampleNames(data)[n],labelTranscripts = T)
}
get_gene_transcripts <- function(x){
  df <- fpkm[grep(paste("\\b",x,"\\b",sep = ""), results_transcripts$geneNames),]
  return(df)
}
plot_fpkm <- function(id,meas){
  df <- data.frame(fpkm[id,])
  df$samples <- meas
  names(df) <- c("fpkm", "group")
  ggplot(data=df) + geom_boxplot(aes(x = group, y = fpkm, fill = group)) + 
    ggtitle(paste("Transcript ", id, sep = "")) +
    theme_bw()
}
plot_introns <- function(id,meas){
  df <- data.frame(introns[id,])
  df$samples <- meas
  names(df) <- c("fpkm", "group")
  ggplot(data=df) + geom_boxplot(aes(x = group, y = fpkm, fill = group)) + 
    ggtitle(paste("Transcript ", id, sep = "")) +
    theme_bw()
}
closestColor = function(x, colscale){
  choices = rev(heat.colors(length(colscale)))
  diffs = abs(x-colscale)
  return(choices[which.min(diffs)])
}

PlotTranscripts.mod <- function (gene, gown, samples = NULL, colorby = "transcript", 
                                 meas = "FPKM", legend = TRUE, labelTranscripts = FALSE, main = NULL, 
                                 blackBorders = TRUE, log = FALSE, logbase = 2, customCol = NULL, 
                                 customOrder = NULL) {
  if (class(gown) != "ballgown") 
    stop("gown must be a ballgown object")
  if (gown@RSEM & colorby == "exon") {
    stop(.makepretty("RSEM objects do not yet have exon-level measurements,\n            so must color by transcript."))
  }
  if (!gown@RSEM & meas == "TPM") {
    stop("only RSEM objects have TPM measurements.")
  }
  if (is.null(samples)) {
    samples = sampleNames(gown)[1]
    if (colorby != "none") 
      message(paste("defaulting to sample", samples))
  }
  if (!all(samples %in% sampleNames(gown))) {
    stop(.makepretty("all entries of \"samples\" must be part of \"gown\". Use\n            sampleNames(gown) to check."))
  }
  stopifnot(colorby %in% c("transcript", "exon", "none"))
  if (colorby == "transcript") {
    stopifnot(meas %in% c("cov", "FPKM", "TPM"))
  }
  if (colorby == "exon") {
    emeas = c("rcount", "ucount", "mrcount", "cov", "cov_sd", 
              "mcov", "mcov_sd")
    stopifnot(meas %in% emeas)
  }
  if (!(gown@meas == "all" | meas %in% gown@meas)) {
    stop(paste("gown does not contain", meas, "measurements."))
  }
  if (colorby == "none") 
    legend = FALSE
  if (!is.null(customCol) & (colorby != "transcript")) {
    stop("Custom coloring is only available at transcript level currently")
  }
  if (!is.null(customCol) & legend) {
    stop("legend must be FALSE if you provide custom colors")
  }
  n = length(samples)
  westval = ifelse(labelTranscripts, 4, 2)
  if (n > 1) {
    numrows = floor(sqrt(n))
    numcols = ceiling(n/numrows)
    par(mfrow = c(numrows, numcols), mar = c(5, westval, 
                                             4, 2), oma = c(0, 0, 2, 0))
  }
  else {
    par(mar = c(5, westval, 4, 2))
  }
  ma = IRanges::as.data.frame(structure(gown)$trans)
  if (names(ma)[2] != "group_name") {
    stop("IRanges::as.data.frame has changed. Please report as issue.")
  }
  thetranscripts = indexes(gown)$t2g$t_id[indexes(gown)$t2g$g_id == 
                                            gene]
  if (!is.null(customOrder)) {
    if (!all(sort(customOrder) == sort(thetranscripts))) {
      stop(.makepretty("customOrder must include each transcript in gene\n                exactly once."))
    }
  }
  if (substr(ma$group_name[1], 1, 2) == "tx") {
    warning("your ballgown object was built with a deprecated version of\n            ballgown - would probably be good to re-build!")
    thetranscripts = paste0("tx", thetranscripts)
  }
  if (!is.null(customCol) & (length(customCol) != length(thetranscripts))) {
    stop("You must have the same number of custom colors as transcripts")
  }
  gtrans = ma[ma$group_name %in% thetranscripts, ]
  xax = seq(min(gtrans$start), max(gtrans$end), by = 1)
  numtx = length(unique(thetranscripts))
  ymax = ifelse(legend, numtx + 1.5, numtx + 1)
  if (length(unique(gtrans$seqnames)) > 1) {
    stop("gene spans multiple chromosomes (?)")
  }
  if (length(unique(gtrans$strand)) > 1) {
    stop("gene contains exons from both strands (?)")
  }
  if (colorby != "none") {
    g_id = texpr(gown, "all")$gene_id
    if (colorby == "transcript") {
      smalldat = texpr(gown, meas)[which(g_id == gene), 
                                   ]
      t_id = texpr(gown, "all")$t_id[which(g_id == gene)]
    }
    if (colorby == "exon") {
      e_id_full = eexpr(gown, "all")$e_id
      smalldat = eexpr(gown, meas)[which(e_id_full %in% 
                                           gtrans$id), ]
      e_id = e_id_full[which(e_id_full %in% gtrans$id)]
    }
    if (numtx == 1) {
      snames = names(smalldat)
      smalldat = matrix(smalldat, nrow = 1)
      colnames(smalldat) = snames
    }
    if (log) {
      smalldat = log(smalldat + 1, base = logbase)
    }
    maxcol = quantile(as.matrix(smalldat), 0.99)
    colscale = seq(0, maxcol, length.out = 200)
    introntypes = c("ucount", "rcount", "mrcount")
    color.introns = meas %in% introntypes & gown@meas %in% 
      c("all", introntypes)
  }
  else {
    color.introns = TRUE
  }
  for (s in 1:n) {
    plot(xax, rep(0, length(xax)), ylim = c(0, ymax), type = "n", 
         xlab = "genomic position", yaxt = "n", ylab = "")
    if (n > 1) {
      title(samples[s])
    }
    colName = paste(meas, samples[s], sep = ".")
    if (colorby != "none") {
      colIndex = which(colnames(smalldat) == colName)
    }
    if (is.null(customOrder)) {
      transcript_loop = unique(gtrans$group_name)
    }
    else {
      transcript_loop = customOrder
    }
    for (tx in transcript_loop) {
      txind = which(transcript_loop == tx)
      if (colorby == "transcript") {
        if (!is.null(customCol)) {
          mycolor = customCol[txind]
        }
        else {
          mycolor = closestColor(smalldat[, colIndex][which(t_id == 
                                                              tx)], colscale)
        }
        stopifnot(length(mycolor) > 0)
      }
      else if (colorby == "none") {
        mycolor = "gray70"
      }
      gtsub = gtrans[gtrans$group_name == tx, ]
      gtsub = gtsub[order(gtsub$start), ]
      for (exind in 1:dim(gtsub)[1]) {
        if (colorby == "exon") {
          mycolor = closestColor(smalldat[, colIndex][which(e_id == 
                                                              gtsub$id[exind])], colscale)
          stopifnot(length(mycolor) > 0)
        }
        borderCol = ifelse(blackBorders, "black", mycolor)
        polygon(x = c(gtsub$start[exind], gtsub$start[exind], 
                      gtsub$end[exind], gtsub$end[exind]), y = c(txind - 
                                                                   0.4, txind + 0.4, txind + 0.4, txind - 0.4), 
                col = mycolor, border = borderCol)
        if (exind != dim(gtsub)[1]) {
          if (!color.introns) {
            lines(c(gtsub$end[exind], gtsub$start[exind + 
                                                    1]), c(txind, txind), lty = 2, col = "gray60")
          }
          if (color.introns) {
            intronindex = which(expr(gown)$intron$start == 
                                  gtsub$end[exind] + 1 & expr(gown)$intron$end == 
                                  gtsub$start[exind + 1] - 1 & expr(gown)$intron$chr == 
                                  unique(gtsub$seqnames) & expr(gown)$intron$strand == 
                                  unique(gtsub$strand))
            icolumnind = which(names(expr(gown)$intron) == 
                                 colName)
            icol = closestColor(expr(gown)$intron[intronindex, 
                                                  icolumnind], colscale)
            lines(c(gtsub$end[exind] + 10, gtsub$start[exind + 
                                                         1] - 10), c(txind, txind), lwd = 3, col = icol)
            lines(c(gtsub$end[exind], gtsub$start[exind + 
                                                    1]), c(txind + 0.1, txind + 0.1), lwd = 0.5, 
                  col = "gray60")
            lines(c(gtsub$end[exind], gtsub$start[exind + 
                                                    1]), c(txind - 0.1, txind - 0.1), lwd = 0.5, 
                  col = "gray60")
          }
        }
      }
    }
    if (legend) {
      leglocs = seq(min(xax) + 1, max(xax) - 1, length = length(colscale) + 
                      1)
      for (i in 1:length(colscale)) {
        polygon(x = c(leglocs[i], leglocs[i], leglocs[i + 
                                                        1], leglocs[i + 1]), y = c(ymax - 0.3, ymax, 
                                                                                   ymax, ymax - 0.3), border = NA, col = rev(heat.colors(length(colscale)))[i])
      }
      text(x = seq(min(xax) + 1, max(xax) - 1, length = 10), 
           y = rep(ymax + 0.1, 10), labels = round(colscale, 
                                                   2)[seq(1, length(colscale), length = 10)], 
           cex = 0.5)
      if (log) {
        text(x = median(xax), y = ymax - 0.5, labels = paste("log expression, by", 
                                                             colorby), cex = 0.5)
      }
      else {
        text(x = median(xax), y = ymax - 0.5, labels = paste("expression, by", 
                                                             colorby), cex = 0.5)
      }
    }
    if (labelTranscripts) {
      axis(side = 2, at = c(1:numtx), labels = transcript_loop, 
           cex.axis = 0.75, las = 1)
    }
  }
  if (!is.null(main)) {
    title(main, outer = (n > 1))
  }
  else {
    if (n == 1) {
      title(paste0(gene, ": ", samples))
    }
    else {
      title(gene, outer = TRUE)
    }
  }
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
introns <- iexpr(data, meas = "rcount")
whole_intron_table = iexpr(data, 'all')
transcript_id_by_exon = indexes(data)$e2t$t_id[match(unique(indexes(data)$e2t$e_id), indexes(data)$e2t$e_id)] 
transcript_id_by_gene = indexes(data)$t2g$t_id #transcript/gene mapping 
geneID = indexes(data)$t2g$g_id[match(transcript_id_by_exon, transcript_id_by_gene)] 

pheatmap(cor(introns[rowSums(introns)>0,]))
pheatmap(cor(exons[rowSums(exons)>0,]))


design <- model.matrix(~meta$groups) 
fit <- lmFit(eexpr(data), design = design, method = "ls", correlation = TRUE) 
ex <- diffSplice(fit, geneid = geneID) 
ts <- topSplice(ex,coef = ncol(fit$design),FDR = 0.05,number = 10000)
head(ts)

plotSplice(fit,FDR = 0.05,coef = ncol(fit$design),geneid = geneID,genecolname = 1)




grep_and_plot_fpkm("Hnrnph1", "transgenity")
grep_and_plot_transcripts("Hnrnph1",15,FALSE)
dev.off()
pheatmap(get_gene_transcripts("Hnrnph1"))
plot_fpkm(54735, meta$groups)
plot_introns(1483, meta$groups)
i <- get_gene_transcripts("Rpl35a")[,0]

grep_and_plot_transcripts("Ltbp4",11,FALSE)

plot_introns(124356, meta$groups)
plot_fpkm(124356, meta$groups)


p <- prcomp(t(introns), center = TRUE)
ggbiplot(p, groups = meta$groups, var.axes = F,ellipse = T) + theme_bw() + ggtitle("Introns,PCA")

p <- prcomp(t(fpkm), center = TRUE)
ggbiplot(p, groups = meta$groups, var.axes = F,ellipse = T) + theme_bw() + ggtitle("FPKM, PCA")
