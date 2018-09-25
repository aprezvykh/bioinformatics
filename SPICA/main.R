options(warn = F)
library(ballgown)
library(DEXSeq)
library(SAJR)
library(gplots)
library(DescTools)
library("msa")
library(igraph)
library(regsplice)
grep_and_plot_fpkm <- function(x,gv){
  u <- unique(data$ballgownres$res.trans[grep(paste("\\b",x,"\\b",sep = ""), data$ballgownres$res.trans$geneNames),]$geneIDs)
  plotMeans(u, data$ballgown, groupvar=gv, meas='FPKM', colorby='transcript',labelTranscripts = TRUE)
}
closestColor = function(x, colscale){
  choices = rev(heat.colors(length(colscale)))
  diffs = abs(x-colscale)
  return(choices[which.min(diffs)])
}

plotTranscripts.mod <- function (gene, gown, samples = NULL, colorby = "transcript", 
                                 meas = "FPKM", legend = TRUE, labelTranscripts = FALSE, main = NULL, 
                                 blackBorders = TRUE, log = FALSE, logbase = 2, customCol = NULL, 
                                 customOrder = NULL, lim.x) {
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
    color.introns = FALSE
  }
  for (s in 1:n) {
    plot(xax, rep(0, length(xax)), ylim = c(0, ymax), type = "n", 
         xlab = "genomic position", yaxt = "n", ylab = "",xlim = lim.x)
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
col.pan <- colorpanel(100, "blue", "white", "red")
dir <- "~/transcriptomes/SPICA/"
setwd(dir)
data <- readRDS("data/all_merged_data.rds")

data$cash$ens <- mapIds(org.Mm.eg.db, 
                                        keys=as.character(data$cash$AccID), 
                                        column="ENSEMBL", 
                                        keytype="SYMBOL",
                                        multiVals="first")
#############

gene <- "Arhgef4"
ens_id <- unique(data$cash[grep(paste("\\b", gene, "\\b", sep = ""), data$cash$AccID),]$ens)

hts_expr <- data$htseq$counts[grep(ens_id, rownames(data$htseq$counts)),]
sajr_expr <- data$sajr.expr$cnts[grep(ens_id, rownames(data$sajr.expr$cnts)),]
sajr_ir <- data$sajr.rawdata$ir[grep(ens_id, rownames(data$sajr.rawdata$ir)),]

g <- data$ballgownres$res.trans[grep(paste("\\b", gene, "\\b", sep = ""),data$ballgownres$res.trans$geneNames),]$id

sajr_ir[is.na(sajr_ir)] <- 0
pheatmap(sajr_ir)

exons_expr <- data$ballgownres$exons[g,]
introns_expr <- data$ballgownres$introns[g,]

dexseq_exon <- data$dexseq_exon_fc[grep(ens_id, rownames(data$dexseq_exon_fc)),]
dexseq_introns_only <- data$dexseq_intron_fc[grep("i", rownames(data$dexseq_intron_fc)),]
dexseq_intron <- dexseq_introns_only[grep(ens_id, rownames(dexseq_introns_only)),]

pheatmap(dexseq_intron)
###dexseq_sub part
dx_sub_e = data$dexseq_exon[geneIDs(data$dexseq_exon) %in% ens_id,]
dx_sub_i = data$dexseq_introm[geneIDs(data$dexseq_introm) %in% ens_id,]

dx_sub_e <- estimateSizeFactors(dx_sub_e)
dx_sub_e <- estimateDispersions(dx_sub_e)

formulaFullModel = ~ files + transgenity + transgenity:exon + files:exon
formulaReducedModel = ~ files + exon

dx_sub_e = testForDEU(dx_sub_e,
                      fullModel = formulaFullModel,
                      reducedModel = formulaReducedModel)

dx_sub_e = estimateExonFoldChanges(dx_sub_e,fitExpToVar="transgenity")
dxr_e = DEXSeqResults(dx_sub_e)
gene_transcripts <- unique(unlist(dxr_e$transcripts, recursive = T))

dir.create("tmp")

gtf_gene <- data$gtf[grep(ens_id, data$gtf$gene_id),]
gtf_gene <- gtf_gene[grep('\\bgene\\b', gtf_gene$type),]
system(paste("touch ", 
             ens_id, 
             ".gene.fasta", 
             sep = ""))

system(paste("echo ", "'", ">", ens_id, "' ", ">> ",ens_id, ".gene.fasta", sep = ""))

system(paste('samtools faidx data/ref/refgenome.fasta ', 
      as.character(gtf_gene$seqnames), 
      ":", as.character(gtf_gene$start), 
      "-", as.character(gtf_gene$end),
      " | sed 1d", 
      " >> ",
      ens_id, 
      ".gene.fasta",
      sep = ""
      ))

for(f in gene_transcripts){
    print(f)
    transcript_sub <- data$gtf[grep(f, data$gtf$transcript_id),]
    transcript_sub <- transcript_sub[grep("\\btranscript", transcript_sub$type),]
    tr <- paste(transcript_sub$seqnames,":",transcript_sub$start,"-",transcript_sub$end, sep = "")
    system(paste("touch ", ens_id, ".transcripts.fasta", sep = ""))
    system(paste("echo ", "'", ">", f, "' ", ">> ",ens_id, ".transcripts.fasta", sep = ""))
    system(paste("samtools faidx data/ref/refgenome.fasta ", tr, "| sed 1d"," >> ", ens_id, ".transcripts.fasta",sep = ""))
}

gene_fasta <- grep("gene", list.files(dir), value = T)
trans_fasta <- grep("trans", list.files(dir), value = T)
system(paste("makeblastdb -in ", trans_fasta, " -dbtype nucl", sep = ""))

system(paste("blastn -db ", 
      trans_fasta, 
      " -query ", 
      gene_fasta, 
      " -num_threads 16",
#      " -qcov_hsp_perc 100",
      " -outfmt 6", 
      " >> ",
      ens_id,
      ".blasted",
      sep = ""))


blast.path <- paste(ens_id, ".blasted", sep = "")

df <- read.delim(blast.path)

names(df) <- c("qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen", 
                "qstart", 
                "qend", 
                "sstart", 
                "send", 
                "evalue", 
                "bitscore")

df <- df[df$pident == 100,]
df <- df[df$bitscore > 100,]


comb <- t(combn(nrow(df),2))
comb <- data.frame(comb)

non_hom_genes <- vector()
for(f in 1:nrow(comb)){
  o1 <- seq(df$qstart[as.numeric(comb[f,][1])], df$qend[as.numeric(comb[f,][2])])
  o2 <- seq(df$sstart[as.numeric(comb[f,][1])], df$send[as.numeric(comb[f,][2])])
    if(length(intersect(o1,o2))<1){
      non_hom_genes <- append(paste(df$sseqid[as.numeric(comb[f,][1])], "-", df$sseqid[as.numeric(comb[f,][2])], sep = ""), non_hom_genes)
      print(length(intersect(o1,o2)))
    }
}
non_hom_genes

df$idx <- seq(1:nrow(df))

ggplot(data=df, aes(x = qstart, y = idx)) + geom_segment(aes(x = qstart, xend = qend, y = idx, yend = idx)) + 
  geom_text(aes(label = sseqid),hjust = 1) + theme_bw()

system("rm *.fasta *.blasted *.nhr *.nin *.nsq")
dev.off()


plotDEXSeq(dxr_e, ens_id, displayTranscripts=TRUE, legend=TRUE,
           cex.axis=1.3, cex=1.1, lwd=1,fitExpToVar = "transgenity")



#int.gene <- c("Arhgef4")
#ENSMUST00000159747 - up
#ENSMUST00000047664 - down
v1 <- data$gtf[grep("ENSMUST00000159747", data$gtf$transcript_id),]
v2 <- data$gtf[grep("ENSMUST00000047664", data$gtf$transcript_id),]

ggplot(data=v1) + geom_segment(aes(x = start, xend = end, y = type, yend = type, color = type))
ggplot(data=v2) + geom_segment(aes(x = start, xend = end, y = type, yend = type, color = type))
dev.off()

ids <- data$ballgownres$res.trans[grep("\\bArhgef4\\b", data$ballgownres$res.trans$geneNames),]$id
ids <- as.character(ids)

pheatmap(data$ballgownres$fpkm[ids,data$meta$model == "fus"])

plotTranscripts("MSTRG.330", 
                gown = data$ballgown, 
                samples = "s_tg2-3",
                blackBorders = F,
                labelTranscripts = T)


######


write.csv(data$cash,"~/transcriptomes/reads/intron_retention/all_cash.csv")
