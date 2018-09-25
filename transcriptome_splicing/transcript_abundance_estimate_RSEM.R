library(annotables)
library(stringr)
library(seqinr)
estimate.transcript.abundance <- function(i){
  sums <- rowSums(subcpm[subcpm$ensgene == i,][,1:samples.in.set])
  vars <- as.character(rowVars(subcpm[subcpm$ensgene == i,][,1:samples.in.set]))
  vv <- vector()
  for(f in sums){
    v <- 100*(f/sum(sums))
    # print(100*(f/sum(sums)))
    vv <- append(v,vv)
  }
  vv <- data.frame(vv)
  vv$trans <- subcpm[subcpm$ensgene == i,]$enstrans
  vv$gene <- i
  vv$vars <- vars
  names(vv) <- c("abund", "trans", "gene", "var")
  return(vv)
}

col.pan <- colorpanel(100, "blue", "white", "red")
dir <- "~/transcriptomes/reads/intron_retention/fus/"
setwd("~/transcriptomes/reads/intron_retention/fus/")
iso <- read.delim("Run_2018-09-20_11-18-06.isoforms.rsem.quantification.tsv")
genes <- read.delim("Run_2018-09-20_11-18-06.genes.rsem.quantification.tsv")
iso$X <- as.character(iso$X)
rownames(iso) <- iso$X
iso$X <- NULL
cash <- read.csv("~/transcriptomes/reads/intron_retention/all_cash.csv")
gtf <- as.data.frame(rtracklayer::import("~/transcriptomes/ref/gff/Mus_musculus.GRCm38.90.gtf"))

samples <- data.frame(samples = colnames(iso),
                      condition = c("wt1", "wt1", "wt1", "wt1", "wt1", 
                                    "wt3", "wt3", "wt3", "wt3", "wt3", 
                                    "tg1", "tg1", "tg1", "tg1", "tg1", 
                                    "tg2", "tg2", "tg2", "tg2", 
                                    "tg3", "tg3", "tg3", "tg3", "tg3"))

design <- model.matrix(~samples$condition)

dge <- DGEList(counts = iso,
               samples = samples$samples,
               group = samples$condition)

dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
cpm <- data.frame(cpm(dge,log = F))

tt <- exactTest(dge,pair = c("tg2", "tg3"))

res <- topTags(tt, 
               n = nrow(dge$counts), 
               adjust.method = "BH")

res <- as.data.frame(res)
res <- res[which(res$FDR<0.05),]

res$rsem.id <- rownames(res)
rownames(res) <- seq(1:nrow(res))

res <- data.frame(separate(data = res, col = rsem.id, into = c("enstrans", "gene.name"), sep = "_"))
res <- res %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

res$symbol <- mapIds(org.Mm.eg.db, 
                     keys=as.character(res$ensgene), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

res$name <- mapIds(org.Mm.eg.db, 
                   keys=as.character(res$ensgene), 
                   column="GENENAME", 
                   keytype="ENSEMBL",
                   multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db, 
       keys=as.character(res$ensgene), 
       column="ENTREZID", 
       keytype="ENSEMBL",
       multiVals="first")


##isoform abundance
cpm <- as.data.frame(cpm(dge))
rsem.ids <- rownames(cpm)
rownames(cpm) <- seq(1:nrow(cpm))

cpm <- cpm[,grep("tg2|tg3|rsem.id", samples$condition)]
samples.in.set <- ncol(cpm)
cpm$rsem.id <- rsem.ids
cpm <- data.frame(separate(data = cpm, col = rsem.id, into = c("enstrans", "gene.name"), sep = "_"))
cpm <- cpm %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

df <- data.frame(table(t(cpm$ensgene)))
df <- df[df$Freq>1,]
subcpm <- cpm[cpm$ensgene %in% df$Var1,]
subcpm <- subcpm[rowSums(subcpm[,1:samples.in.set])>10,]

tr.ab <- NULL
tr.ab <- lapply(unique(subcpm$ensgene), estimate.transcript.abundance)

dfb <- data.frame()

for(f in tr.ab){
  df <- data.frame(f)
  dfb <- rbind(df, dfb)
}

dfb[is.na(dfb),]$abund <- 0
###dbg
dfb <- dfb[order(dfb$trans),]
dfb$abund <- as.numeric(dfb$abund)
dfb$var<- as.numeric(dfb$var)
dfb$sign.var <- dfb$var/dfb$abund

int <- intersect(res$enstrans, dfb$trans)
res <- res[order(res$enstrans),]

dfb.abund <- dfb[dfb$trans %in% int,]
res.abund <- res[res$enstrans %in% int,]
res.abund$abund <- dfb.abund$abund
res.abund$sign.var <- dfb.abund$sign.var
View(res.abund)


fq <- data.frame(table(t(res$ensgene)))
fq <- fq[fq$Freq > 1,]
fq$Var1 <- as.character(fq$Var1)


mul.iso.genes <- vector()
for(f in fq$Var1){
  df <- res[res$ensgene == f,]
  lfcs <- as.character(df$logFC)
  lfcs <- as.numeric(lfcs)
  if(length(lfcs[lfcs<0]) > 0 & length(lfcs[lfcs>0]) > 0){
    print(f)
    mul.iso.genes <- append(f, mul.iso.genes)
  }
}

mul.iso.genes
mig <- res[res$ensgene %in% mul.iso.genes,]

ii <- intersect(mig$ensgene, cash[cash$st == "tg2tg3",]$ens)


dfb <- dfb[order(dfb$trans),]
dfb$abund <- as.numeric(dfb$abund)
dfb$var<- as.numeric(dfb$var)
dfb$sign.var <- dfb$var/dfb$abund

int <- intersect(mig$enstrans, dfb$trans)
mig <- mig[order(mig$enstrans),]

dfb.abund <- dfb[dfb$trans %in% int,]
mig.abund <- mig[mig$enstrans %in% int,]
mig.abund$abund <- dfb.abund$abund
mig.abund$sign.var <- dfb.abund$sign.var

View(mig.abund)
View(mig.abund[mig.abund$abund>20,])

splice.mig.genes <- intersect(mig.abund$ensgene, cash$ens)

spliced.iso <- mig.abund[mig.abund$ensgene %in% splice.mig.genes,]
spliced.cash <- cash[cash$ens %in% splice.mig.genes,]

View(spliced.iso)


##########################
setwd(dir)
system("rm -r tmp/")
dir.create("tmp")
setwd("tmp")
ens_id <- "ENSMUSG00000025006"
gtf_gene <- gtf[grep(ens_id, gtf$gene_id),]
gtf_gene <- gtf_gene[grep('\\bgene\\b', gtf_gene$type),]
gene_transcripts <- spliced.iso[spliced.iso$ensgene == ens_id,]$enstrans
transcript.coordinates <- vector()


for(f in gene_transcripts){
  print(f)
  transcript_sub <- gtf[grep(f, gtf$transcript_id),]
  transcript_sub <- transcript_sub[grep("\\btranscript", transcript_sub$type),]
  tr <- paste(transcript_sub$seqnames,":",transcript_sub$start,"-",transcript_sub$end, sep = "")
  transcript.coordinates <- append(tr, transcript.coordinates)
  system(paste("touch ", ens_id, ".transcripts.fasta", sep = ""))
  system(paste("echo ", "'", ">", f, "' ", ">> ",ens_id, ".transcripts.fasta", sep = ""))
  system(paste("samtools faidx ~/transcriptomes/reads/intron_retention/fus/ref/refgenome.fa ", tr, "| sed 1d"," >> ", ens_id, ".transcripts.fasta",sep = ""))
}

names(transcript.coordinates) <- gene_transcripts

trans_fasta <- grep("trans", list.files(getwd()), value = T)
system(paste("makeblastdb -in ", trans_fasta, " -dbtype nucl", sep = ""))

system(paste("blastn -db ", 
             trans_fasta, 
             " -query ", 
             trans_fasta, 
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

non.overlap.blast <- df[which(df$qseqid != df$sseqid),]
non.overlap.blast <- non.overlap.blast[non.overlap.blast$pident == 100,]
View(non.overlap.blast)

library(seqinr)
trans.fasta <- read.fasta("ENSMUSG00000025006.transcripts.fasta")
splice.event.location <- as.character(cash[grep(ens_id, cash$ens),]$Location)
gene.gtf.features <- gtf[grep(ens_id, gtf$gene_id),]

trans.overlap <- vector()

for(i in 1:nrow(gene.gtf.features)){
    sq <- seq(gene.gtf.features$start[i], gene.gtf.features$end[i])
    if(!identical(gtf[grep(40344356, sq),]$transcript_id, character(0))){
      print(gtf[grep(40344356, sq),]$transcript_id)
      trans.overlap <- append(gtf[grep(40344356, sq),]$transcript_id, trans.overlap)
    }
}

trans.overlap <- trans.overlap[complete.cases(trans.overlap)]
intersect(unique(trans.overlap), res.abund[grep(ens_id, res.abund$ensgene),]$enstrans)



####
View(spliced.cash)
big.common.table <- data.frame()
spliced.cash$ens <- as.character(spliced.cash$ens)

scash <- cash[cash$st == "tg2tg3",]
spliced.scash <- spliced.cash[spliced.cash$st == "tg2tg3",]
unique(spliced.scash$ens)[2]

for(f in unique(spliced.scash$ens)){
    ens_id <- f
    print(ens_id)
    splice.event.location <- as.character(scash[grep(ens_id, scash$ens),]$Location)
    gene.gtf.features <- gtf[grep(ens_id, gtf$gene_id),]
    gene_transcripts <- spliced.iso[spliced.iso$ensgene == ens_id,]$enstrans
    transcript.coordinates <- vector()
    for(i in gene_transcripts){
      print(i)
      transcript_sub <- gtf[grep(i, gtf$transcript_id),]
      transcript_sub <- transcript_sub[grep("\\btranscript", transcript_sub$type),]
      tr <- paste(transcript_sub$seqnames,":",transcript_sub$start,"-",transcript_sub$end, sep = "")
      transcript.coordinates <- append(tr, transcript.coordinates)
    }
    names(transcript.coordinates) <- gene_transcripts
    spl.start <- strsplit(strsplit(splice.event.location, ":")[[1]][2], "-")[[1]][1]
    spl.end <- strsplit(strsplit(splice.event.location, ":")[[1]][2], "-")[[1]][2]
    
    transcript.coordinates <- data.frame(transcript.coordinates)
    transcript.coordinates <- separate(transcript.coordinates,col = 1,into = c("stop", "start", by = "-"))
    names(transcript.coordinates) <- c("chr", "start", "stop")
    tr.with.spl.events <- vector()
    for(p in 1:nrow(transcript.coordinates)){
      i.var <- intersect(seq(spl.start, spl.end), seq(transcript.coordinates$start[p], transcript.coordinates$stop[p]))
      print(p)
      print(i.var)
        if(!identical(i.var, integer(0))){
        tr.with.spl.events <- append(rownames(transcript.coordinates)[p], tr.with.spl.events)
      } 
    }
    common.table <- spliced.iso[spliced.iso$enstrans %in% tr.with.spl.events,]
    common.table$Loc <- spliced.scash[spliced.scash$ens %in% ens_id,]$Location
    common.table$dPSI <- spliced.scash[spliced.scash$ens %in% ens_id,]$delta_PSI
    common.table$FDR <- spliced.scash[spliced.scash$ens %in% ens_id,]$FDR
    common.table$type <- spliced.scash[spliced.scash$ens %in% ens_id,]$SplicingType
    common.table$div <- spliced.scash[spliced.scash$ens %in% ens_id,]$div
    big.common.table <- rbind(common.table, big.common.table)

}
View(big.common.table)



