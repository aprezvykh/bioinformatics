#!/usr/bin/Rscript
library(annotables)
library(stringr)
library(seqinr)
library(ggbiplot)
library(edgeR)
library(limma)
estimate.transcript.abundance <- function(i){
  sums <- rowSums(subcpm[subcpm$ensgene == i,][,1:samples.in.set])
  vars <- as.character(rowVars(subcpm[subcpm$ensgene == i,][,1:samples.in.set]))
  vv <- vector()
  for(f in sums){
    v <- 100*(f/sum(sums))
    # print(100*(f/sum(sums)))
    vv <- append(v,vv)
  }
  vv <- rev(vv)
  vv <- data.frame(vv)
  vv$trans <- subcpm[subcpm$ensgene == i,]$enstrans
  vv$gene <- i
  vv$vars <- vars
  names(vv) <- c("abund", "trans", "gene", "var")
  return(vv)
}


parse_coords <- function(df){
  s <- str_split_fixed(df$Location, ":", 2)
  chr <- s[,1]
  start <- str_split_fixed(s[,2],"-",2)[,1]
  stop <- str_split_fixed(s[,2],"-",2)[,2]
  df$chr <- chr
  df$start <- as.numeric(start)
  df$stop <- as.numeric(stop)
  df$div <- df$stop - df$start
  return(df)
}
cutoff_sign <- function(df){
  zz <- df[which(df$FDR<0.05),]
  zz <- zz[which(abs(zz$delta_PSI)>0.1),]
  return(zz)
}



#get.overlap.transcripts <- function(){
#  setwd(dir)
#  system("rm -r tmp/")
#  dir.create("tmp")
#  setwd("tmp")
#  ens_id <- "ENSMUSG00000025006"
#  gtf_gene <- gtf[grep(ens_id, gtf$gene_id),]
#  gtf_gene <- gtf_gene[grep('\\bgene\\b', gtf_gene$type),]
#  gene_transcripts <- spliced.iso[spliced.iso$ensgene == ens_id,]$enstrans
#  transcript.coordinates <- vector()
  
  
#  for(f in gene_transcripts){
#    print(f)
#    transcript_sub <- gtf[grep(f, gtf$transcript_id),]
#    transcript_sub <- transcript_sub[grep("\\btranscript", transcript_sub$type),]
#    tr <- paste(transcript_sub$seqnames,":",transcript_sub$start,"-",transcript_sub$end, sep = "")
#    transcript.coordinates <- append(tr, transcript.coordinates)
#    system(paste("touch ", ens_id, ".transcripts.fasta", sep = ""))
#    system(paste("echo ", "'", ">", f, "' ", ">> ",ens_id, ".transcripts.fasta", sep = ""))
#    system(paste("samtools faidx ~/transcriptomes/reads/intron_retention/fus/ref/refgenome.fa ", tr, "| sed 1d"," >> ", ens_id, ".transcripts.fasta",sep = ""))
#  }
  
#  names(transcript.coordinates) <- gene_transcripts
  
#  trans_fasta <- grep("trans", list.files(getwd()), value = T)
#  system(paste("makeblastdb -in ", trans_fasta, " -dbtype nucl", sep = ""))
  
#  system(paste("blastn -db ", 
#               trans_fasta, 
#               " -query ", 
#               trans_fasta, 
#               " -num_threads 16",
#               #      " -qcov_hsp_perc 100",
#               " -outfmt 6", 
#               " >> ",
#               ens_id,
#               ".blasted",
#               sep = ""))
  
  
#  blast.path <- paste(ens_id, ".blasted", sep = "")
  
#  df <- read.delim(blast.path)
  
#  names(df) <- c("qseqid",
#                 "sseqid",
#                 "pident",
#                 "length",
#                 "mismatch",
#                 "gapopen", 
#                 "qstart", 
#                 "qend", 
#                 "sstart", 
#                 "send", 
#                 "evalue", 
#                 "bitscore")
  
#  df <- df[df$pident == 100,]
#  df <- df[df$bitscore > 100,]
  
#  non.overlap.blast <- df[which(df$qseqid != df$sseqid),]
#  non.overlap.blast <- non.overlap.blast[non.overlap.blast$pident == 100,]

#  library(seqinr)
#  trans.fasta <- read.fasta("ENSMUSG00000025006.transcripts.fasta")
#  splice.event.location <- as.character(cash[grep(ens_id, cash$ens),]$Location)
#  gene.gtf.features <- gtf[grep(ens_id, gtf$gene_id),]
  
#  trans.overlap <- vector()
  
#  for(i in 1:nrow(gene.gtf.features)){
#    sq <- seq(gene.gtf.features$start[i], gene.gtf.features$end[i])
#    if(!identical(gtf[grep(40344356, sq),]$transcript_id, character(0))){
#      print(gtf[grep(40344356, sq),]$transcript_id)
#      trans.overlap <- append(gtf[grep(40344356, sq),]$transcript_id, trans.overlap)
#    }
#  }
  
#  trans.overlap <- trans.overlap[complete.cases(trans.overlap)]
#  intersect(unique(trans.overlap), res.abund[grep(ens_id, res.abund$ensgene),]$enstrans)
  
#}

###setting directory
dir <- "~/transcriptomes/reads/intron_retention/sod_glia/"
setwd(dir)

###reading and preparing source data
iso <- read.delim("rsem/Run_2018-10-02_12-32-49.isoforms.rsem.quantification.tsv")
genes <- read.delim("rsem/Run_2018-10-02_12-32-49.genes.rsem.quantification.tsv")


iso$SRR646651.rsem.isoforms.results <- NULL
iso$SRR646671.rsem.isoforms.results <- NULL

genes$SRR646651.rsem.genes.results <- NULL
genes$SRR646671.rsem.genes.results <- NULL

cash <- read.csv("cash.csv")
cash <- cash[,!grepl("^[tg]*[wt]", names(cash))]

gtf <- as.data.frame(rtracklayer::import("~/transcriptomes/ref/gff/Mus_musculus.GRCm38.90.gtf"))

iso$X <- as.character(iso$X)
rownames(iso) <- iso$X
iso$X <- NULL

#########ATTTTENCION!!!!!!!
###group_cnt and group_case varibales should be included
###in samples data frame, names must be the same. 

group_cnt <- "wt1"
group_case <- "tg1"

samples <- data.frame(samples = colnames(iso),
                      condition = c("tg1", "tg1", 
                                    "tg2", "tg2", "tg2", 
                                    "tg3", "tg3", "tg3",
                                    "wt1", "wt1"))

design <- model.matrix(~samples$condition)

dge <- DGEList(counts = iso,
               samples = samples$samples,
               group = samples$condition)

dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
cpm <- data.frame(cpm(dge,log = F))

tt <- exactTest(dge,pair = c(group_cnt, group_case))

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

write.csv(res, paste(group_cnt, "_", group_case, ".res.csv", sep  = ""))

##isoform abundance
grep.cpm.regexp <- paste(group_cnt, "|", group_case, "|", "rsem.ids", sep = "")
cpm <- as.data.frame(cpm(dge))
rsem.ids <- rownames(cpm)
rownames(cpm) <- seq(1:nrow(cpm))

cpm <- cpm[,grep(grep.cpm.regexp, samples$condition)]
samples.in.set <- ncol(cpm)
cpm$rsem.id <- rsem.ids
cpm <- data.frame(separate(data = cpm, col = rsem.id, into = c("enstrans", "gene.name"), sep = "_"))
cpm <- cpm %>% 
  dplyr::inner_join(grcm38_tx2gene, by = c("enstrans" = "enstxp"))

subcpm <- cpm

tr.ab <- NULL
tr.ab <- lapply(unique(subcpm$ensgene), estimate.transcript.abundance)

dfb <- data.frame()

for(f in tr.ab){
  df <- data.frame(f)
  dfb <- rbind(df, dfb)
  print(nrow(dfb))
}



#dfb[is.na(dfb),]$abund <- 0
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

cash.stattest <- paste(group_cnt, group_case, sep = "")
mul.iso.genes
mig <- res[res$ensgene %in% mul.iso.genes,]

ii <- intersect(mig$ensgene, cash[cash$st == cash.stattest,]$ens)


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


splice.mig.genes <- intersect(mig.abund$ensgene, cash$ens)

spliced.iso <- mig.abund[mig.abund$ensgene %in% splice.mig.genes,]


##########################
cash
spliced.cash <- cash[cash$ens %in% res$ensgene,]
big.common.table <- data.frame()
spliced.cash
spliced.cash$ens <- as.character(spliced.cash$ens)
scash <- NULL
scash <- cash[cash$st == cash.stattest,]
spliced.scash <- spliced.cash[spliced.cash$st == cash.stattest,]
spliced.iso <- res[res$ensgene %in% spliced.scash$ens,]

#nrow(spliced.scash)
unique(spliced.scash$ens)


for(f in unique(spliced.scash$ens)){
    ens_id <- f
    print(ens_id)
    splice.event.location <- as.character(scash[grep(ens_id, scash$ens),]$Location)
        for(k in splice.event.location){
        print(k)
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
        spl.start <- strsplit(strsplit(k, ":")[[1]][2], "-")[[1]][1]
        spl.end <- strsplit(strsplit(k, ":")[[1]][2], "-")[[1]][2]
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
        if(nrow(spliced.iso[spliced.iso$enstrans %in% tr.with.spl.events,]) == 0){
              next
        }
        common.table <- spliced.iso[spliced.iso$enstrans %in% tr.with.spl.events,]
        common.table$dpsi <- spliced.scash[grep(k, spliced.scash$Location),]$delta_PSI
        common.table$type <- spliced.scash[grep(k, spliced.scash$Location),]$SplicingType
        common.table$div <- spliced.scash[grep(k, spliced.scash$Location),]$div
        common.table$spl.FDR <- spliced.scash[grep(k, spliced.scash$Location),]$FDR
        common.table$spl.loc <- spliced.scash[grep(k, spliced.scash$Location),]$Location
        common.table$nex <- spliced.scash[grep(k, spliced.scash$Location),]$Exon
        big.common.table <- rbind(common.table, big.common.table)
                }
}

big.common.table$micro <- ifelse(big.common.table$div<28, "micro", "normal")

big.common.table.abund <- data.frame()

for(f in big.common.table$enstrans){
  print(f)
  ab <- dfb[grep(f, dfb$trans),]$abund
  sas <- big.common.table[grep(f, big.common.table$enstrans),]
  sas$abund <- ab
  big.common.table.abund <- rbind(sas, big.common.table.abund)
}

res.filename <- paste(group_cnt, "_", group_case, ".corr.csv", sep = "")
write.csv(big.common.table.abund, paste("~/transcriptomes/reads/intron_retention/sod_glia/", res.filename, sep = ""))

View(big.common.table.abund)

