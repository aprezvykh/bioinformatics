library(GenomicRanges)
rg <- rtracklayer::import('~/Documents/194.226.21.15/Rezvykh.stuff/RepMask/dvir-all-r1.06.gtf')
#df <- rg[which(rg$type == "mRNA"),]
df <- as.data.frame(rg)

for (f in unique(df$seqnames)){
  sub <- df[which(df$seqnames == f),]
  ss <- sub[which(sub$type == "start_codon"),]
    for (i in 1:nrow(ss)){
      a <- as.numeric(ss$start[i+1]) - as.numeric(ss$start[i])
      i <- i+1
      if(is.na(a) == TRUE){
        next
      } else if (is.null(a) == TRUE){
        next
      } else if (length(a) < 1){
        next
      }
      if(abs(as.numeric(a)) > 200000 | abs(as.numeric(a) < 300000)){
      print(paste(f, "-", ss$start[i+1], ":",ss$end[i], sep = ""))
      } else {
        next
      }
    }
}


