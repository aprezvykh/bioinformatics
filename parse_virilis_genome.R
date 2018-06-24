library(GenomicRanges)
rg <- rtracklayer::import('~/Documents/194.226.21.15/Rezvykh.stuff/RepMask/dvir-all-r1.06.gtf')
df <- rg[which(rg$type == "mRNA"),]
df <- as.data.frame(df)



for (f in unique(df$seqnames)){
  sub <- df[which(df$seqnames == f),]
  ss <- sub[which(sub$type == "mRNA"),]
    for (i in 1:nrow(ss)){
      a <- as.numeric(ss$start[i+1]) - as.numeric(ss$start[i])
      i <- i+1
      if(abs(a) > 200000){
        print(paste(f, "-", df$start, ":",df$end, sep = ""))
      }
    }
}
