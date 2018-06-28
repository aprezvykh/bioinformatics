library(GenomicRanges)
rg <- rtracklayer::import('~/Documents/194.226.21.15/Rezvykh.stuff/RepMask/dvir-all-r1.06.gtf')
#df <- rg[which(rg$type == "mRNA"),]
df <- as.data.frame(rg)

z <- vector()
for (f in unique(df$seqnames)){
  sub <- df[which(df$seqnames == f),]
  ss <- sub[which(sub$type == "start_codon"),]
    for (i in 1:nrow(ss)){
     a <- as.numeric(ss$start[i+1]) - as.numeric(ss$start[i])
     i <- i+1
     z <- append(z,
                 ifelse(a>200000, paste(f, ":",ss$start[i+1],"-",ss$start[i],sep = ""), 
                                        print("SAS")))
  }
}


z <- sub("SAS", "", z)
z <- z[!z == ""]

z <- z[!grepl("NA", z)]
z <-z[complete.cases(z)]
z


scaffold_12875:14091566-14439860
