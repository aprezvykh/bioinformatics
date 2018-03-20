
library(seqbias)
library(Rsamtools)
library(ggplot2)

dir <- c("~/R/x86_64-pc-linux-gnu-library/3.2/seqbias/aikar/")
setwd(dir)
l <- grep("bam", list.files(dir), value = TRUE)
ref_fn <- system.file("aikar/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa", package = "seqbias" )
ref_f <- FaFile(ref_fn)
ref_seqs <- scanFaIndex(ref_f)
I <- random.intervals(ref_seqs, n = 5, m = 100000)
seqs <- scanFa(ref_f, I)
neg_idx <- as.logical(I@strand == '-')
seqs[neg_idx] <- reverseComplement(seqs[neg_idx])


for (f in l){
  print(f)
  reads_fn <- system.file(paste("aikar/", f, sep = ""), package = "seqbias")
  counts <- count.reads(reads_fn, I, binary = T)
  
  
  freqs <- kmer.freq(seqs, counts)
  if( require(ggplot2) ) {
    P <- qplot( x = pos,
                y = freq,
                ylim = c(0.15,0.4),
                color = seq,
                data  = freqs,
                geom  = "line" )
    P <- P + facet_grid( seq ~ . )
    png(paste(f, "non_fitted", ".png", sep = ""))
    print(P)
    dev.off()
  } else {
    par( mar = c(5,1,1,1), mfrow = c(4,1) )
    with( subset( freqs, seq == "a" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "a", type = 'l' ))
    with( subset( freqs, seq == "c" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "c", type = 'l' ))
    with( subset( freqs, seq == "g" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "g", type = 'l' ))
    with( subset( freqs, seq == "t" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "t", type = 'l' ))
  } 
  
  
  sb <- seqbias.fit(ref_fn, reads_fn, L = 5, R = 15)
  bias <- seqbias.predict(sb, I)
  counts.adj <- mapply(FUN = `/`, counts, bias, SIMPLIFY = F)
  freqs.adj <- kmer.freq(seqs, counts.adj)
  
  if( require(ggplot2) ) {
    P <- qplot( x = pos,
                y = freq,
                ylim = c(0.15,0.4),
                color = seq,
                data  = freqs.adj,
                geom  = "line" )
    P <- P + facet_grid( seq ~ . )
    png(paste(f, "fitted", ".png", sep = ""))
    print(P)
    dev.off()
  } else {
    par( mar = c(5,1,1,1), mfrow = c(4,1) )
    with( subset( freqs.adj, seq == "a" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "a", type = 'l' ))
    with( subset( freqs.adj, seq == "c" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "c", type = 'l' ))
    with( subset( freqs.adj, seq == "g" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "g", type = 'l' ))
    with( subset( freqs.adj, seq == "t" ),
          plot( freq ~ pos, ylim = c(0.15,0.4), sub = "t", type = 'l' ))
  }
  seqbias.save(sb, paste(f, "_fitted_model.yaml"))
  
}
