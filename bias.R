source("https://bioconductor.org/biocLite.R")
biocLite("seqbias")
library("seqbias")
library(Rsamtools)
library(devtools)


ref_fn <- "~/BAM/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
ref_f <- FaFile(ref_fn)
open.FaFile(ref_f)
ref_seqs <- scanFaIndex(ref_f)
I <- random.intervals(ref_seqs, n = 5, m = 100000)
seqs <- scanFa(ref_f, I)
