source("https://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite("seqbias")
library("seqbias")
library("Rsamtools")

readBAM <- function(bamFile){
  
  bam <- scanBam(bamFile)
  

  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field

  return(bam_df)
}

reads_fn <- readBAM("~/BAM/1.bam")

ref_fn <- read.fasta(file = system.file("~/BAM/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa", package = "seqbias")
ref_f <- FaFile(ref_fn)

open.FaFile(ref_f)
reads_fn <- system.file("~/BAM/1.bam", package = "seqbias")
