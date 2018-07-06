library(ggplot2)
library(reshape2)
library(rtracklayer)
library(dplyr)
library(gridExtra)
library(GenomicAlignments)
library(lattice)
library(ggbiplot)
library(mgcv)
library(data.table)
all.scaffolds <- read.table("~/genomes/12_genomes/Dvir/dvir-all-chromosome-r1.06.fasta.fai")
all.scaffolds <- as.character(all.scaffolds$V1)
setwd('~/chip-seq/Aravin/')
x <- "H3K9_1.nonunique.sorted.bam"



bam <- readGAlignments(x)
cov <- coverage(bam)
xcov <- as.numeric(cov$`scaffold_12966`)
xcov <- as.data.frame(xcov)
xcov$idx <- seq(1:nrow(xcov))





