#!/usr/bin/Rscript
args <- commandArgs()

library("Biostrings", quietly = T)
print("Inserting spacer indexes!")

dir <- args[6]

file <- "potential_ngg.fasta.parsed"
df <- read.delim(file, stringsAsFactors = F, header = F)
numbers <- seq(1:length(df[grep(">", df$V1),]))
pams <- paste(">PAM_", numbers, sep = "")  
sequences <- df[grep(">", df$V1, invert = T),]
sequences <- unique(unlist(strsplit(sequences, " ")))
system("touch ngg.headers.fasta")
cat(paste(pams, '\n', sequences, sep = ''),sep = '\n',file = "ngg.headers.fasta")
