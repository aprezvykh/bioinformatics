library(rtracklayer)
library(Biostrings)


gtf <- rtracklayer::import("~/dvir_gtf_parse/dvir-all-r1.06.gtf")
df <- as.data.frame(gtf)

f <- readDNAStringSet(filepath = '~/dvir_gtf_parse/dvir-all-miRNA-r1.06.fasta')

mir.transcripts <- data.frame(0)

for (z in 1:319){
s <- as.data.frame(strsplit(names(f[z]), ";"))
ssub <- sub("ID=", "", s[3,])
trans <- sub(" ", "", ssub)
print(trans)
mir.transcripts <- rbind(trans, mir.transcripts)
}

names(mir.transcripts) <- c("trans")



df.parsed <- df[which(df$transcript_id %in% mir.transcripts$trans),]
rtracklayer::export(df.parsed, "dvir_mir_parsed.gtf", format = "GTF")
  