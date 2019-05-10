df <- read.delim("~/Karpov.sac/homolog.genes/blast.outfmt6")
names(df) <- c("qseqid",
               "sseqid",
               "pident",
               "length",
               "mismatch",
               "gapopen", 
               "qstart", 
               "qend", 
               "sstart", 
               "send", 
               "evalue", 
               "bitscore")

dfnu <- df[df$qseqid != df$sseqid,]
dfnu <- dfnu[dfnu$pident > 98,]
dfnu <- dfnu[dfnu$length > 500,]

dfnu <- dfnu[order(dfnu$bitscore, decreasing = T),]

as.character(unique(dfnu$qseqid))

dfnu[dfnu$qstart == 1 & dfnu$sstart == 1,]
