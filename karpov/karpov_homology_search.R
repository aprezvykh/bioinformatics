source("https://bioconductor.org/biocLite.R")
biocLite("org.Sc.sgd.db")
library("org.Sc.sgd.db")

fa <- read.fasta("~/Karpov.sac/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",as.string = T)
gtf <- data.frame(import("~/Karpov.sac/Saccharomyces_cerevisiae.R64-1-1.93.gtf"))

chr.ids <- names(fa)
as.character(gregexpr("ggttagg",fa$I,ignore.case = T)[[1]])[1]



df <- read.delim("~/Karpov.sac/genes/genes.blast")
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

keytypes(org.Sc.sgd.db)
columns(org.Sc.sgd.db)

dfnu$qseqid.genename <- mapIds(org.Sc.sgd.db, 
                               keys=as.character(dfnu$qseqid), 
                               column="DESCRIPTION", 
                               keytype="ENSEMBL",
                               multiVals="first")
View(dfnu)
