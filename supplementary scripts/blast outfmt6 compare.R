dir <- c("~/blast outfmt6 compare/pep/hp1d/")
setwd(dir)
g <- grep("blastp", list.files(dir), value = T)

df <- data.frame()
for (f in g){
  df.inter <- read.table(f)
  df.inter$V1 <- paste(f)
  df <- rbind(df.inter, df)
}

names(df) <- c("query", "subject", "pident", "length",
               "mismatch", "gapopen", "qstart", "qend",
               "sstart", "send", "evalue", "bitscore")

