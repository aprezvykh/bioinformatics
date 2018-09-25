dir <- "~/transcriptomes/reads/intron_retention/merge/"
g <- grep("\\b.bam\\b", list.files(dir), value = T)
g <- g[grep("bai", g,invert = T)]
paste(g,collapse = ",")
paste(seq(1:36),collapse = ",")
