dir <- "~/transcriptomes/reads/intron_retention/all_counts/"
setwd(dir)
g <- grep("counts", list.files(dir), value = T)

for(f in g){
    print(f)
    df <- read.delim(f,header = F)
    ndf <- data.frame(df$V1, df$V7)
    names(ndf) <- NULL
    write.table(ndf, paste(f, ".parsed.counts", sep = ""),row.names=FALSE,sep="\t", quote = FALSE,col.names = F)
}

