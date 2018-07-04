setwd("~/genomes/circos/outfmt/")
df <- read.delim("split.out", header = F)
links <- data.frame(df$V1, df$V7, df$V8,
                    df$V2, df$V9, df$V10)
names(links) <- c("query", "qstart", "qend",
                  "subject", "sstart", "send")

#links$query <- as.character(links$query)
#links$subject <- as.character(links$query)

big.df <- data.frame()
for (i in seq(from=2, to=nrow(links), by = 2)){
  f.row <- data.frame(paste("link", i+10, sep = ""), links$query[i-1], links$qstart[i-1], links$qstart[i-1], paste("id=", i+20, sep = ""))
  s.row <- data.frame(paste("link", i+10, sep = ""), links$subject[i], links$sstart[i], links$qend[i],paste("id=", i+20, sep = ""))
  names(f.row) <- c("link","chr", "start", "end","id")
  names(s.row) <- c("link","chr", "start", "end","id")
  sub.df <- f.row
  sub.df <- rbind(s.row, f.row)
  big.df <- rbind(sub.df, big.df)
  }

write.table(big.df,file = "my_links.txt")
getwd()
write.table(links,file = "my_links_3.txt", row.names = F, col.names = F)
