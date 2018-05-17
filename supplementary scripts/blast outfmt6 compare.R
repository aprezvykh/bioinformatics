library(ggplot2)
dir <- c("~/my_outfmt/")
setwd(dir)
g <- grep("outfmt6", list.files(dir), value = T)

for (f in g){
    print(f)
    df <- read.table(f)
    names(df) <- c("query", "subject", "pident", "length",
                   "mismatch", "gapopen", "qstart", "qend",
                   "sstart", "send", "evalue", "bitscore")

    g1 <- ggplot(data=df) + geom_point(aes(x = pident, y = length, col = bitscore, size = mismatch)) + 
                      geom_smooth(aes(x = pident, y = length))
    ggsave(filename = paste(f, ".png", sep = ""), device = "png")
    
    df.ordered <- df[order(df$sstart),]
    df.ordered$num <- seq(1:nrow(df.ordered))
    g2 <- ggplot(data=df.ordered) + geom_line(aes(x=num, y=sstart))
    g3 <- ggplot(data=df.ordered) + geom_smooth(aes(x=num, y=bitscore))
}





