setwd("~/genomes/IR.w.flanks.20kb/")
df <- read.table("other_pts/1.fasta.out")

names(df) <- c("score", "div", "del", "ins", 
               "qseq", "qbegin", "qend", "qleft",
               "strand", "repeat_id", "class",
               "sbegin","send", "sleft")


df$sbegin <- gsub("\\(|\\)", "", df$sbegin)
df$sleft <- gsub("\\(|\\)", "", df$sleft)


df$sbegin <- as.numeric(df$sbegin)
df$send <- as.numeric(df$send)
df$sleft <- as.numeric(df$sleft)
df$sum <- df$sbegin + df$send + df$sleft
df$perc <- 100*((df$send - df$sbegin)/df$sum)

df <- data.frame(df)
df$class <- as.character(df$class)

unique(df$class)



df[grepl("DNA", df$class),]$class <- "DNA"
df[grepl("LINE", df$class),]$class <- "LINE"
df[grepl("LTR", df$class),]$class <- "LTR"
df[grepl("Simple", df$class),]$class <- "Satellite"
df[grepl("Low", df$class),]$class <- "Satellite"

df$sum.inv <- ifelse(test = df$perc > 0, "Straight", "Inverted")

ggplot(data = df) + 
  geom_segment(aes(x=qbegin, 
                   xend=qend,
                   y=perc, 
                   yend=perc,
                   color = class), 
               size = 0.1 ,
               arrow = arrow(length = unit(0.3, "cm"), 
                             ends = ifelse(df$sum.inv == "Straight", "first", "last"), type = "closed")) + 
  theme_bw()

