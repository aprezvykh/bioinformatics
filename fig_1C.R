setwd("~/genomes/IR.w.flanks.20kb/")
df <- read.table("rmask.out")

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


myb.3 <- 20001
myb.5 <- 22879
ran.3 <- 143218
ran.5 <- 192819

df[grepl("DNA", df$class),]$class <- "DNA"
df[grepl("LINE", df$class),]$class <- "LINE"
df[grepl("LTR", df$class),]$class <- "LTR"
df[grepl("Simple", df$class),]$class <- "Satellite"
df[grepl("Low", df$class),]$class <- "Satellite"

df <- df[!grepl("Unknown", df$class),]

df$sum.inv <- ifelse(test = df$perc > 0, "Straight", "Inverted")

ggplot(data = df) + 
  geom_segment(aes(x=qbegin, 
                   xend=qend,
                   y=perc, 
                   yend=perc,
                   color = class), 
               size = 0.1 ,
               arrow = arrow(length = unit(0.3, "cm"), 
                               ends = ifelse(df$sum.inv == "Straight", "first", "last"), type = "closed",angle = 30)) + 
  theme_bw() + 
  geom_segment(aes(x = myb.3, xend = myb.5, y = 0, yend = 0), size = 1.5) + 
  geom_segment(aes(x = ran.3, xend = ran.5, y = 0, yend = 0), size = 1.5) + 
  geom_hline(yintercept=0) + 
  geom_text(aes(x = (myb.3+myb.5)/2, y = -10), label = "Myb") + 
  geom_text(aes(x = (ran.3+ran.5)/2, y = -10), label = "Ranbp16") + 
  scale_x_continuous(breaks = c(0,50000,100000,150000,200000), labels = c("0kb", "50kb","100kb","150kb","200kb"),name = "Genomic region") + 
  scale_y_continuous("% of query coverage")
  
