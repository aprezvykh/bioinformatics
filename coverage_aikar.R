setwd("~/coverage/")
x <- read.delim(file = "Transcripts Read Coverage table.tsv", header = T)

larv <- grep("adult", colnames(x))
y <- x[,larv]

matplot(y, type = "l")

plot(x$X1.3rd.insta_splitted, type = "l")
