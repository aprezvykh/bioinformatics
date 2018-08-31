#!/usr/bin/env Rscript
library("txtplot")
a <- read.table("data.tsv")
#a$V4 <- sub("T", "", a$V4)
a$V4 <- as.numeric(a$V4)
txtplot(a$V4)

