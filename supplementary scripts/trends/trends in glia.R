library(xlsx)
d12 <- as.data.frame(read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg2/glia/deg_glia_12.csv"))
d23 <- as.data.frame(read.csv("~/GitHub/counts/ALS Mice/new filtering/tg2-tg3/glia/deg_glia_23.csv"))

a <- as.data.frame(intersect(d12$NA., d23$NA.))
df <- data.frame()
tr1 <- data.frame()
tr2 <- data.frame()
d12$NA.

for (f in a$`intersect(d12$NA., d23$NA.)`){
a <- grep(paste(f), d12$NA.)
df <- d12[a,]
tr1 <- rbind(df, tr1)
}

for (f in a$`intersect(d12$NA., d23$NA.)`){
a <- grep(paste(f), d23$NA.)
df <- d23[a,]
tr2 <- rbind(df, tr2)
}


tr <- data.frame(tr1$X.1, tr1$logFC, tr2$logFC)
