dir <- c('~/coverage.mouse/')
setwd(dir)
files <- grep("txt", list.files(dir), value = TRUE)
n = length(files)
df <- data.frame()
a <- data.frame()
for (f in files){
  a <- read.table(file = paste(f), header = TRUE)
  a$Percentile <- NULL
  a <- as.data.frame(a)
  rownames(a) <- paste(f)
  df <- rbind(a, df)
}

leg <- rownames(df)
df <- as.data.frame(t(df))

pdf("Coverage.pdf", height = 10, width = 10, family = "Helvetica")
matplot(df,type = "l", pch=1, col = 1:n, xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8) #plot
legend("topleft", legend = leg, col=1:n, pch=1)
dev.off()
