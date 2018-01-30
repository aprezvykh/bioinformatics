dir <- c('~/coverage/')
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


plot(df$`8.bam.geneBodyCoverage.txt`, type = "l")


sas <- data.frame()


for (f in 1:ncol(df)){
  s <- data.frame(df[,f], seq(1:100))
  colnames(s) <- c("intensity", "index")
  m <- which.max(s$intensity)
  sas <- rbind(m, sas)
  print(m)
}


mean(sas$X91L)







pallete = c("#F46D43", "#66C2A5", "#cd8845", "#3288BD", "#a8bf32", "#5E4FA2", "#D53E4F", "#d6d639", "#8ed384", "#9E0142", "#ebba2f")
density.cols = colorRampPalette(pallete)(dim(df)[2])
#pdf("Coverage.pdf", height = 10, width = 10, family = "Helvetica")
matplot(df,type = "l", pch=1, col = 1:n, xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8) #plot
legend("topleft", legend = leg, col=1:n, pch=1)
#dev.off()



