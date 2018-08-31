#!/usr/bin/Rscript

args <- commandArgs()
dir <- args[6]
#dir <- ("~/transcriptomes/reads/AIKAR.17.05.18/bias/")

if(nchar(dir) < 1){
  print("Typical usage: gene_coverage_plot /path/to/folder/with/files/. Specify folder")
}

print(dir)
setwd(dir)
g <- grep("txt", list.files(dir), value = T)

if (length(g) < 1){
  print("No files specified! Check your folder!")
  break
}

print(paste("Directory is" , dir, sep = " "))
data <- data.frame()

for (f in g){
  df <- read.delim(f)
  if (ncol(df) < 101){
    print(paste(f, "Is not output of geneBody_coverage.py! Skipping", sep = " "))
    next
  }
  data <- rbind(df, data)
  print(paste(f, "Binded!", sep = " "))
}

rownames(data) <- data$Percentile
data$Percentile <- NULL
data <- data.frame(t(data))
pdf(file = "3'-Bias Coverage Plot.pdf", width = 20, height = 20, family = "Helvetica")
matplot(data, type = "l", pch=1,col = 1:6)
legend("topright", inset=0.01, legend=colnames(data), col=c(1:6),pch=1,
       bg= ("white"), horiz=F)
dev.off()
print("DONE!")


