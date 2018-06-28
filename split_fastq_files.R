dir <- "~/transcriptomes/reads/AIKAR.dmelanogaster/"
setwd(dir)
x <- grep("splitted", list.files(dir), value = TRUE)
for (f in x){
  z <- grep(strtrim(f, width = 10), list.files(dir), value=TRUE)
  s <- subset(z, !grepl("splitted", z))
  a <- paste("cat",s[1],s[2],s[3],s[4],">",f)
  print("CHECK YOUR DATA!")
  print(a)
  system(paste(a))
}
