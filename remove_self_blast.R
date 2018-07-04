dir <- "~/genomes/circos/outfmt/"
setwd(dir)
l <- grep("outfmt6", list.files(dir), value = TRUE)
for (i in 1:length(l)){
    f <- strtrim(strsplit(l[i], "_")[[1]][1],4)
    s <- strtrim(strsplit(l[i], "_")[[1]][2],4)
    if(f == s){
      system(paste("rm ", l[i], sep = ""))
    }
    
}


