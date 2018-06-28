dir <- "~/transcriptomes/reads/virilis_101/"
setwd(dir)

x <- grep("1", list.files(dir), value = TRUE)

for (f in x){
setwd(paste(dir, f, "/", "Files/", sep = ""))
system(paste("sudo ", "touch ", f, ".fastq", sep = ""))
}
