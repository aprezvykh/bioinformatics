rep.dir <- "~/genomes/comp.transposones/rep/"
fa.dir <- "~/genomes/comp.transposones/"
rep <- grep(".fa", list.files(rep.dir), value = TRUE)
fa <- grep(".fa", list.files(fa.dir), value = TRUE)
df <- data.frame(rep,fa)
setwd(fa.dir)
for (i in seq(1:9)){
  print(i)
  system(paste(" ~/sofware/RepeatMasker/./RepeatMasker -lib rep/", 
               rep[i],
               " ",
               fa[i],
               " -pa 64 -e wublast -gff",
               sep = ""))
}

