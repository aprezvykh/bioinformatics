library(xlsx)
glia <- read.csv("~/bioinformatics/counts/ALS Mice/compare_new/microglia.csv")
moto <- read.csv("~/bioinformatics//counts/ALS Mice/compare_new/motoneurons.csv")

compare <- data.frame(glia$X, glia$sum, moto$sum)
compare <- subset(compare, ((compare$glia.sum + compare$moto.sum)/2) > 8)

compare$glia.affinity <- compare$glia.sum/compare$moto.sum
compare$moto.affinity <- compare$moto.sum/compare$glia.sum

gl_aff <- data.frame(compare$glia.X, compare$glia.affinity)
moto_aff <- data.frame(compare$glia.X, compare$moto.affinity)

gl_notaff <- subset(gl_aff, gl_aff$compare.glia.affinity < 5)
moto_notaff <- subset(moto_aff, moto_aff$compare.moto.affinity < 5)

gl_aff <- subset(gl_aff, gl_aff$compare.glia.affinity > 5)
moto_aff <- subset(moto_aff, moto_aff$compare.moto.affinity > 5)

common_genes <- data.frame()

exp <- read.xlsx(file = "~/bioinformatics/counts/ALS Mice/experimental/results/tg_1-tg_3/Results edgeR.xlsx", sheetIndex = 2)
x <- read.xlsx(file = "~/bioinformatics/counts/ALS Mice/motoneurons marker.xlsx", sheetIndex = 1)

in1 <- intersect(moto_aff$compare.glia.X, exp$NA.)

sig_moto <- as.data.frame(intersect(in1, x$NA.))
sig_glia <- as.data.frame(intersect(gl_aff$compare.glia.X, exp$NA.))

deg_moto <- data.frame()
deg_glia <- data.frame()
for (f in sig_moto$`intersect(in1, x$NA.)`){
  a <- grepl(paste(f), exp$NA.)
  r <- exp[a,]
  deg_moto <- rbind(r, deg_moto)
}

for (f in sig_glia$`intersect(gl_aff$compare.glia.X, exp$NA.)`){
  a <- grepl(paste(f), exp$NA.)
  r <- exp[a,]
  deg_glia <- rbind(r, deg_glia)
}

deg_moto <- deg_moto[complete.cases(deg_moto), ]
deg_glia<- deg_glia[complete.cases(deg_glia), ]

write.csv(deg_moto, file = "deg_moto.csv")
write.csv(deg_glia, file = "deg_glia.csv")


