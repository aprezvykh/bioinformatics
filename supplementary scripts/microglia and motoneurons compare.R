library(xlsx)
glia <- read.csv("~/GitHub/counts/ALS Mice/new filtering/microglia_expression_profile.csv")
moto <- read.csv("~/GitHub//counts/ALS Mice/new filtering/moto_expression_profile.csv")

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

exp <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg2/manifestation_deg_12.csv")
x <- read.xlsx(file = "~/GitHub/counts/ALS Mice/old diff/motoneurons marker.xlsx", sheetIndex = 1)

in1 <- intersect(moto_aff$compare.glia.X, exp$NA.)

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


sig_moto <- as.data.frame(intersect(in1, x$NA.))
sig_glia <- as.data.frame(intersect(gl_aff$compare.glia.X, exp$NA.))

a <- outersect(sig_moto$`intersect(in1, x$NA.)`, exp$NA.)
b <- outersect(sig_glia$`intersect(gl_aff$compare.glia.X, exp$NA.)`, exp$NA.)

sig_out <- as.data.frame(intersect(a,b))

deg_moto <- data.frame()
deg_glia <- data.frame()
deg_out <- data.frame()

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

for (f in sig_out$`intersect(a, b)`){
  a <- grepl(paste(f), exp$NA.)
  r <- exp[a,]
  deg_out <- rbind(r, deg_out)
}

setwd("~/GitHub/counts/ALS Mice/new filtering/tg1-tg2/")

deg_moto <- deg_moto[complete.cases(deg_moto), ]
deg_glia<- deg_glia[complete.cases(deg_glia), ]
deg_out<- deg_out[complete.cases(deg_out), ]
write.csv(deg_moto, file = "deg_moto_12.csv")
write.csv(deg_glia, file = "deg_glia_12.csv")
write.csv(deg_out, file = "deg_out_12.csv")
