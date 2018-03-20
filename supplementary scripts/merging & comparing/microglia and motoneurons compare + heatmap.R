### BASIC FILTERING
library(xlsx)
library(ggplot2)
library(gplots)
glia <- read.csv("~/counts/ALS Mice/new filtering/microglia_expression_profile.csv")
moto <- read.csv("~/counts/ALS Mice/new filtering/moto_expression_profile.csv")

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

exp <- as.data.frame(read.xlsx("~/GitHub/counts/ALS Mice/new filtering-0.5/exp_tg1_tg2.xlsx", sheetIndex = 2))
x <- read.xlsx(file = "~/GitHub/counts/ALS Mice/old diff/motoneurons marker.xlsx", sheetIndex = 1)

exp <- exp[complete.cases(exp), ]
flt <- as.data.frame(read.xlsx("~/GitHub/counts/ALS Mice/tg1-tg2/manifestation_genes_only_12.xlsx", sheetIndex = 1))
r <- data.frame()
df <- data.frame()
names(flt) <- c("id")

for (f in flt$id){
  a <- grepl(paste(f), exp$NA.)
  df <- exp[a,]
  r <- rbind(df, r)
  
}

exp <- r
in1 <- intersect(moto_aff$compare.glia.X, r$NA.)
write.xlsx(r, "Filtered_only_13.xlsx", sheetName = "12")

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

setwd("~/GitHub/counts/ALS Mice/tg1-tg2/")

deg_moto <- deg_moto[complete.cases(deg_moto), ]
deg_glia<- deg_glia[complete.cases(deg_glia), ]
deg_out<- deg_out[complete.cases(deg_out), ]
write.csv(deg_moto, file = "deg_moto_13.csv")
write.csv(deg_glia, file = "deg_glia_13.csv")
write.csv(deg_out, file = "deg_out_13.csv")

#################


man <- read.csv("~/GitHub/counts/ALS Mice/new filtering/tg1-tg3/manifestation_deg_13.csv")

df <- NULL

df_glia <- data.frame()
for (f in man$NA.){
  a <- grep(paste(f), glia$X)
  r <- glia[a,]
  df_glia <- rbind(r, df_glia)
}
df_glia$sum <- NULL
df_glia$Symbol <- NULL
df_glia$Name <- NULL

df_moto <- data.frame()
for (f in man$NA.){
  a <- grep(paste(f), moto$X)
  r <- moto[a,]
  df_moto <- rbind(r, df_moto)
}
df_moto$sum <- NULL
df_moto$Symbol <- NULL
df_moto$Name <- NULL

hm <- data.frame()

hm <- df_glia
col.pan <- colorpanel(100, "blue", "white", "red")
hm <- cbind(df_moto, hm)
rownames(hm) <- hm$X

hm$X <- NULL
hm$X <- NULL


hm <- as.matrix(hm)
hm <- t(scale(t(hm)))
setwd("~/")
pdf(file = "Compare 2.pdf", width = 10, height = 10)
heatmap.2(hm, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6), main = "Genes differential expression")
dev.off()
#hm$moto <- (hm$moto_1.counts + hm$moto_2.counts)/2
#hm$glia <- (hm$mouse_1.counts +
#              hm$mouse_2.counts +
#                hm$mouse_3.counts +
#                  hm$mouse_4.counts +
#                    hm$mouse_5.counts + 
#                      hm$mouse_6.counts)/6

#fit_hm <- data.frame(hm$moto, hm$glia)
#rownames(fit_hm) <- rownames(hm)
#fit_hm <- as.matrix(fit_hm)
#fit_hm <- t(scale(t(fit_hm)))

#hm <- as.data.frame(hm)
#hm <- data.frame(hm$moto, hm$glia)
#hm <- as.matrix(hm)


intersect(gl_aff$compare.glia.X, moto_aff$compare.glia.X)

directory <- '~/counts/ALS Mice/experimental/'
setwd(directory)
sampleFiles <- grep('mouse',list.files(directory),value=TRUE)
sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                     'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                     'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2', 
                     'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')


sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
y <- readDGE(files = sampleTable$sampleName, group = sampleTable$condition, labels = sampleTable$fileName)
df <- data.frame(y$counts)


colnames(df)
g <- df[(rownames(df) %in% gl_aff$compare.glia.X),]
m <- df[(rownames(df) %in% moto_aff$compare.glia.X),]

cs.g <- as.data.frame(colSums(g))
cs.m <- as.data.frame(colSums(m))
cs.split <- data.frame(cs.g$`colSums(g)`, cs.m$`colSums(m)`)
rownames(cs.split) <- rownames(cs.g)






c1.mean.g <- as.data.frame(t(lapply(df[,grep("Control.1", colnames(df))], mean)))
