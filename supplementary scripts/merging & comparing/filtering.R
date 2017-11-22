library("xlsx")
setwd("~/GitHub//counts/ALS Mice/new filtering/tg1-tg2/")
man <- read.xlsx("~/GitHub//counts/ALS Mice/new filtering/tg1-tg2/manifestation_genes_only_12.xlsx", sheetIndex = 1)  
diff <- read.xlsx("~/GitHub/counts/ALS Mice/experimental/results/tg_1-tg_2/Results edgeR.xlsx", sheetIndex = 2)
man <- as.data.frame(man$Ensembl.ID)
manifestation <- data.frame()
names(man) <- c("genes")

df <- data.frame()

for (f in man$genes){
  a <- grep(paste(f), diff$NA.)
  df <- diff[a,]
  manifestation <- rbind(df, manifestation)
}

et_annot <- as.data.frame(manifestation)
write.csv(et_annot, file = "~/GitHub/counts/ALS Mice/new filtering/tg1-tg2/manifestation_deg_12.csv")
