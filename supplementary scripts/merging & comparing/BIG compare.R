cd4 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/cd4.csv")
cd4$avg <- cd4$mouse_SRR391076.counts
cd22 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/cd22.csv")
moto <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/moto.csv")
glia <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/microglia.csv")
tg1 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/tg1.csv")
tg2 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/tg_2.csv")
tg3 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/tg_3.csv")
ntg1 <- read.csv("~/GitHub/counts/ALS Mice/immune cells compare/control_1.csv")

dflist <- list("cd4", "cd22", "moto","glia", "tg1", "tg2", "tg3", "ntg1")

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
cd4 <- cd4[order(cd4$avg, decreasing = TRUE),]
cd22 <- cd22[order(cd22$avg, decreasing = TRUE),]
moto <- moto[order(moto$avg, decreasing = TRUE),]
glia <- glia[order(glia$avg, decreasing = TRUE),]
tg1 <- tg1[order(tg1$avg, decreasing = TRUE),]
tg2 <- tg2[order(tg2$avg, decreasing = TRUE),]
tg3 <- tg3[order(tg3$avg, decreasing = TRUE),]
ntg1 <- ntg1[order(ntg1$avg, decreasing = TRUE),]



cd4 <- cd4[complete.cases(cd4),]
cd22 <- cd22[complete.cases(cd22),]
moto <- moto[complete.cases(moto),]
glia <- glia[complete.cases(glia),]
tg1 <- tg1[complete.cases(tg1),]
tg2 <- tg2[complete.cases(tg2),]
tg3 <- tg3[complete.cases(tg3),]
ntg1 <- ntg1[complete.cases(ntg1),]

i <- intersect(intersect(cd4$X, cd22$X), intersect(cd22$X, glia$X))
all.common <- intersect(i, moto$X)

cd4_common <- data.frame()
cd22_common <- data.frame()
moto_common <- data.frame()
glia_common <- data.frame()





for (f in all.common){
a <- grep(paste(f), cd4$X)
df <- cd4[a,]
cd4_common <- rbind(df, cd4_common)
}

for (f in all.common){
  a <- grep(paste(f), cd22$X)
  df <- cd22[a,]
  cd22_common <- rbind(df, cd22_common)
}

for (f in all.common){
  a <- grep(paste(f), moto$X)
  df <- moto[a,]
  moto_common <- rbind(df, moto_common)
}

for (f in all.common){
  a <- grep(paste(f), glia$X)
  df <- glia[a,]
  glia_common <- rbind(df, glia_common)
}




cd4 <- cd4[seq(1:1000),]
cd22 <- cd22[seq(1:1000),]
moto <- moto[seq(1:1000),]
glia <- glia[seq(1:1000),]
tg1 <- tg1[seq(1:1000),]
tg2 <- tg2[seq(1:1000),]
tg3 <- tg3[seq(1:1000),]
ntg1 <- ntg1[seq(1:1000),]


df.all.avg <- data.frame(all.common, cd4_common$avg, cd22_common$avg, moto_common$avg, glia_common$avg)
names(df.all.avg) <- c("gene", "cd4", "cd22", "moto", "glia")
write.csv(df.all.avg, file = "ALLCOMMONAVERAGE.csv")


aff <- data.frame()

for (f in seq(1:nrow(df.all.avg))){
  a <- data.frame(df.all.avg[f,])
  
  if ((df.all.avg$cd4/df.all.avg$cd22) >= 5 && (df.all.avg$cd4/df.all.avg$moto) >= 5 && (df.all.avg$cd4/df.all.avg$glia) >= 5){
    aff <- rbind 
    
  } else {
    print("false")
}
  
}



df.all.avg$`cd4/22` <- df.all.avg$cd4/df.all.avg$cd22
df.all.avg$`cd4/moto` <- df.all.avg$cd4/df.all.avg$moto
df.all.avg$`cd4/glia` <- df.all.avg$cd4/df.all.avg$glia

df.all.avg$`cd22/4` <- df.all.avg$cd22/df.all.avg$cd4
df.all.avg$`cd22/moto` <- df.all.avg$cd22/df.all.avg$moto
df.all.avg$`cd22/glia` <- df.all.avg$cd22/df.all.avg$glia

df.all.avg$`moto/cd4` <- df.all.avg$moto/df.all.avg$cd4
df.all.avg$`moto/cd22` <- df.all.avg$moto/df.all.avg$cd22
df.all.avg$`moto/glia` <- df.all.avg$moto/df.all.avg$glia

df.all.avg$`glia/cd4` <- df.all.avg$glia/df.all.avg$cd4
df.all.avg$`glia/cd22` <- df.all.avg$glia/df.all.avg$cd22
df.all.avg$`glia/moto` <- df.all.avg$glia/df.all.avg$moto


aff_cd4 <- data.frame()
aff_cd22 <- data.frame()
aff_glia <- data.frame()
aff_moto <- data.frame()



g <- data.frame(df.all.avg$`cd4/22`, df.all.avg$`cd4/moto`, df.all.avg$`cd4/glia`)
rownames(g) <- df.all.avg$gene

for (s in seq(1:nrow(g))){
a <- data.frame(g[s,])

  if (a$df.all.avg..cd4.22. > 5 && a$df.all.avg..cd4.moto. < 1 && a$df.all.avg..cd4.glia. < 1){
  aff_cd4 <- rbind(a, aff_cd4)
} else {
  print("no!")
}

}

g <- data.frame(df.all.avg$`cd22/4`, df.all.avg$`cd22/moto`, df.all.avg$`cd22/glia`)
rownames(g) <- df.all.avg$gene


for (s in seq(1:nrow(g))){
  a <- data.frame(g[s,])
  
  if (a$df.all.avg..cd22.4. > 5 && a$df.all.avg..cd22.moto. < 1 && a$df.all.avg..cd22.glia. < 1){
    aff_cd22 <- rbind(a, aff_cd22)
  } else {
    print("no!")
  }
  
}

g <- data.frame(df.all.avg$`moto/cd4`, df.all.avg$`moto/cd22`, df.all.avg$`moto/glia`)
rownames(g) <- df.all.avg$gene

a <- g[1,]
for (s in seq(1:nrow(g))){
  a <- data.frame(g[s,])
  
  if (a$df.all.avg..moto.cd4. > 5 && a$df.all.avg..moto.cd22. < 1 && a$df.all.avg..moto.glia. < 1){
    aff_moto <- rbind(a, aff_moto)
  } else {
    print("no!")
  }
  
}

g <- data.frame(df.all.avg$`glia/cd4`, df.all.avg$`glia/cd22`, df.all.avg$`glia/moto`)
rownames(g) <- df.all.avg$gene

a <- g[1,]
for (s in seq(1:nrow(g))){
  a <- data.frame(g[s,])
  
  if (a$df.all.avg..glia.cd4. > 5 && a$df.all.avg..glia.cd22. < 1 && a$df.all.avg..glia.moto. < 1){
    aff_glia <- rbind(a, aff_glia)
  } else {
    print("no!")
  }
  
}

library(org.Mm.eg.db)

aff_cd4$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(aff_cd4), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

aff_cd4$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(aff_cd4), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

aff_cd22$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(aff_cd22), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

aff_cd22$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(aff_cd22), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")


aff_glia$Symbol <- mapIds(org.Mm.eg.db, 
                         keys=row.names(aff_glia), 
                         column="SYMBOL", 
                         keytype="ENSEMBL",
                         multiVals="first")

aff_glia$Name <- mapIds(org.Mm.eg.db, 
                       keys=row.names(aff_glia), 
                       column="GENENAME", 
                       keytype="ENSEMBL",
                       multiVals="first")

aff_moto$Symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(aff_moto), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

aff_moto$Name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(aff_moto), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")



in1 <- intersect(rownames(aff_cd22), tg1$X)
in3 <- intersect(rownames(aff_cd22), tg3$X)

finaldiff <- data.frame(outersect(in1, in3))
rownames(finaldiff) <- finaldiff$outersect.in1..in3.

finaldiff$Symbol <- mapIds(org.Mm.eg.db, 
                          keys=row.names(finaldiff), 
                          column="SYMBOL", 
                          keytype="ENSEMBL",
                          multiVals="first")

finaldiff$Name <- mapIds(org.Mm.eg.db, 
                        keys=row.names(finaldiff), 
                        column="GENENAME", 
                        keytype="ENSEMBL",
                        multiVals="first")


