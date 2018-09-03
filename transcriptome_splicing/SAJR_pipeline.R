library(pheatmap)
library(devtools)
library(org.Mm.eg.db)
library(ggbiplot)
library(RColorBrewer)
library(gplots)
library(GenomicRanges)
library(SAJR)
library(genefilter)
library(reshape2)
library(gridExtra)
library(plotrix)
library(GO.db)

col.pan <- colorpanel(100, "blue", "white", "red")
setwd("~/Documents/intron_retention/")
data <- readRDS("data.rds")
data.f <- readRDS("data_filtered.rds")
all.alts <- readRDS("alts.rds")
gdata <- readRDS("gdata.rds")


#ref <- read.csv("mazin_ref.csv")

plot_ir <- function(gene,gr){
  sub <- data.f$ir[which(rownames(data.f$ir) == gene),]
  ens <- gsub("\\..*","",gene)
  sub[is.na(sub)] <- 0
  sub <- data.frame(sub)
  sub$group <- rownames(sub)
  sub$group <- gr
  names(sub) <- c("ir", "group")
  g1 <- ggplot(data=sub) + geom_boxplot(aes(x = group, y = ir, fill = group)) + 
    theme_bw() + 
    ggtitle(paste("Inclusion rate - ",gene, sep = ""))
  gd <- data.frame(t(gdata[grep(ens,rownames(gdata$cnts)),]$cnts))
  names(gd) <- c("gene")
  gd$group <- gr
  names(gd) <- c("expr", "gr")
  g2 <- ggplot(data=gd) + geom_boxplot(aes(x = gr, y = expr, fill = gr)) + 
    theme_bw() + 
    ggtitle(paste("Expression, corellation coefficient - ", round(cor(sub$ir, gd$expr),2)))
  grid.arrange(g1,g2)
}
make.stattest.and.annotate <- function(sajr.obj,matrix,stattest.id,control.id,case.id,psi.cutoff,pval.cutoff){
  st <- meta[unlist(lapply(names(meta),grepl,stattest.id),recursive = T)]
  g <- grep(names(st), names(sajr.obj$seg))
  pval <- sajr.obj$seg[,g]
  sub <- data.frame(sajr.obj$ir,pval)
  st <- as.character(unlist(st))
  names(sub) <- c(st,names(pval))
  sub[is.na(sub)] <- 0
  unique(st)
  sub$dpsi <- rowMeans(sub[,grep(case.id,names(sub))])-rowMeans(sub[,grep(control.id,names(sub))])
  if(nrow(sub)<1){
    warning("No significant genes at all! Check your data")
  }
  sub <- sub[which(abs(sub$dpsi) > psi.cutoff),]
  if(nrow(sub)<1){
    warning("No significant genes after dPSI cutoff! Check dPSI cutoff value!")
  }
  sub <- sub[which(sub[,grep("adj", names(sub))] < pval.cutoff),]
  if(nrow(sub)<1){
    warning("No significant genes after p-value cutoff! Check p-value!")
  }
  
  sub$ens_gene <- gsub("\\..*","",rownames(sub))
  
  sub$symbol <- mapIds(org.Mm.eg.db, 
                       keys=as.character(sub$ens_gene), 
                       column="SYMBOL", 
                       keytype="ENSEMBL",
                       multiVals="first")
  
  sub$name <- mapIds(org.Mm.eg.db, 
                     keys=as.character(sub$ens_gene), 
                     column="GENENAME", 
                     keytype="ENSEMBL",
                     multiVals="first")
  
  
  sub$entrez <- mapIds(org.Mm.eg.db, 
                       keys=as.character(sub$ens_gene), 
                       column="ENTREZID", 
                      keytype="ENSEMBL",
                       multiVals="first")
  
  sub$GOID <-     mapIds(org.Mm.eg.db, 
                         keys=as.character(sub$ens_gene), 
                         column="GO", 
                         keytype="ENSEMBL",
                         multiVals="first")
  
  sub$term <- mapIds(GO.db, 
                     keys=sub$GOID, 
                     column="TERM", 
                     keytype="GOID",
                     multiVals="first")
  
  sub$term <- as.character(sub$term)
  
#  alt <- data.f$seg[rownames(data.f$seg) %in% rownames(sub),]
#  sub <- cbind(alt,sub)
  View(sub)
  write.csv(sub, paste(stattest.id, "_", case.id, "-vs-", control.id, ".csv", sep = ""))
}
plot_boxplots <- function(data, group_vec){
  df <- data.frame(as.character(colMeans(data$ir)),stringsAsFactors = F)
  names(df) <- c("mean")
  df$mean <- as.numeric(as.character(df$mean))
  df$sds <- as.character(apply(data$ir, 2, sd))
  df$sds <- as.numeric(as.character(df$sds))
  df$gr <- group_vec
  names(df) <- c("mean", "SD", "group")
  for_m <- data.frame(df$mean,df$group)
  names(for_m) <- c("mean","group")
  g1 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) +
    ggtitle("Inclusion ratio by groups") + 
    theme_bw()
  
  
  df <- data.frame(as.character(colMeans(data$i)),stringsAsFactors = F)
  names(df) <- c("mean")
  df$mean <- as.numeric(as.character(df$mean))
  df$sds <- as.character(apply(data$i, 2, sd))
  df$sds <- as.numeric(as.character(df$sds))
  df$gr <- group_vec
  names(df) <- c("mean", "SD", "group")
  for_m <- data.frame(df$mean,df$group)
  names(for_m) <- c("mean","group")
  g2 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) + 
    ggtitle("Inclusion reads, mean by group") + theme_bw()
  
  
  df <- data.frame(as.character(colMeans(data$e)),stringsAsFactors = F)
  names(df) <- c("mean")
  df$mean <- as.numeric(as.character(df$mean))
  df$sds <- as.character(apply(data$e, 2, sd))
  df$sds <- as.numeric(as.character(df$sds))
  df$gr <- group_vec
  names(df) <- c("mean", "SD", "group")
  for_m <- data.frame(df$mean,df$group)
  names(for_m) <- c("mean","group")
  g3 <- ggplot(data=for_m) + geom_boxplot(aes(x = group, y = mean, fill = group), show.legend = F) + 
    ggtitle("Exclusion reads, mean by group") + theme_bw()
  grid.arrange(g1,g2,g3)
}
all_vars_pca <- function(data,group_vec){
  big.df <- data.frame()
  df <- data$ir
  big.df <- rbind(df, big.df)
  df <- data$i
  big.df <- rbind(df, big.df)
  df <- data$e
  big.df <- rbind(df, big.df)
  p <- prcomp(t(big.df), center = TRUE)
  ggbiplot(p, var.axes = F, groups = group_vec, ellipse = T) + theme_bw() + ggtitle("All vars PCA")
}
get_corr_per_group <- function(sajr_data, sajr_gdata,me){
  u <- unique(me)
  v <- vector()
  for(f in u){
    c <- cor(colMeans(sajr_data$ir[,grep(f,me)],na.rm = T),colMeans(sajr_gdata$cnts[,grep(f,me)],na.rm = T))
    v <- append(c,v)
  }
  names(v) <- u
  v <- as.data.frame(v)
  v$n <- rownames(v)
  g <-ggplot(data=v) + geom_bar(aes(x = n, y = v, fill = n), stat = "identity") + 
    theme_bw() + 
    ggtitle("Coefficient of corellation between IR and counts")
  return(g)
  
}

###Annotate events
data.f$seg$sites <- sub("ad","cassete_exon",data.f$seg$sites)
data.f$seg$sites <- sub("aa","alternative_acceptor",data.f$seg$sites)
data.f$seg$sites <- sub("dd","alternative_donor",data.f$seg$sites)
data.f$seg$sites <- sub("da","retainted_intron",data.f$seg$sites)
View(data.f$seg)
###COMMON CORALLATION AND STUFF
samples <- seq(1:24)

transgenity <- c("tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt")

age <- c("60d", "60d", "60d", "60d", "60d", 
         "90d", "90d", "90d", "90d", 
         "120d", "120d", "120d", "120d", "120d", 
         "60d", "60d", "60d", "60d", "60d", 
         "120d", "120d", "120d", "120d", "120d")

groups <- c('tg1', 'tg1', 'tg1', 'tg1', 'tg1', 
            'tg2', 'tg2', 'tg2', 'tg2', 
            'tg3', 'tg3', 'tg3', 'tg3', 'tg3', 
            'wt1', 'wt1', 'wt1', 'wt1', 'wt1', 
            'wt3', 'wt3', 'wt3', 'wt3', 'wt3')


meta <- list(samples=samples,transgenity=transgenity, age = age, groups = groups)

n <- meta$groups
data.f$ir[is.na(data.f$ir)] <- 0
p <- prcomp(t(data.f[rowSums(data.f$ir)>0,]$ir), center = TRUE)
ggbiplot(p, var.axes = F, groups = meta$groups, ellipse = T) + theme_bw() + ggtitle("Inclusion rate, PCA")

p <- prcomp(t(data.f[rowSums(data.f$i)>0,]$i), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Inclusion reads, PCA")

p <- prcomp(t(data.f[rowSums(data.f$e)>0,]$e), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Exclusion reads, PCA")


pheatmap(cor(data.f[rowSums(data.f$ir)>0,]$ir, method = "pearson"),color = col.pan,main = "Inclusion ratio, Spearman corr.")
pheatmap(cor(data.f[rowSums(data.f$i)>0,]$i, method = "pearson"),color = col.pan,main = "Inclusion reads, Spearman corr.")
pheatmap(cor(data.f[rowSums(data.f$e)>0,]$e, method = "pearson"),color = col.pan,main = "Exclusion reads, Spearman corr.")

df <- data.frame(colMeans(data.f$i), colMeans(data.f$e))
df$gr <- n
names(df) <- c("i", "e", "groups")
ggplot(data=df) + geom_point(aes(x = i, y = e, fill = groups, color = groups)) + theme_bw()

###Counts
p <- prcomp(t(gdata$cnts), center = TRUE)
ggbiplot(p, var.axes = F, groups = n, ellipse = T) + theme_bw() + ggtitle("Counts, PCA")
pheatmap(cor(gdata$cnts[rowSds(gdata$cnts)>1,], method = "pearson"), color = col.pan,main = "Inclusion ratio, Spearman corr.")

rpkm <- calcRPKM(data,gdata,F)
nrow(rpkm$rpkm)
rpkm <- rpkm[rowSums(rpkm$rpkm)>0,]$rpkm
rpkm[is.na(rpkm)] <- 0
pheatmap(cor(rpkm, method = "pearson"), color = col.pan,main = "Inclusion ratio, Spearman corr.")

###STATISTICAL TESTS

glm.function <- terms(x ~ transgenity*age)
data.f.glm = fitSAGLM(data.f,glm.function,meta)
data.f.pv = calcSAPvalue(data.f.glm)

data.f$seg <- cbind(data.f.pv, data.f$seg)
View(data.f$seg)

###pvalue.adjustment
data.f$seg$transgenity.adj <- p.adjust(data.f$seg$transgenity, method = "BH")
data.f$seg$age.adj <- p.adjust(data.f$seg$age, method = "BH")
data.f$seg$`transgenity:age.adj` <- p.adjust(data.f$seg$`transgenity:age`, method = "BH")


rn <- rownames(data.frame(data.f[which(data.f$seg$age.adj < 0.05),]))

#pval adjusted plots

plot(data.f$seg$transgenity,data.f$seg$transgenity.adj,xlab = "raw",ylab = "adjusted,BH",pch = ".")
lines(lowess(data.f$seg$transgenity.adj~data.f$seg$transgenity,f = .09), col = "red")

plot(data.f$seg$age,data.f$seg$age.adj,xlab = "raw",ylab = "adjusted,BH",pch = ".")
lines(lowess(data.f$seg$age.adj~data.f$seg$age,f = .09), col = "red")


make.stattest.and.annotate(sajr.obj = data.f,
              matrix = meta,
              stattest.id = "transgenity:age",
              control.id = "wt",
              case.id = "tg",
              psi.cutoff = 0.2,
              pval.cutoff = 0.05)


get_corr_per_group(data,gdata,meta$transgenity)
plot_ir("ENSMUSG00000024766.s10", meta$groups)
View(data.f$seg)
