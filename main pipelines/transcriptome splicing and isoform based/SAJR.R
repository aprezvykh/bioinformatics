#!/usr/bin/Rscript
setwd("~/transcriptomes/reads/intron_retention/sod_moto/sajr/")
library(SAJR)
library(edgeR)
#pdf('output.pdf')
data = loadSAData(ann.gff='a.gff',c(1,2,3,4))
data = setSplSiteTypes(data,'a.gff')
table(data$seg$sites)
data.f = data[data$seg$type %in% c('ALT','INT') & data$seg$position %in% c('LAST','INTERNAL','FIRST') & apply(data$i+data$e>=10,1,sum)==2 & apply(data$ir,1,sd,na.rm=TRUE) > 0,]
# split genes into alternatives
all.alts = makeAlts(data$seg,'a.gff',remove.exn.ext = F)
par(mfrow=c(4,4))
plotAllAlts(all.alts,to.plot = 16)
no.int.ret.alt = filterAlts(all.alts,TRUE)
plotAllAlts(no.int.ret.alt,to.plot = 16)
par(mfrow=c(1,1))
plotAllAlts(all.alts, min.cnt = TRUE)
alt.summary(all.alts)


mod = list(f=factor(c('control','control','treatment','treatment')))
data.f.glm = fitSAGLM(data.f,terms(x ~ f),mod,overdisp = F)
df <- data.frame(data.f)
# calculate p-value. Since there are no replicates, use binomial distribution
# newer use 'overdisp = FALSE' in real work.
data.f.pv = calcSAPvalue(data.f.glm)
data.f.pv[1:10,]
#plot histogramm of p-values
hist(data.f.pv[,2])
data.f.pv[1:10,]
#make BH correction
data.f.pv[,2] = p.adjust(data.f.pv[,2],method='BH')
#choose significant ones:
data.sign = data.f[data.f.pv[,2] <=0.05,]
length(data.sign)
# print top ten (by amplitude)
data.sign[order(abs(data.sign$ir[,1]-data.sign$ir[,2]),decreasing = TRUE),][1:10,]

df <- data.frame(data.sign)
df <- df[complete.cases(df),]


df$gene<- mapIds(org.Mm.eg.db, 
                  keys=df$seg.gene_id, 
                  column="GENENAME", 
                  keytype="ENSEMBL",
                  multiVals="first")

write.csv(df, "results/sign_splice.csv")

dir.create("alt_splice_png")
setwd("alt_splice_png")
for (f in df$seg.gene_id){
  if(is.na(f)){
    next
  }  else {
    png(paste(f,".png",sep = ""))
    plotAlt(segs = all.alts[which(all.alts$gene_id == paste(f)),]$segs,
            all.alts[which(all.alts$gene_id == paste(f)),]$ints,
            col.exn = "red",
            col.int = "green",
            col.alt = "blue",
            col.junc = "green",
            main = paste(f,"\n",df[which(df$seg.gene_id == f),]$gene))
    dev.off()
  }
}

# check selfonsistency
setwd("../results")
plotCorrHM(data$ir)
plotMDS(data$ir)
plotMDS(data$ir[data$seg$sites == 'ad',],pch=rep(c(7,19),each=4),col=rep(c('red','blue'),each=5),main='Cassette exons')

