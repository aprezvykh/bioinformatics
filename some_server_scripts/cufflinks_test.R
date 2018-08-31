source("https://bioconductor.org/biocLite.R")
biocLite("ballgown")
library(ballgown)
## ----loadmethods, echo=FALSE, message=FALSE, warning=FALSE---------------
library(methods)

## ----installme, eval=FALSE-----------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("ballgown")

## ----makebgobj, message=FALSE--------------------------------------------
library(ballgown)
data_directory = system.file('extdata', package='ballgown') # automatically finds ballgown's installation directory
# examine data_directory:
data_directory

# make the ballgown object:
bg = ballgown(dataDir=data_directory, samplePattern='sample', meas='all')


## ----struct--------------------------------------------------------------
structure(bg)$exon
structure(bg)$intron
structure(bg)$trans

## ----getexpr-------------------------------------------------------------
transcript_fpkm = texpr(bg, 'FPKM')
transcript_cov = texpr(bg, 'cov')
whole_tx_table = texpr(bg, 'all')
exon_mcov = eexpr(bg, 'mcov')
junction_rcount = iexpr(bg)
whole_intron_table = iexpr(bg, 'all')
gene_expression = gexpr(bg)

## ----pData---------------------------------------------------------------
pData(bg) = data.frame(id=sampleNames(bg), group=rep(c(1,0), each=10))

## ----indexex-------------------------------------------------------------
exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)
phenotype_table = pData(bg)

plotTranscripts(gene='XLOC_000454', gown=bg, samples='sample12', 
                meas='FPKM', colorby='transcript', 
                main='transcripts from gene XLOC_000454: sample 12, FPKM')

## ----plotTranscripts2, fig.height=10, fig.width=10, fig.cap=""-----------
plotTranscripts('XLOC_000454', bg, 
                samples=c('sample01', 'sample06', 'sample12', 'sample19'), 
                meas='FPKM', colorby='transcript')

## ----plotMeans, fig.height=5, fig.width=10, fig.cap=""-------------------
plotMeans('XLOC_000454', bg, groupvar='group', meas='FPKM', colorby='transcript')

## ----stattest------------------------------------------------------------
stat_results = stattest(bg, feature='transcript', meas='FPKM', covariate='group')
head(stat_results)

## ----stattest_time-------------------------------------------------------
pData(bg) = data.frame(pData(bg), time=rep(1:10, 2)) #dummy time covariate
timecourse_results = stattest(bg, feature='transcript', meas='FPKM', covariate='time', timecourse=TRUE)

## ----stattest_adjust-----------------------------------------------------
group_adj_timecourse_results = stattest(bg, feature='transcript', meas='FPKM', covariate='time', 
                                        timecourse=TRUE, adjustvars='group')

## ----customMod-----------------------------------------------------------
# create example data:
set.seed(43)
sex = sample(c('M','F'), size=nrow(pData(bg)), replace=TRUE)
age = sample(21:52, size=nrow(pData(bg)), replace=TRUE)

# create design matrices:
mod = model.matrix(~ sex + age + pData(bg)$group + pData(bg)$time)
mod0 = model.matrix(~ pData(bg)$group + pData(bg)$time)

# run differential expression tests:
adjusted_results = stattest(bg, feature='transcript', meas='FPKM', mod0=mod0, mod=mod)
head(adjusted_results)

## ----cluster-------------------------------------------------------------
clusterTranscripts(gene='XLOC_000454', gown=bg, k=2, method='kmeans')

## ----clusterviz, fig.height=8, fig.width=8, fig.cap=""-------------------
plotLatentTranscripts(gene='XLOC_000454', gown=bg, k=2, method='kmeans', returncluster=FALSE)

## ----collapse------------------------------------------------------------
agg = collapseTranscripts(gene='XLOC_000454', gown=bg, k=2, method='kmeans')
stattest(gowntable=agg$tab, pData=pData(bg), feature='transcript_cluster', 
         covariate='group', libadjust=FALSE)

## ----sessioninfo, results='markup'---------------------------------------
sessionInfo()
