source("https://bioconductor.org/biocLite.R")
biocLite()
library(devtools)
devtools::install_github("rstats-db/RSQLite")
devtools::install_github("rstats-db/RSQLite@b-%23223-icc")

library(cummeRbund)

cuff_data<-readCufflinks(system.file("extdata", package="cummeRbund", lib.loc = .libPaths()[1]), rebuild = TRUE)
cuff_data 


csDensity(genes(cuff_data))
csScatter(genes(cuff_data), "10", "11")
csVolcano(genes(cuff_data), "10", "11")
mygene < - getGene(cuff_data, "regucalcin")
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))



