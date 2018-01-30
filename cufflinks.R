source("http://www.bioconductor.org/biocLite.R")
remove.packages("cummeRbund")
biocLite("cummeRbund")
install.packages("Rcpp")
install.packages("digest")
install_version("RSQLite", version = "1.1-2", repos = "http://cran.us.r-project.org")
install_version("RSQLite", version = "2.0", repos = "http://cran.us.r-project.org")
install_github("r-dbi/RSQLite")
install.packages("RSQLite")
install.packages("yaml")
library(devtools)
library(Rcpp)
library("RSQLite")
library(cummeRbund)
library(reshape2)
library(ggplot2)
library(plyr)
library(fastcluster)
library(rtracklayer)
library(Gviz)
library(BiocGenerics)

cuff_data<-readCufflinks(system.file("extdata", package="cummeRbund"), rebuild = T)
cuff_data 


csDensity(genes(cuff_data))
csScatter(genes(cuff_data), "10", "11")
csVolcano(genes(cuff_data), "10", "11")
mygene < - getGene(cuff_data, "regucalcin")
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))



