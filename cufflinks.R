source("http://www.bioconductor.org/biocLite.R")
setwd("~/cufflinks/10_11_diffout.celegans/")
biocLite("cummeRbund")
install.packages("Rcpp")
install_version("RSQLite", version = "1.1-2", repos = "http://cran.us.r-project.org")
install_version("RSQLite", version = "2.0", repos = "http://cran.us.r-project.org")
install.packages("RSQLite")
library(devtools)
library(Rcpp)
require(devtools)
library("RSQLite")
library(cummeRbund)

cuff_data <- readCufflinks("~/cufflinks/10_11_diffout.celegans/", rebuild = TRUE)

csDensity(genes(cuff_data))
csScatter(genes(cuff_data), "10", "11")
csVolcano(genes(cuff_data), "10", "11")
mygene < - getGene(cuff_data, "regucalcin")
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))


