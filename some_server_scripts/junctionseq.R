install.packages("~/plotrix_3.7-1.tar.gz", repos = NULL, type="source")

source("https://bioconductor.org/biocLite.R")
biocLite("JunctionSeq")
library("JunctionSeq")
install.packages("http://hartleys.github.io/JunctionSeq/install/JctSeqData_LATEST.tar.gz", repos = NULL, type = "source")


decoder.file <- system.file("extdata/annoFiles/decoder.bySample.txt",
                            package="JctSeqData",
                            mustWork=TRUE);

decoder <- read.table(decoder.file,
                      header=TRUE,
                      stringsAsFactors=FALSE);

gff.file <- system.file("extdata/cts/withNovel.forJunctionSeq.gff.gz",
                      package="JctSeqData",
                      mustWork=TRUE)

countFiles.noNovel <- system.file(paste0("extdata/cts/",
                                         decoder$sample.ID,
                                         "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt.gz"),
                                  package="JctSeqData", mustWork=TRUE)

countFiles <- system.file(paste0("extdata/cts/",
                                 decoder$sample.ID,
                                 "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
                          package="JctSeqData", mustWork=TRUE)

gff.file <- system.file(
  "extdata/tiny/withNovel.forJunctionSeq.gff.gz",
  package="JctSeqData",
  mustWork=TRUE)

countFiles <- system.file(paste0("extdata/tiny/",
                                 decoder$sample.ID,
                                 "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"),
                          package="JctSeqData", mustWork=TRUE)


jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 32,
                               analysis.type = "junctionsAndExons")

design <- data.frame(condition = factor(decoder$group.ID))

geneID.to.symbol.file <- system.file(
                  "extdata/annoFiles/ensid.2.symbol.txt",
                  package="JctSeqData",
                  mustWork=TRUE)

jscs = readJunctionSeqCounts(countfiles = countFiles,
                             samplenames = decoder$sample.ID,
                             design = design,
                             flat.gff.file = gff.file,
                             gene.names = geneID.to.symbol.file)


jscs <- estimateJunctionSeqSizeFactors(jscs)
jscs <- estimateJunctionSeqDispersions(jscs, nCores = 8)
jscs <- fitJunctionSeqDispersionFunction(jscs)
jscs <- testForDiffUsage(jscs, nCores = 8)
jscs <- estimateEffectSizes( jscs, nCores = 8)
buildAllPlots(jscs=jscs,
              outfile.prefix = "./plots/",
              use.plotting.device = "png",
              FDR.threshold = 0.01
)
plotDispEsts(jscs)
