library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb <- keepSeqlevels(txdb, "chr16")
seqlevelsStyle(txdb) <- "NCBI"



txf_ucsc <- convertToTxFeatures(txdb)
txf_ucsc <- txf_ucsc[txf_ucsc %over% gr]

sgf_ucsc <- convertToSGFeatures(txf_ucsc)

path <- system.file("extdata", package = "SGSeq")
si$file_bam <- file.path(path, "bams", si$file_bam)

sgfc_ucsc <- analyzeFeatures(si, features = txf_ucsc)


df <- plotFeatures(sgfc_ucsc, geneID = 1)
sgfc_pred <- analyzeFeatures(si, which = gr)
head(rowRanges(sgfc_pred))
