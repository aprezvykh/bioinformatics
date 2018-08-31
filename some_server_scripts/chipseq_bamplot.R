#!/usr/bin/Rscript
library(here, quietly = TRUE)
library(GenomicAlignments, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
library(ggplot2, quietly = TRUE)

current.dir <- here()
args = commandArgs()
aln_file <- args[6]
transcript <- args[7]
gtf <- args[8]

plot_name <- system(paste("basename", aln_file), intern = TRUE)

a <- system(paste("cat ", gtf, " | grep ", transcript, " | grep mRNA", sep = ""), intern = TRUE)
chr <- strsplit(a, "\t")[[1]][1]
start <- as.numeric(strsplit(a, "\t")[[1]][4])
end <- as.numeric(strsplit(a, "\t")[[1]][5])
plotCov <- function(bampath, mychr, mystart, myend, ymin, ymax) {
  ga <- readGAlignments(bampath, use.names=TRUE, param=ScanBamParam(which=GRanges(mychr, IRanges(mystart, myend))))
  pos <- as.numeric(coverage(ga[strand(ga)=="+"])[[mychr]])[mystart:myend]
  pos <- data.frame(Chr=rep(mychr, length(pos)), Strand=rep("+", length(pos)), Position=mystart:myend, Coverage=pos)
  neg <- as.numeric(coverage(ga[strand(ga)=="-"])[[mychr]])[mystart:myend]
  neg <- data.frame(Chr=rep(mychr, length(neg)), Strand=rep("-", length(neg)), Position=mystart:myend, Coverage=-neg)
  covdf <- rbind(pos, neg)
  ggplot(covdf, aes(Position, Coverage, fill=Strand)) + 
    geom_bar(stat="identity", position="identity") + 
    ylim(ymin,ymax) + 
    ggtitle(paste(plot_name, "\n",
                  transcript, "\n",
                  "Coordinates: ", chr, ":", start, "-", end, sep = "")) + 
    theme_bw()
}

png(paste(plot_name, "-",transcript, ".png", sep = ""))
paste(plot_name, "-",transcript, ".png", sep = "")
plotCov(bampath=aln_file, 
        mychr=chr, 
        mystart=start, 
        myend=end, 
        ymin=-50, 
        ymax=50)
dev.off()
