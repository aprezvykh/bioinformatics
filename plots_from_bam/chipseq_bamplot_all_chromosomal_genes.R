#!/usr/local/bin/Rscript
library(here, quietly = TRUE)
library(GenomicAlignments, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
library(ggplot2, quietly = TRUE)

current.dir <- here()
args = commandArgs()
aln_file <- args[6]
scaffold <- args[7]
gtf <- args[8]
out_dir <- args[9]

#aln_file <- "~/Documents/ chip/194.226.21.15/Rezvykh.Chip/H3K9-1.SRR1533726_1.sorted.bam"
#scaffold <- "scaffold_13050"
#gtf <- "~/Documents/genomes/d.virilis/dvir-all-r1.06.gtf"
#out_dir <- "~/Scripts/out_test/"
plot_name <- system(paste("basename", aln_file), intern = TRUE)

a <- system(paste("cat ", gtf, " | grep ", scaffold, " | grep mRNA", sep = ""), intern = TRUE)
print(paste(length(a), "genes found!"), sep  = " ")

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
setwd(out_dir)
for (f in a){
    print(f)
    chr <- strsplit(f, "\t")[[1]][1]
    start <- as.numeric(strsplit(f, "\t")[[1]][4])
    end <- as.numeric(strsplit(f, "\t")[[1]][5])
    transcript <- strsplit(f, ";")[[1]][2]
    print(paste("Plotting ",gsub('\"', "",strsplit(transcript, "\"")[[1]][2], fixed = TRUE), "...", sep = ""))
    p <- plotCov(bampath=aln_file, 
            mychr=chr, 
            mystart=start, 
            myend=end, 
            ymin=-50, 
            ymax=50)
    ggsave(p, filename = paste(plot_name, "-",gsub('\"', "",strsplit(transcript, "\"")[[1]][2], fixed = TRUE), ".png", sep = ""))
    print(paste(plot_name, "plotted!", sep = " "))
}
getwd()
