library(RcisTarget)
library(RcisTarget.mm9.motifDatabases.20k)
library(DT)
library(reshape2)
library(htmlwidgets)
library(visNetwork)
install.packages("visNetwork")
data("mm9_10kbpAroundTss_motifRanking")
data("mm9_direct_motifAnnotation")
dir <- c("~/GitHub/counts/ALS Mice/new filtering-0.5/tg1-tg2/glia/")
file <- grep("deg", list.files(dir), value = TRUE)
setwd(dir)
et_annot <- read.csv(paste(file))
rownames(et_annot) <- et_annot$NA.

geneList1 <- et_annot$symbol
geneLists <- list(geneListName=geneList1)
motifRankings <- mm9_10kbpAroundTss_motifRanking

motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                         motifAnnot_direct=mm9_direct_motifAnnotation)
motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo
result <- datatable(resultsSubset[,-c("enrichedGenes", "TF_inferred"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))


saveWidget(result, file="TF enrichment tg12-glia.html")

motifs_AUC <- calcAUC(geneLists, motifRankings, nCores=1)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, motifAnnot_direct=mm9_direct_motifAnnotation)
signifMotifNames <- motifEnrichmentTable$motif[1:5]
incidenceMatrix <- getSignificantGenes(geneLists, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")


motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))

web <- visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)

saveWidget(web, file="TF enrichment web.html")
