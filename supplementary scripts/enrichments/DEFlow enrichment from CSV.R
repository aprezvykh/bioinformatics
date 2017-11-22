## GO & REACTOME ENRICHMENT FROM CSV FILE.
## NOT WORKING WITH DROSOPHILA
library("gage")
library("gageData")
library("pathview")
library("dplyr")
library("ReactomePA")
library("reactome.db")
library("topGO")
library("clusterProfiler")
library("xlsx")
library("AnnotationDbi")
library("org.Mm.eg.db")
library("GO.db")
gs_size <- 1
kegg_plots <- TRUE
data(go.sets.mm)
data(go.subs.mm)

library("dplyr")

dir <- c("~/GitHub/counts/ALS Mice/tg1-tg2/other/")
file <- grep("deg", list.files(dir), value = TRUE)
setwd(dir)
et_annot <- read.csv(paste(file))

rownames(et_annot) <- et_annot$X
et_annot_high <- subset(et_annot, et_annot$logFC > 0)
et_annot_low <- subset(et_annot, et_annot$logFC < 0)
foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

gomfsets = go.sets.mm[go.subs.mm$MF]
gomfres = gage(foldchanges, gsets=gomfsets, same.dir=TRUE,set.size = c(5, 100), rank.test = TRUE)
lapply(gomfres, head)
gomfres <- as.data.frame(gomfres)
gomfres <- gomfres[complete.cases(gomfres), ]
write.xlsx(gomfres, file = "GO.xlsx", sheetName = "GO_MF", append = TRUE)

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
gobpres <- as.data.frame(gobpres)
gobpres <- gobpres[complete.cases(gobpres), ]
write.xlsx(gobpres, file = "GO.xlsx", sheetName = "GO_BP", append = TRUE)

goccsets = go.sets.mm[go.subs.mm$CC]
goccres = gage(foldchanges, gsets=goccsets, same.dir=TRUE)
lapply(goccres, head)
goccres <- as.data.frame(goccres)
goccres <- goccres[complete.cases(goccres), ]
write.xlsx(goccres, file = "GO.xlsx", sheetName = "GO_CC", append = TRUE)


dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "mouse", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "All reactome", append = TRUE)
head(as.data.frame(x))
dev.off()
dev.off()

par(mar=c(1,1,1,1))
pdf(file = "barplot.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

pdf(file = "enrichmap.pdf", width = 12, height = 17, family = "Helvetica")
enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7, n = 20, font.size = 20)
dev.off()

pdf(file = "cnetplot.pdf", width = 12, height = 17, family = "Helvetica")
cnetplot(x, foldChange = foldchanges, categorySize="pvalue", showCategory = 20)
dev.off()

#HIGH
df_high <- et_annot_high$entrez
x <- enrichPathway(gene=df_high, organism = "mouse", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Up", append = TRUE)

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW

df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "mouse", minGSSize=gs_size, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Down", append = TRUE)

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()



###KEGG expression profile (without lfc, but with generatio)

kk <- enrichKEGG(gene = df_high, organism = "mmu", pvalueCutoff = 0.2)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=3,  font.size = 9)
dev.off()

kk <- enrichKEGG(gene = df_low, organism = "mmu", pvalueCutoff = 0.2)
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=3,  font.size = 9)
dev.off()


data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]

keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggreswr <- as.data.frame(keggres)
write.xlsx(keggreswr, file = "KEGG pathview.xlsx", sheetName = "KEGG", append = TRUE)

if (kegg_plots == TRUE){
  keggrespathways = data.frame(id=rownames(keggres$greater),
                               keggres$greater) %>% 
    tbl_df() %>% 
    filter(row_number()<=20) %>% 
    .$id %>% 
    as.character()
  keggresids = substr(keggrespathways, start=1, stop=8)
  
  plot_pathway = function(pid){
    pathview(gene.data=foldchanges, 
             pathway.id=pid, 
             species="mmu", 
             new.signature=FALSE)
  }
  detach("package:dplyr", unload=TRUE)
  dir.create("kegg")
  setwd("kegg")
  tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
}



