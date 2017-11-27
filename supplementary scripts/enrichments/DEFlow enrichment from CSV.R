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

dir <- c("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/")
file <- grep("tg", list.files(dir), value = TRUE)
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

setwd(dir)

et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))
goana_up <- goana(de = et_annot_high$entrez, species = "Mm")
go_up_30 <- topGO(goana_up, n=30)
go_up_30$perc = (go_up_30$DE/go_up_30$N)*100
go_up_100 <- topGO(goana_up, n=100)
go_up_100$perc = (go_up_100$DE/go_up_100$N)*100
go_up_500 <- topGO(goana_up, n=500)
go_up_500$perc = (go_up_500$DE/go_up_500$N)*100
go_up_500$perc <- round(go_up_500$perc, digits = 4)

goana_down <- goana(de = et_annot_low$entrez, species = "Mm")
go_down_30 <- topGO(goana_down, n=30)
go_down_30$perc = (go_down_30$DE/go_down_30$N)*100
go_down_100 <- topGO(goana_down, n=100)
go_down_100$perc = (go_down_100$DE/go_down_100$N)*100
go_down_500 <- topGO(goana_down, n=500)
go_down_500$perc = (go_down_500$DE/go_down_500$N)*100
go_down_500$perc <- round(go_down_500$perc, digits = 4)




write.xlsx(go_up_30, file = "Goana GO tests, upreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_up_100, file = "Goana GO tests, upreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_up_500, file = "Goana GO tests, upreg.xlsx", sheetName = "top500", append = TRUE)

write.xlsx(go_down_30, file = "Goana GO tests, downreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_down_100, file = "Goana GO tests, downreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_down_500, file = "Goana GO tests, downreg.xlsx", sheetName = "top500", append = TRUE)

keg_com <- kegga(de = et_annot$entrez, species="Mm")
tk_common <- topKEGG(keg_com, n=100)
tk_common <- subset(tk_common, tk_common$P.DE < 0.05)
tk_common$perc <- (tk_common$DE/tk_common$N)*100
tk_common <- tk_common[order(tk_common$perc, decreasing = TRUE),]
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Mm")
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Mm")
tk_down <- topKEGG(keg_down, n=100)
tk_down <- subset(tk_down, tk_down$P.DE < 0.05)
tk_down$perc <- (tk_down$DE/tk_down$N)*100
tk_down <- tk_down[order(tk_down$perc, decreasing = TRUE),]
write.xlsx(tk_down, file = "kegga.xlsx", sheetName = "Downreg", append = TRUE)
getwd()
rownames(tk_common) <- substring(rownames(tk_common), 6)


dir.create("another kegg plots")
setwd("another kegg plots")


for (f in rownames(tk_common)){
  plot_pathway(f)
}


cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
cpm <- cpm[complete.cases(cpm),]
rownames(cpm) <- cpm$X
cpm$X <- NULL
cpm$Symbol <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="SYMBOL", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$Name <- mapIds(org.Mm.eg.db, 
                   keys=row.names(cpm), 
                   column="GENENAME", 
                   keytype="ENSEMBL",
                   multiVals="first")

cpm$entrez <- mapIds(org.Mm.eg.db, 
                     keys=row.names(cpm), 
                     column="ENTREZID", 
                     keytype="ENSEMBL",
                     multiVals="first")

cpm$GOID <-     mapIds(org.Mm.eg.db, 
                       keys=row.names(cpm), 
                       column="GO", 
                       keytype="ENSEMBL",
                       multiVals="first")

cpm$term <- mapIds(GO.db, 
                   keys=cpm$GOID, 
                   column="TERM", 
                   keytype="GOID",
                   multiVals="first")
cpm$term <- as.character(cpm$term)



setwd(dir)
dir.create("GO boxplots")
setwd("GO boxplots")
all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez
GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Mm.eg.db", nodeSize = 2)

allGO <- genesInTerm(GOdata)
goana_straight <- goana(de = et_annot$entrez, species = "Mm")
go_all_1000 <- topGO(goana_straight, n=1000)
go_all_1000 <- go_all_1000[!(go_all_1000$Ont %in% remove), ] 
go_all_1000$perc = (go_all_1000$DE/go_all_1000$N)*100
go_all_1000 <- subset(go_all_1000, go_all_1000$P.DE < 0.05)
go_all_1000 <- subset(go_all_1000, go_all_1000$perc > 20)
go_all_1000 <- subset(go_all_1000, go_all_1000$N < 100)
go_all_1000 <- subset(go_all_1000, go_all_1000$DE <= 10)
go_all_1000 <- subset(go_all_1000, go_all_1000$DE > 5)
go_all_1000 <- go_all_1000[order(go_all_1000$P.DE, decreasing = TRUE),]
go_all_1000 <- go_all_1000[seq(1:20),]
go_all_1000 <- go_all_1000[complete.cases(go_all_1000),]

for (i in 1:nrow(go_all_1000)){
  str <- go_all_1000[i,]
  f <- rownames(str)
  plot.name <- str$Term
  plot.ont <- str$Ont
  plot.perc <- round(str$perc, digits = 2)
  plot.pval <- round(str$P.DE, digits = 10)
  z <- allGO[[f]]
  z <- as.data.frame(z)
  thm <- cpm[(cpm$entrez %in% z$z),]
  rownames(thm) <- thm$Symbol
  thm$Symbol <- NULL
  thm$Name <- NULL
  thm$GOID <- NULL
  thm$term <- NULL
  thm$entrez <- NULL
  colnames(thm) <- sampleCondition
  thm <- log2(thm)
  thm <-as.data.frame(t(thm))
  thm$Ñondition <- rownames(thm)
  thm <- melt(thm, id.vars = "Ñondition")
  
  g <- ggplot(thm, aes(x = variable, y = value)) +
    geom_boxplot(data = thm, aes(fill = Ñondition), alpha = 0.5) + 
    scale_x_discrete(name = "Experimental Groups") + 
    scale_y_continuous(name = "logCPM") + 
    theme_bw() +
    geom_signif(comparisons = list(c("tg_2", "tg_3")), map_signif_level = TRUE) +
    #facet_wrap(~variable, scales="free") + 
    theme(axis.title.x=element_blank())
  g <- g + ggtitle(paste(f, ":", plot.name, "\n", "Ontology:", plot.ont, "\n", "Percent of significant genes: ", plot.perc, "%", "\n", "Term p-value:", plot.pval)) + theme(plot.title = element_text(hjust = 0.5, size = 10)) 
  ggsave(filename = paste(i, "_", plot.name, ".pdf", sep = ""), plot = g, width = 20, height = 20, units = "cm")
  
}  

setwd(dir)
getwd()
