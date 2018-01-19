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
library(edgeR)
gs_size <- 1
kegg_plots <- TRUE
data(go.sets.mm)
data(go.subs.mm)

library("dplyr")

dir <- c("~/counts/AIKAR/results/")
file <- grep("cel", list.files(dir), value = TRUE)
et_annot <- read.xlsx(paste(file), sheetIndex = 1)



rownames(et_annot) <- et_annot$NA.
et_annot_high <- subset(et_annot, et_annot$logFC > 0)
et_annot_low <- subset(et_annot, et_annot$logFC < 0)
foldchanges = et_annot$logFC
names(foldchanges) = et_annot$entrez

dir.create("early")
setwd("early")

dfa <- as.character(et_annot$entrez)
x <- enrichPathway(gene=dfa, organism = "celegans", minGSSize=1, readable = TRUE )
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
x <- enrichPathway(gene=df_high, organism = "celegans", minGSSize=1, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Up", append = TRUE)

par(mar=c(1,1,1,1))
pdf(file = "barplot_high.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()

#LOW
df_low <- et_annot_low$entrez
x <- enrichPathway(gene=df_low, organism = "celegans", minGSSize=1, readable = TRUE )
write.xlsx(x, "Reactome.xlsx", sheetName = "Down", append = TRUE)

par(mar=c(1,1,1,1))
pdf(file = "barplot_low.pdf", width = 12, height = 17, family = "Helvetica")
barplot(x, showCategory=30,  font.size = 9)
dev.off()



###KEGG expression profile (without lfc, but with generatio)
kk <- enrichKEGG(gene = et_annot_high$uniprot, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_upreg", append = TRUE)
pdf(file = "KEGG_upreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=3,  font.size = 9)
dev.off()

kk <- enrichKEGG(gene = et_annot_low$uniprot, organism = "cel", pvalueCutoff = 0.05, keyType = "uniprot")
write.xlsx(kk, file = "KEGG.xlsx", sheetName = "KEGG_downreg", append = TRUE)
pdf(file = "KEGG_downreg.pdf", width = 12, height = 17, family = "Helvetica")
barplot(kk, showCategory=3,  font.size = 9)
dev.off()


plot_pathway = function(pid){
    pathview(gene.data=foldchanges, 
             pathway.id=pid, 
             species="mmu", 
             new.signature=FALSE)
  }



et_annot_high <- as.data.frame(subset(et_annot, logFC > 0))
et_annot_low <- as.data.frame(subset(et_annot, logFC < 0))


goana_up <- goana(de = as.vector(et_annot_high$entrez), species = "Ce")
go_up_30 <- topGO(goana_up, n=30)
go_up_30$perc = (go_up_30$DE/go_up_30$N)*100
go_up_100 <- topGO(goana_up, n=100)
go_up_100$perc = (go_up_100$DE/go_up_100$N)*100
go_up_500 <- topGO(goana_up, n=500)
go_up_500$perc = (go_up_500$DE/go_up_500$N)*100
go_up_500$perc <- round(go_up_500$perc, digits = 4)

goana_down <- goana(de = as.vector(et_annot_low$entrez), species = "Ce")
go_down_30 <- topGO(goana_down, n=30)
go_down_30$perc = (go_down_30$DE/go_down_30$N)*100
go_down_100 <- topGO(goana_down, n=100)
go_down_100$perc = (go_down_100$DE/go_down_100$N)*100
go_down_500 <- topGO(goana_down, n=500)
go_down_500$perc = (go_down_500$DE/go_down_500$N)*100
go_down_500$perc <- round(go_down_500$perc, digits = 4)



getwd()
write.xlsx(go_up_30, file = "Goana GO tests, upreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_up_100, file = "Goana GO tests, upreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_up_500, file = "Goana GO tests, upreg.xlsx", sheetName = "top500", append = TRUE)

write.xlsx(go_down_30, file = "Goana GO tests, downreg.xlsx", sheetName = "top30", append = TRUE)
write.xlsx(go_down_100, file = "Goana GO tests, downreg.xlsx", sheetName = "top100", append = TRUE)
write.xlsx(go_down_500, file = "Goana GO tests, downreg.xlsx", sheetName = "top500", append = TRUE)






keg_com <- kegga(de = as.vector(et_annot$entrez), species="Ce", convert = TRUE)
tk_common <- topKEGG(keg_com, n=100)
tk_common <- subset(tk_common, tk_common$P.DE < 0.05)
tk_common$perc <- (tk_common$DE/tk_common$N)*100
tk_common <- tk_common[order(tk_common$perc, decreasing = TRUE),]
write.xlsx(tk_common, file = "kegga.xlsx", sheetName = "all Kegg", append = TRUE)

keg_up <- kegga(de = et_annot_high$entrez, species="Ce", convert = TRUE)
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = et_annot_low$entrez, species="Ce", convert = TRUE)
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

