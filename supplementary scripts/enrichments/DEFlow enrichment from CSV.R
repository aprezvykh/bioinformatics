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

dir <- c("~/counts/AIKAR.final/")
setwd(dir)
et_annot <- read.xlsx("filtered diffexpression ET_tt_LFC0,5.xlsx", sheetIndex = 2)



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
             species="cel", 
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

keg_up <- kegga(de = as.vector(et_annot_high$entrez), species="Ce", convert = TRUE)
tk_up <- topKEGG(keg_up, n=100)
tk_up <- subset(tk_up, tk_up$P.DE < 0.05)
tk_up$perc <- (tk_up$DE/tk_up$N)*100
tk_up <- tk_up[order(tk_up$perc, decreasing = TRUE),]
write.xlsx(tk_up, file = "kegga.xlsx", sheetName = "Upreg", append = TRUE)


keg_down <- kegga(de = as.vector(et_annot_low$entrez), species="Ce", convert = TRUE)
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

setwd("~/counts/AIKAR.final/early/")

all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez
GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)

allGO <- genesInTerm(GOdata)

s <- subset(go_up_500, go_up_500$perc > 20 & go_up_500$P.DE < 0.001)
s <- s[which(s$Ont == "BP"),]
s <- s[which(nchar(s$Term) < 50),]
s <- s[order(s$DE, decreasing = FALSE),]
s <- s[seq(1:30),]
s <- s[complete.cases(s),]
png(file = "Go terms upreg.png", width = 1024, height = 768)
#pdf(file = "Go terms upreg.pdf", width = 12, height = 17, family = "Helvetica")
g_u <- ggplot(s, aes(x = reorder(Term, DE), y = DE, fill = P.DE)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(breaks = s$Term, name = "Significant Upreguearlyd Terms") + 
  coord_flip() +
  scale_y_continuous(name = "Number of differentialy expressed genes") +
  theme_bw() + 
  theme(axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"))
g_u
dev.off()

s <- subset(go_down_500, go_down_500$perc > 20 & go_down_500$P.DE < 0.001)
s <- s[which(s$Ont == "BP"),]
s <- s[which(nchar(s$Term) < 50),]
s <- s[order(s$DE, decreasing = FALSE),]
s <- s[seq(1:30),]
s <- s[complete.cases(s),]
png(file = "Go terms downreg.png", width = 1024, height = 768)
#pdf(file = "Go terms downreg.pdf", width = 12, height = 17, family = "Helvetica")
g_d <- ggplot(s, aes(x = reorder(Term, DE), y = DE, fill = P.DE)) + 
  geom_bar(stat="identity") + 
  scale_x_discrete(breaks = s$Term, name = "Significant Downreguearlyd Terms") + 
  coord_flip() +
  scale_y_continuous(name = "Number of differentialy expressed genes") +
  theme_bw() + 
  theme(axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="plain"))
g_d
dev.off()

getwd()

### GO with foldchanges
go_fc <- goana(de = as.vector(et_annot$entrez), species = "Ce")
go_fc_res <- topGO(go_fc, n=1000, ontology = "BP")
df <- data.frame()
for (f in rownames(go_fc_res)){
  x <- allGO[[paste(f)]]
  z <- et_annot[(et_annot$entrez %in% x),]
  mean.logfc<- mean(z$logFC)
  mean.cpm <- mean(z$logCPM)
  b <- go_fc_res[which(rownames(go_fc_res) == f),]
  b$logFC <- mean.logfc
  b$cpm <- mean.cpm
  df <- rbind(b, df)
}
df <- df[complete.cases(df),]
write.xlsx(df, file = "goana tests with foldchanges.xlsx", sheetName = "Top 1000 GO terms, BP", append = TRUE)


go_fc <- goana(de = as.vector(et_annot$entrez), species = "Ce")
go_fc_res <- topGO(go_fc, n=1000, ontology = "MF")
df <- data.frame()
for (f in rownames(go_fc_res)){
  x <- allGO[[paste(f)]]
  z <- et_annot[(et_annot$entrez %in% x),]
  mean.logfc<- mean(z$logFC)
  mean.cpm <- mean(z$logCPM)
  b <- go_fc_res[which(rownames(go_fc_res) == f),]
  b$logFC <- mean.logfc
  b$cpm <- mean.cpm
  df <- rbind(b, df)
}
df <- df[complete.cases(df),]
write.xlsx(df, file = "goana tests with foldchanges.xlsx", sheetName = "Top 1000 GO terms, MF", append = TRUE)


go_fc <- goana(de = as.vector(et_annot$entrez), species = "Ce")
go_fc_res <- topGO(go_fc, n=1000, ontology = "CC")
df <- data.frame()
for (f in rownames(go_fc_res)){
  x <- allGO[[paste(f)]]
  z <- et_annot[(et_annot$entrez %in% x),]
  mean.logfc<- mean(z$logFC)
  mean.cpm <- mean(z$logCPM)
  b <- go_fc_res[which(rownames(go_fc_res) == f),]
  b$logFC <- mean.logfc
  b$cpm <- mean.cpm
  df <- rbind(b, df)
}
df <- df[complete.cases(df),]
write.xlsx(df, file = "goana tests with foldchanges.xlsx", sheetName = "Top 1000 GO terms, CC", append = TRUE)



### BP
all_genes = c(et_annot$logFC)
names(all_genes) <- et_annot$entrez

GOdata = new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.bp <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 30", append = TRUE)
write.xlsx(top100.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 100", append = TRUE)
write.xlsx(top500.res.bp, file = "topGO eliminations.xlsx", sheetName = "BP, top 500", append = TRUE)



dir.create("GO graphs")
setwd("GO graphs")

pdf(file = "top 5 GO graph, BP.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph, BP.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph, BP.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph, BP.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

setwd("~/counts/AIKAR.final/early/")

#MF

GOdata = new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.mf <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 30", append = TRUE)
write.xlsx(top100.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 100", append = TRUE)
write.xlsx(top500.res.mf, file = "topGO eliminations.xlsx", sheetName = "MF, top 500", append = TRUE)


setwd("GO graphs")

pdf(file = "top 5 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph, MF.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

setwd("~/counts/AIKAR.final/early/")


##CC

GOdata = new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(s) s < 
               0.05, description = "Test", annot = annFUN.org, mapping = "org.Ce.eg.db", nodeSize = 2)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

top30.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                         classicKS = resultKS, elimKS = resultKS.elim,
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

top100.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

top500.res.cc <- GenTable(GOdata, classicFisher = resultFisher,
                          classicKS = resultKS, elimKS = resultKS.elim,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 500)

write.xlsx(top30.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 30", append = TRUE)
write.xlsx(top100.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 100", append = TRUE)
write.xlsx(top500.res.cc, file = "topGO eliminations.xlsx", sheetName = "cc, top 500", append = TRUE)


setwd("GO graphs")

pdf(file = "top 5 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

pdf(file = "top 15 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 15, useInfo = 'all')
dev.off()

pdf(file = "top 30 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 30, useInfo = 'all')
dev.off()

pdf(file = "top 50 GO graph, cc.pdf", width = 10, height = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')
dev.off()

setwd("~/counts/AIKAR.final/")

