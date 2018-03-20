### DISEASE ASSOCIATION

#if (disease_association == TRUE){
et_annot_high <- as.data.frame(subset(et_annot, logFC > logfchigh_cutoff))
et_annot_low <- as.data.frame(subset(et_annot, logFC < logfclow_cutoff))

setwd(results.dir)
if (analyze_all_samples == TRUE){
  setwd("all")
} else {
  setwd(stattest)
}

table <- read.delim(file = '~/GitHub/data/curated_gene_disease_associations.tsv')
table <- as.data.frame(table)

diseases_up <- data.frame()
diseases_down <- data.frame()
for (f in et_annot_high$symbol){
  sub <- NULL
  sub <- table[grepl(paste(f), table$geneSymbol, ignore.case = TRUE),]
  sub <- sub[order(sub$score, decreasing = TRUE),]
  sub <- sub[seq(1:5),]
  sub <- as.data.frame(sub$diseaseName)
  sub <- transpose(sub)
  sub$gene <- paste(f)
  diseases_up <- rbind(sub, diseases_up)
}
# write.xlsx(diseases_up, file = "Diseases association by Disgenet.xlsx", sheetName = "upreg", append = TRUE)

for (f in et_annot_low$symbol){
  sub <- NULL
  sub <- table[grepl(paste(f), table$geneSymbol, ignore.case = TRUE),]
  sub <- sub[order(sub$score, decreasing = TRUE),]
  sub <- sub[seq(1:5),]
  sub <- as.data.frame(sub$diseaseName)
  sub <- transpose(sub)
  sub$gene <- paste(f)
  diseases_down<- rbind(sub, diseases_down)
}
# write.xlsx(diseases_down, file = "Diseases association by Disgenet.xlsx", sheetName = "downreg", append = TRUE)


up <- as.data.frame(table(unlist(diseases_up)))
up <- up[order(up$Freq, decreasing = TRUE),]
up <- up[seq(1:diseases_set),]
names(up) <- c("Disease", "Frequency")
write.xlsx(up, file = "Top Diseases by Disgenet.xlsx", sheetName = "upreg", append = TRUE)

down <- as.data.frame(table(unlist(diseases_down)))
down <- down[order(down$Freq, decreasing = TRUE),]
down <- down[seq(1:diseases_set),]
names(down) <- c("Disease", "Frequency")
write.xlsx(down, file = "Top Diseases by Disgenet.xlsx", sheetName = "downreg", append = TRUE)
}


### PANTHER.DB

#if (panther_analysis == TRUE){
pan_up <- et_annot_high
pan_down <- et_annot_low

pan_up$goslim <- mapIds(PANTHER.db, 
                        keys=et_annot_high$entrez, 
                        column="GOSLIM_ID", 
                        keytype="ENTREZ",
                        multiVals="first")

pan_up$pathway <- mapIds(PANTHER.db, 
                         keys=et_annot_high$entrez, 
                         column="PATHWAY_ID", 
                         keytype="ENTREZ",
                         multiVals="first")
pan_down$goslim <- mapIds(PANTHER.db, 
                          keys=et_annot_low$entrez, 
                          column="GOSLIM_ID", 
                          keytype="ENTREZ",
                          multiVals="first")

pan_down$pathway <- mapIds(PANTHER.db, 
                           keys=et_annot_low$entrez, 
                           column="PATHWAY_ID", 
                           keytype="ENTREZ",
                           multiVals="first")


go_pan_up <- as.data.frame(table(unlist(pan_up$goslim)))
go_pan_up <- go_pan_up[order(go_pan_up$Freq, decreasing = TRUE),]
go_pan_up <- go_pan_up[seq(1:go_terms_set),]
names(go_pan_up) <- c("GO term", "Frequency")

go_pan_down <- as.data.frame(table(unlist(pan_down$goslim)))
go_pan_down <- go_pan_down[order(go_pan_down$Freq, decreasing = TRUE),]
go_pan_down <- go_pan_down[seq(1:go_terms_set),]
names(go_pan_down) <- c("GO term", "Frequency")

pth_pan_up <- as.data.frame(table(unlist(pan_up$pathway)))
pth_pan_up <- pth_pan_up[order(pth_pan_up$Freq, decreasing = TRUE),]
pth_pan_up <- pth_pan_up[seq(1:pathways_set),]
names(pth_pan_up) <- c("Pathway ID", "Frequency")

pth_pan_down <- pth_pan_down[order(pth_pan_down$Freq, decreasing = TRUE),]
pth_pan_down <- pth_pan_down[seq(1:pathways_set),]
names(pth_pan_down) <- c("Pathway ID", "Frequency")

write.xlsx(go_pan_up, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(go_pan_down, file = "GOSlim Terms by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)

write.xlsx(pth_pan_up, file = "Top Pathways by PANTHER.xlsx", sheetName = "UP", append = TRUE)
write.xlsx(pth_pan_down, file = "Top Pathways by PANTHER.xlsx", sheetName = "DOWN", append = TRUE)
}
