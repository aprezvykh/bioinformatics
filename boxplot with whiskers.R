library(ggsignif)
sampleCondition <- c('Control-1', 'Control-1', 'Control-1', 'Control-1', 'Control-1', 
                     'Control-3', 'Control-3', 'Control-3', 'Control-3', 'Control-3', 
                     'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 'Tg-1', 
                     'Tg-2', 'Tg-2', 'Tg-2', 'Tg-2',
                     'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3', 'Tg-3')
                    
row.names.remove <- c("NA.1")
col.pan <- colorpanel(100, "blue", "white", "red")
cpm <- read.csv("~/counts/ALS Mice/experimental/results/all/overall logCPM.csv")
cpm <- cpm[!(row.names(cpm) %in% row.names.remove), ] 
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

x <- read.xlsx("~/counts/ALS Mice/experimental/results/FUS.xlsx", sheetIndex = 2)
r <- grep("ENSMUSG00000030795", rownames(cpm), ignore.case = TRUE)
thm <- cpm[r,]
rownames(thm) <- thm$Symbol
plot.name <- as.character(rownames(thm))
plot.description <- as.character(thm$Name)
plot.term <- as.character(thm$term)
thm$Symbol <- NULL
thm$Name <- NULL
thm$GOID <- NULL
thm$term <- NULL
thm$entrez <- NULL
colnames(thm) <- sampleCondition
thm <-as.data.frame(t(thm))
thm$Condition <- rownames(thm)
thm$Group <- thm$Condition
names(thm) <- c("gene", "Condition", "Group")
names(thm) <- c("Condition", "Group", "gene")

g <- ggplot(thm, aes(x = Condition, y = gene)) + 
        geom_boxplot(aes(fill = Group), alpha = 0.5, coef = 1) + 
        stat_boxplot(geom = "errorbar", width = 0.5) + 
        geom_line(aes(x = Condition, y = mean(gene))) + 
        scale_x_discrete(name = "Experimental Groups") + 
        scale_y_continuous(name = "Log10(Counts per million)") + 
        theme_bw() + 
        geom_signif(comparisons = list(c("Tg-1", "Tg-3")), map_signif_level = TRUE) +
        ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        geom_line(aes(y = Condition))
g                     

ggsave(filename = paste(plot.name, "pdf", sep = "."), plot = g)


g <- ggplot(thm, aes(x = Condition, y = gene)) + 
  geom_boxplot(aes(fill = Group), alpha = 0.5) + 
  stat_boxplot(geom = "errorbar", width = 0.5) + 
  scale_x_discrete(name = "Experimental Groups") + 
  scale_y_continuous(name = "Log10(Counts per million)") + 
  theme_bw() + 
  geom_signif(comparisons = list(c("Control-3", "Tg-1")), map_signif_level = TRUE) +
  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) +
  ggtitle(paste("Gene official symbol: ", "delta-FUS(1-359)", "\n", "Gene name:", "Fused in sarcoma", "\n", "Direct GO term:", "Nucleus")) + 
  theme(plot.title = element_text(hjust = 0.5))
g                     
ggsave(filename = paste("delta-FUS", "pdf", sep = "."), plot = g)
