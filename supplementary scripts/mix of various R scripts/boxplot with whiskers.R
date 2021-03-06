library(ggplot2)
library(ggsignif)
library(org.Mm.eg.db)
library(GO.db)
library(gplots)
library(reshape)
setwd("~/counts/ALS Mice/")
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


setwd("~/counts/ALS Mice/experimental/results/all/lipids.new/")


library(xlsx)
lipids <- read.xlsx("lipids.xlsx", sheetIndex = 1)

for (f in lipids$genes){
print(f)
r <- grep(paste(f), rownames(cpm), ignore.case = TRUE)
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

g <- ggplot(thm, aes(x = Condition, y = gene)) + 
        geom_boxplot(aes(fill = Group), alpha = 0.5, coef = 1) + 
        stat_boxplot(geom = "errorbar", width = 0.5) + 
        geom_line(aes(x = Condition, y = mean(gene))) + 
        scale_x_discrete(name = "Experimental Groups") + 
        scale_y_continuous(name = "Log10(Counts per million)") + 
        theme_bw() + 
        geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
        geom_signif(comparisons = list(c("Tg-1", "Tg-2")), map_signif_level = TRUE) + 
        geom_signif(comparisons = list(c("Control-1", "Tg-1")), map_signif_level = TRUE) + 
  
  
        #geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) + 
        ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
        theme(plot.title = element_text(hjust = 0.5))

        ggsave(filename = paste(plot.name, "png", sep = "."), plot = g, height = 20, width = 20, units = "cm")

}




###dekta-FUS

#g <- ggplot(thm, aes(x = Condition, y = gene)) + 
#  geom_boxplot(aes(fill = Group), alpha = 0.5) + 
#  stat_boxplot(geom = "errorbar", width = 0.5) + 
# scale_x_discrete(name = "Experimental Groups") + 
#  scale_y_continuous(name = "Log10(Counts per million)") + 
#  theme_bw() + 
#  geom_signif(comparisons = list(c("Control-3", "Tg-1")), map_signif_level = TRUE) +
#  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) +
#  ggtitle(paste("Gene official symbol: ", "delta-FUS(1-359)", "\n", "Gene name:", "Fused in sarcoma", "\n", "Direct GO term:", "Nucleus")) + 
#  theme(plot.title = element_text(hjust = 0.5))
#g                     
#ggsave(filename = paste("delta-FUS", "pdf", sep = "."), plot = g)

x <- c("ENSMUSG00000079547", "ENSMUSG00000061232", "ENSMUSG00000041538", 
      "ENSMUSG00000060550", "ENSMUSG00000073409", "ENSMUSG00000053835")
thm <- cpm[(rownames(cpm) %in% x),]
rownames(thm) <- thm$Symbol
thm$Symbol <- NULL
thm$Name <- NULL
thm$GOID <- NULL
thm$term <- NULL
thm$entrez <- NULL
colnames(thm) <- sampleCondition
thm <-as.data.frame(t(thm))
thm$�ondition <- rownames(thm)
thm <- melt(thm, id.vars = "�ondition")

g <- ggplot(thm, aes(x = variable, y = value)) +
  geom_boxplot(data = thm, aes(fill = �ondition), alpha = 0.5) + 
  scale_x_discrete(name = "Experimental Groups") + 
  scale_y_continuous(name = "logCPM") + 
  theme_bw() +
  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) +
  #facet_wrap(~variable, scales="free") + 
  theme(axis.title.x=element_blank()) + 
  ggtitle("MHC class II genes") + theme(plot.title = element_text(hjust = 0.5, size = 10)) 
g
ggsave(filename = paste(i, "_", plot.name, ".pdf", sep = ""), plot = g, width = 20, height = 20, units = "cm")




#### BARPLOT w MEAN
x <- c("ENSMUSG00000031591")
thm <- cpm[(rownames(cpm) %in% x),]
rownames(thm) <- thm$Symbol
plot.name <- as.character(rownames(thm))
plot.description <- as.character(thm$Name)
plot.term <- as.character(thm$term)
thm$Symbol <- NULL
thm$Name <- NULL
thm$GOID <- NULL
thm$term <- NULL
thm$entrez <- NULL
#colnames(thm) <- sampleCondition
thm <-as.data.frame(t(thm))
thm$�ondition <- rownames(thm)
thm$group <- rownames(thm)
names(thm) <- c("int", "group")

g <- ggplot(thm) + geom_bar(aes(x = group, y = int, fill = group), stat = "identity", alpha= 0.5) + 
              scale_x_discrete(name = "Experimental Groups") + 
              scale_y_continuous(name = "logCPM") + 
              theme_bw() +
              #geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE) +
              theme(axis.title.x=element_blank()) + 
              ggtitle(paste("Gene official symbol: ", plot.name, "\n", "Gene name:", plot.description, "\n", "Direct GO term:", plot.term)) + 
              theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste(plot.name, "-separate", "pdf", sep = "."), plot = g, height = 20, width = 20, units = "cm")
              