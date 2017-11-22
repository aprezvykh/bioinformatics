x <- read.xlsx("genes for boxplots.xlsx", sheetIndex = 2)

for (f in x$Glia){
  src <- c("Microglia")
  r <- grep(paste(f), rownames(cpm), ignore.case = TRUE)
  thm <- cpm[r,]
  rownames(thm) <- thm$Symbol
  plot.name <- rownames(thm)
  plot.description <- thm$Name
  plot.term <- thm$term
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
    geom_boxplot(data = thm, aes(fill = Group), alpha = 0.5) + 
    scale_x_discrete(name = "Experimental Groups") + 
    scale_y_continuous(name = "Counts per million") + 
    theme_bw() + 
    geom_signif(comparisons = list(c("Tg-2", "Tg-3"), c("Tg-1", "Tg-2")), map_signif_level = TRUE)
  
  
  g + ggtitle(paste("Gene official symbol: ", plot.name, "\n" , "Gene name: ", plot.description, "\n", "Tissue: ", src, "\n", "GO Direct term: ", plot.term)) + theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste(plot.name, "-", src, ".pdf", sep = ""), plot = g)
  
}


# "Tissue: ", src, "\n", 


###TRANSGENIC FUS
#(paster after plot.description)
thm <- read.xlsx("dFUS.xlsx", sheetIndex = 1)
thm$Condition <- sampleCondition
pdf(file = "delta-FUS.pdf", width = 10, height = 10, family = "Helvetica")
g <- ggplot(thm, aes(x = Condition, y = gene)) + 
  geom_boxplot(data = thm, aes(fill = Condition), alpha = 0.5) + 
  scale_x_discrete(name = "Experimental Groups") + 
  scale_y_continuous(name = "Reads per million") + 
  theme_bw() +
  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE, annotations = "***")

g + ggtitle(paste("Gene official symbol: ", "delta-FUS","\n", "Gene name: ", plot.description, "\n", "Tissue: ", src, "\n", "GO Direct term: ", plot.term)) + theme(plot.title = element_text(hjust = 0.5)) 
dev.off()


####MORE THAN TWO GENES IN PLOT
b <- read.xlsx("~/counts/ALS Mice/experimental/results/fc1_w_kegga/tg_2-tg_3/Results edgeR.xlsx", sheetIndex = 2)
f <- b$NA.
a <- grepl(paste(f, collapse = "|"), rownames(cpm), ignore.case = TRUE)
thm <- cpm[a,]
rownames(thm) <- thm$Symbol
thm$Symbol <- NULL
thm$Name <- NULL
thm$GOID <- NULL
thm$term <- NULL
colnames(thm) <- sampleTable$condition
thm <- log2(thm)
thm <-as.data.frame(t(thm))
thm$Ñondition <- rownames(thm)
thm <- melt(thm, id.vars = "Ñondition")
pdf(file = "Top Tg2-Tg3 genes.pdf", width = 10, height = 10, family = "Helvetica")
g <- ggplot(thm, aes(x = Ñondition, y = value)) + geom_boxplot() +
  geom_boxplot(data = thm, aes(fill = Ñondition), alpha = 0.5) + 
  scale_x_discrete(name = "Experimental Groups") + 
  scale_y_continuous(name = "logCPM") + 
  theme_bw() +
  geom_signif(comparisons = list(c("Tg-2", "Tg-3")), map_signif_level = TRUE, annotations = "***") +
  facet_wrap(~variable, scales="free") +  
  theme(axis.title.x=element_blank())
g + ggtitle("Top Tg2-Tg3 genes") + theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

