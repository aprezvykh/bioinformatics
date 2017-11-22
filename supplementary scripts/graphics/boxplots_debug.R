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

