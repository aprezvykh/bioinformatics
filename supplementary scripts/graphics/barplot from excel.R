library(ggplot2)
library(xlsx)
library(ggsignif)
a <- read.xlsx("~/coomassie+phos.xlsx", sheetIndex = 4)
colnames(a) <- c("group", "number")
a$group <- c("1 - Control", "2 - A 42", "3 - isoA 42")
a$err <- c(1.3, 2.3, 3.2)
a$number<- a$number*100


pdf("Coomassie.pdf", width = 8, height = 8)
ggplot(a, aes(x = group, y = number)) + 
       geom_bar(stat = "identity", width = 0.2, fill = "black", alpha = 0.5) +
       theme_bw() + 
       scale_y_continuous(name = "Relative matrin-3 Pro-Q Diamond/Matrin-3 Coomassie, %", limits = c(0, 105), expand = c(0, 0)) + 
       scale_x_discrete(name = "Expreimental Groups") +
       geom_signif(comparisons = list(c("1 - Control", "3 - isoA 42")), map_signif_level = FALSE, annotations = "*", y_position = 102, size = 0.5) + 
       ggtitle("") + 
       theme(plot.title = element_text(hjust = 0.5)) + 
       theme(axis.text.x=element_text(colour="black")) + 
       theme(plot.title = element_text(hjust = 0.5)) + 
       theme(axis.text.y=element_text(colour="black")) +
       geom_errorbar(aes(ymin = number-err, ymax = number+err),width=0.1, size=0.5) + 
       coord_cartesian(ylim=c(90,105))

dev.off()        

library(beepr)
while(TRUE){
  beep()
  Sys.sleep(0.05)
}