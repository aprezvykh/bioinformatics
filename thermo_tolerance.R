library(xlsx)
library(ggplot2)
library(ggbiplot)
library(cowplot)
library(ggsignif)
df <- read.xlsx("~/thermotol.xlsx", sheetIndex = 1)
names(df) <- c("Experimental group", "Gender", "39°С", "39.5°С", "sd1", "sd2")
g1 <- ggplot(df, aes(x = `Experimental group`, y = `39°С`, fill = Gender, group = Gender)) + 
             geom_bar(stat = "identity", alpha = 0.7) + 
             geom_errorbar(aes(ymin=`39°С`-sd1, ymax=`39°С`+sd1)) +
             theme_bw() + 
             ggtitle("39°С Thermotolerance") + 
             theme(plot.title = element_text(hjust = 0.5)) + 
             scale_y_continuous("% of alive", limits = c(-3,130)) + 
             geom_signif(comparisons = list(c("Control-M", "12 hz-M")), map_signif_level = TRUE) +
             geom_signif(comparisons = list(c("Control-F", "12 hz-F")), map_signif_level = TRUE) + 
  geom_signif(comparisons = list(c("Control-M", "16 hz-M")), map_signif_level = TRUE) +
  geom_signif(comparisons = list(c("Control-F", "16 hz-F")), map_signif_level = TRUE) + 
  geom_signif(comparisons = list(c("Control-M", "18 hz-M")), map_signif_level = TRUE) +
  geom_signif(comparisons = list(c("Control-F", "18 hz-F")), map_signif_level = TRUE)

  


g1
g2 <- ggplot(df, aes(x = `Experimental group`, y = `39.5°С`, fill = Gender, group = Gender)) + 
             geom_bar(stat = "identity", alpha = 0.7) + 
             geom_errorbar(aes(ymin=`39.5°С`-sd2, ymax=`39.5°С`+sd2)) + 
             theme_bw() + 
             ggtitle("39.5°С Thermotolerance") + 
             theme(plot.title = element_text(hjust = 0.5)) +
             scale_y_continuous("% of alive", limits = c(-3,130)) +
             geom_signif(comparisons = list(c("12 hz-F", "control-M")), map_signif_level = TRUE)

g2
ggsave("thermotol-39.pdf", plot = g1, height = 10)
ggsave("thermotol-39-5.pdf", plot = g2, height = 10)


#prcomp
x <- df
x$Gender <- NULL
x <- x[,2:3]
p <- prcomp(x, scale. = TRUE, center = TRUE)
ggbiplot(p, groups = df$`Experimental group`)
