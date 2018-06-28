library(xlsx)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
col.pan <- brewer.pal(n=12, name = "Paired")
df <- NULL
df <- read.xlsx("~/Sulfur Longevity.xlsx", sheetIndex = 4, header = TRUE)

df$yw_M <- NULL
df$yw_F <- NULL
colnames(df) <- sub("X", "", colnames(df))
s <- melt(data = df,id.vars = "hours")

df.2 <- read.xlsx("~/Sulfur Longevity.xlsx", sheetIndex = 5)
colnames(df.2) <- sub("X", "", colnames(df.2))
df.2 <- melt(data = df,id.vars = "hours")

males <- s[grepl("_M", s$variable),]
females <- s[grepl("_F", s$variable),]

males.2 <- s[grepl("_M", df.2$variable),]
females.2 <- s[grepl("_F", df.2$variable),]


g0 <- ggplot(data = s) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("ALL") + 
  scale_y_continuous(name = "% alive flies")

g3 <- ggplot(data = df.2) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("New exp") + 
  scale_y_continuous(name = "% alive flies")

g4 <- ggplot(data = males.2) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("New exp - males") + 
  scale_y_continuous(name = "% alive flies")
g5 <- ggplot(data = females.2) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("New exp - females") + 
  scale_y_continuous(name = "% alive flies")


g1 <- ggplot(data = males) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("MALES - SUMMARY") + 
  scale_y_continuous(name = "% alive flies")


g2 <- ggplot(data = females) + geom_line(aes(x = hours, y = 100-value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("FEMALES - SUMMARY") + 
  scale_y_continuous(name = "% alive flies")

pdf("Kapl_Meyer_lines_ALIVE.pdf")
g0
g3
g4
g5
g1
g2
dev.off()

pdf("Kapl_Meyer_barplots_ALIVE.pdf")
ggplot(data = s, aes(x = hours, y = 100-value)) + geom_bar(aes(group = variable, fill = variable), color = "black", stat = "identity", position = "dodge") + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,24,36,48,60,72,84,96,120,168)) + 
  scale_y_continuous(name = "% alive flies") + 
  ggtitle("ALL")

ggplot(data = males, aes(x = hours, y = 100-value)) + geom_bar(aes(group = variable, fill = variable), color = "black", stat = "identity", position = "dodge") + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,24,36,48,60,72,84,96,120,168),name = "% alive flies") + 
  ggtitle("MALES") + 
scale_y_continuous(name = "% alive flies")

ggplot(data = females, aes(x = hours, y = 100-value)) + geom_bar(aes(group = variable, fill = variable), color = "black", stat = "identity", position = "dodge") + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,24,36,48,60,72,84,96,120,168),name = "% alive flies") + 
  ggtitle("FEMALES") + 
scale_y_continuous(name = "% alive flies")


ggplot(data = females.2, aes(x = hours, y = 100-value)) + geom_bar(aes(group = variable, fill = variable), color = "black", stat = "identity", position = "dodge") + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,24,36,48,60,72,84,96,120,168),name = "% alive flies") + 
  ggtitle("FEMALES_new only") + 
  scale_y_continuous(name = "% alive flies")


ggplot(data = males.2, aes(x = hours, y = 100-value)) + geom_bar(aes(group = variable, fill = variable), color = "black", stat = "identity", position = "dodge") + 
  theme_bw() + 
  scale_x_continuous(breaks=c(0,24,36,48,60,72,84,96,120,168),name = "% alive flies") + 
  ggtitle("MALES_new only") + 
  scale_y_continuous(name = "% alive flies")


dev.off()
