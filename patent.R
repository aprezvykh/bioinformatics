library(tm)
library(ggplot2)
setwd("~/")
stop.words <- stopwords(kind = "en")
'%!in%' <- function(x,y)!('%in%'(x,y))
df <- read.csv("patent.csv", header = T,stringsAsFactors = F)


# Выявление наиболее часто встречающихся в выборке патентодателей ---------

frq.tab <- data.frame(table(t(df$Заявитель.и.)))
frq.tab <- frq.tab[order(frq.tab$Freq, decreasing = T),][1:10,]
png("patent/1_patent_freq.png", height = 250, width = 250, res = 600, units = "mm")
ggplot(data=frq.tab) + geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity") + 
  theme(text=element_text(family="Liberation Serif", face="bold", size=12,colour = "black")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",size = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.text=element_text(size=12, color = "black", face = "bold"),
        axis.title=element_text(size=12,face="bold", color = "black"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.line.y = element_line(color = "black",size = 0.3)) + 
  xlab("Заявитель") + 
  ylab("Количество патентов")
dev.off()


# выявление наиболее частых изобретателей ---------------------------------


frq.tab <- data.frame(table(t(df$Изобретатель.и.)))
frq.tab <- frq.tab[order(frq.tab$Freq, decreasing = T),][1:10,]
png("patent/1_patent_freq_inventor.png", height = 250, width = 250, res = 600, units = "mm")
ggplot(data=frq.tab) + geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity") + 
  theme(text=element_text(family="Liberation Serif", face="bold", size=12,colour = "black")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",size = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.3),
        axis.text=element_text(size=12, color = "black", face = "bold"),
        axis.title=element_text(size=12,face="bold", color = "black"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.line.y = element_line(color = "black",size = 0.3)) + 
  xlab("Изобретатель") + 
  ylab("Количество патентов")
dev.off()


# Text-mining названий патентов -------------------------------------------
word.freq <- unlist(lapply(df$Название, tolower))
words <- unlist(strsplit(word.freq, " "))
words <- words[words %!in% stop.words]

wf <- data.frame(table(t(words)))
wf <- wf[order(wf$Freq, decreasing = T),]
png("patent/3_word_cloud.png", height = 250, width = 250, res = 600, units = "mm")
wordcloud::wordcloud(words = wf$Var1, freq = wf$Freq,
                     max.words=500, random.order=FALSE, rot.per=0.35, 
                     colors=brewer.pal(8, "Dark2"))
dev.off()

