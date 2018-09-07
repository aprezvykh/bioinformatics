library("xlsx")
setwd("~/")
df <- read.xlsx("Sulfur Longevity.xlsx",sheetIndex = 10)
s <- melt(data = df,id.vars = "hours")
####  

males <- s[grepl("_M", s$variable),]
females <- s[grepl("_F", s$variable),]

####PARAQUAT
pdf("Sulfur-4.03-two_exp.pdf")
ggplot(data = s) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat, ALL") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat, males") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat, females") + 
  scale_y_continuous(name = "alive, %")


###CBS
cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat+CBS, males") + 
  scale_y_continuous(name = "alive, %")

cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat+CBS, females") + 
  scale_y_continuous(name = "alive, %")

#####CSE
cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat+CSE, males") + 
  scale_y_continuous(name = "alive, %")

cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Paraquat+CSE, females") + 
  scale_y_continuous(name = "alive, %")


###CNT
cbs <- df[,grep("G58", colnames(df))]
cbs$hours <- df$hours
cbs.females <- melt(cbs, id.vars = "hours")
ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp2 - Control lines only + Paraquat") + 
  scale_y_continuous(name = "alive, %")


####CONTROL

df <- read.xlsx("Sulfur Longevity.xlsx",sheetIndex = 8)
s <- melt(data = df,id.vars = "hours")
####  

males <- s[grepl("_M", s$variable),]
females <- s[grepl("_F", s$variable),]


ggplot(data = s) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, ALL") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, males") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, females") + 
  scale_y_continuous(name = "alive, %")


###CBS
cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, males") + 
  scale_y_continuous(name = "alive, %")

cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, females") + 
  scale_y_continuous(name = "alive, %")

#####CSE
cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, males") + 
  scale_y_continuous(name = "alive, %")

cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Exp1 - Paraquat, females") + 
  scale_y_continuous(name = "alive, %")


###CNT
cbs <- df[,grep("G58", colnames(df))]
cbs$hours <- df$hours
cbs.females <- melt(cbs, id.vars = "hours")
ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("Control lines only") + 
  scale_y_continuous(name = "alive, %")



####MEANS ALL
df.1 <- read.xlsx("Sulfur Longevity.xlsx",sheetIndex = 8)
df.2 <- read.xlsx("Sulfur Longevity.xlsx",sheetIndex = 10)

df <- data.frame(df.1$hours, df.1$mean, df.2$mean)
names(df) <- c("hours", "control", "paraquat")
s <-  melt(data = df,id.vars = "hours")

ggplot(data = s) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  geom_point(aes(x = hours, y = value, color = variable),size = 3) + 
  theme_bw() + 
  ggtitle("Exp1 and Exp2 means") + 
  scale_y_continuous(name = "alive, %")

dev.off()




####ALL_COMMON
####ALL_COMMON
####ALL_COMMON
####ALL_COMMON
####ALL_COMMON

df <- read.xlsx("Sulfur Longevity.xlsx",sheetIndex = 12)
s <- melt(data = df,id.vars = "hours")
####  

males <- s[grepl("_M", s$variable),]
females <- s[grepl("_F", s$variable),]

####PARAQUAT
pdf("Sulfur-4.03-merged.pdf")
ggplot(data = s) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("All merged, ALL") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("All merged, males") + 
  scale_y_continuous(name = "alive, %")

ggplot(data = females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("All merged, females") + 
  scale_y_continuous(name = "alive, %")


cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("CBS, males") + 
  scale_y_continuous(name = "alive, %")


cbs <- df[,grep("cbs", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("CBS, females") + 
  scale_y_continuous(name = "alive, %")

#####CSE
cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_M <- df$G58492_M
cbs.males <- melt(cbs, id.vars = "hours")
cbs.males <- cbs.males[grepl("_M", cbs.males$variable),]

ggplot(data = cbs.males) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("CSE, males") + 
  scale_y_continuous(name = "alive, %")

cbs <- df[,grep("CSE", colnames(df))]
cbs$hours <- df$hours
cbs$control_F <- df$G58492_F
cbs.females <- melt(cbs, id.vars = "hours")
cbs.females <- cbs.females[grepl("_F", cbs.females$variable),]

ggplot(data = cbs.females) + geom_line(aes(x = hours, y = value, color = variable), stat = "identity", size = 1) + 
  theme_bw() + 
  ggtitle("CSE, females") + 
  scale_y_continuous(name = "alive, %")


dev.off()






