### mean plot
setwd("~/adm/temp.heatmap/")
g <- grep("txt", list.files(), value = TRUE)
df <- data.frame()
for (f in g){
  x <- read.table(paste(f))
  dat <- x$V3
  dat <- gsub("°C", "", dat)
  dat <- gsub("+", "", dat)
  m <- mean(as.numeric(dat))
  df <- rbind(m, df)
  }
names(df) <- c("sas")
png("Temperature plot.png")
plot(df$sas, type = "l", xlab = "Time", ylab = "Temperature", ylim=c(20,80))
dev.off()

###multiple plot
df <- data.frame()
for (f in g){
  x <- read.table(paste(f))
  dat <- x$V3
  dat <- gsub("°C", "", dat)
  dat <- gsub("+", "", dat)
  dat <- as.numeric(dat)
  df <- rbind(dat, df)
  }

png("Temperature plot_multiple.png")
matplot(df, type = c("b"),pch=1, ylim = c(20,70))
dev.off()

