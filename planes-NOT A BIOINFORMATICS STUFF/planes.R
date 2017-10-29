library("ggmap")
library("ggplot2")
library("mapproj")
library("ggmap")
library(maptools)
mylat <- 55.11600477
mylon <- 36.60882175
map <- get_map(location = "Obninsk", zoom = 9)
directory <- "~/GitHub/data/plane-data/"
g <- grep('csv', list.files(directory),value=TRUE)
if (identical(g, character(0)) == TRUE){
  png(filename = "lcd.png")
  plot(x = 1, main = "NO DATA AVALLIBLE")
}

a <- read.csv(g)

names(a) <- c("MsgCount", "HexIdent","Date","Time","Lat","Long","Callsign","Altitude","Speed","Track","Vertical")
u <- unique(a$Callsign)
if (length(u) == 0){
  png(filename = "lcd.png")
  plot(x = 1, main = "NO DATA AVALLIBLE")
}

all <- data.frame()
df <- data.frame()
for (f in u){
  r <- grep(paste(f), a$Callsign)
  s <- data.frame(a[r,])
  s$Date <- NULL
  s$Time <- NULL
  s <- s[complete.cases(s),]
  df <- data.frame(mean(s$Lat), mean(s$Long), mean(s$Altitude),
                   mean(s$Speed), mean(s$Track), mean(s$Vertical))
  df$Callsign <- paste(f)
  names(df) <- c("Lat","Long","Altitude","Speed","Track","Vertical", "Callsign")
  all <- rbind(df, all)

}


png(filename="lcd.png", width = 640, height = 480)
ggmap(map) + geom_point(data = all, aes(x = Long, y = Lat, color = Altitude, size = Speed)) +
geom_text(data = all, aes(x = Long, y = Lat, label = Callsign), hjust = 1)
dev.off()


ggmap(map) + geom_point(data = a, aes(x = Long, y = Lat, color = Altitude, size = Speed))

dir.create("plane graphs")
setwd("~/plane graphs/")

for (f in u){
  d <- print(f)
  r <- grep(paste(f), a$Callsign)
  s <- data.frame(a[r,])
  

}


trackAngle <- function(xy) {
  angles <- abs(c(trackAzimuth(xy), 0) -
                  c(0, rev(trackAzimuth(xy[nrow(xy):1, ]))))
  angles <- ifelse(angles > 180, 360 - angles, angles)
  angles[is.na(angles)] <- 180
  angles[-c(1, length(angles))]
}

a$dist <- ((a$Long-mylon)^2 + (a$Lat - mylat)^2)*100

a <- a[order(a$Callsign),]
mat <- NULL
mat <- as.data.frame(a$Long, a$Lat)
h <- as.data.frame(trackAngle(mat))
