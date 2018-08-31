#!/usr/bin/Rscript
setwd("~/R_scripts/")
loadavg.parse <- function(x){
  s <- x
  ss <- strsplit(s, ",")
  s <- as.numeric(paste(ss[[1]][1], ss[[1]][2], sep = "."))
  return(s)
}
while(TRUE){
    data <- data.frame()
    for(f in seq(1:600)){
      print(f)
      dat <- system("df -BM | grep raid", intern = T)
      free_sp <- as.numeric(sub("M", "",strsplit(dat," ")[[1]][13]))
      up <- system("uptime", intern = T)
      ram <- system("free -m", intern = T)
      df <- data.frame(strsplit(up, " ")[[1]][2],
                       free_sp,
                       loadavg.parse(strsplit(up, " ")[[1]][13]),
                       loadavg.parse(strsplit(up, " ")[[1]][14]),
                       loadavg.parse(strsplit(up, " ")[[1]][15]),
                       strsplit(ram, " ")[[2]][13],
                       strsplit(ram, " ")[[3]][16])
      data <- rbind(df, data)
      Sys.sleep(1)
    }
    names(data) <- c("time", "disk", "up1", "up2", "up3", "ram", "swap")
    name <- paste("alex_", tail(data$time,1), "_log.csv", sep = "")
    name <- sub(":","-",name)
    write.csv(data, paste("logs/",name, sep = ""))
    print(paste(name, " saved!", sep = ""))
}
