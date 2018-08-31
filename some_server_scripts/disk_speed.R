#!/usr/bin/Rscript

while(TRUE){
    v <- vector()
    cols <- as.numeric(system("tput cols", intern = TRUE))
    lines <- as.numeric(system("tput lines", intern = TRUE))
    for(f in 1:10){
      dat_1 <- system("df -BM | grep raid", intern = T)
      free_sp_1 <- as.numeric(sub("M", "",strsplit(dat_1," ")[[1]][13]))
      Sys.sleep(1)
      dat_2 <- system("df -BM | grep raid", intern = T)
      free_sp_2 <- as.numeric(sub("M", "",strsplit(dat_2," ")[[1]][13]))
      free <- system("df -h | grep raid", intern = T)
      all_free <- strsplit(free," ")
      all_free_num <- as.numeric(sub("G","", all_free[[1]][22]))
      speed <- free_sp_2 - free_sp_1
      v <- append(speed, v)
    }
    
    system("clear")
    txtplot::txtplot(v,
                     width = cols, 
                     height = lines, 
                     xlab = paste(as.character(all_free_num), "Gb free, average speed - ", round(mean(v)/6,2), "Mb/m"),
                     ylim = c(-200,200))
}







