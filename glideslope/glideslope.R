library(ggplot2)
library(xlsx)
df <- read.xlsx("~/GitHub/glideslope/glide.xlsx", sheetIndex = 1)

ggplot(df) + geom_point(aes(x = l, y = alt), size = 2) + 
             theme_bw() + 
             geom_segment(aes(x = l, y = alt, xend = (l-speed/100), yend = alt),
             arrow=arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) + 
             geom_segment(aes(x = l, y = alt, xend = l, yend = alt-vert),
             arrow=arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
             geom_text(aes(label = n, x = l, y = alt), vjust = -1) + 
             geom_text(aes(label = speed, x = l, y = alt), hjust = 3, vjust = -1) +
             geom_text(aes(label = vert, x = l, y = alt), hjust = -1, vjust = 1) + 
             geom_segment(aes(x = 0, y = 0, xend = 10, yend = 535), color = "green") + 
             geom_segment(aes(x = 0, y = 0, xend = 10, yend = 465), color = "green") + 
             geom_segment(aes(x = 0, y = 0, xend = 10, yend = 600), color = "red") + 
             geom_segment(aes(x = 0, y = 0, xend = 10, yend = 400), color = "red") + 
             scale_x_continuous(name = "Radial, km") + 
             scale_y_continuous(name = "Altitude, m")

  
             

  