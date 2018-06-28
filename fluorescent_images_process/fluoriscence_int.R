install.packages("rgdal")
library(jpeg) 
library(ggplot2) 
library(raster)
library(rgdal)
library(reshape2)

dir_cnt <- "~/immunofluoriscence/AhR_Tj_control+IND/Контроль/"
dir_exp <- "~/immunofluoriscence/AhR_Tj_control+IND/Индинол/"


list_cnt <- grep("tif", list.files(dir_cnt), value = TRUE)
list_exp <- grep("tif", list.files(dir_exp), value = TRUE)

get_value <- function(img_src){
  r1 <- raster(img_src)
  m <- matrix(2 * pointDistance(expand.grid(-5:5, -5:5), c(0, 0), lonlat=FALSE),
              ncol=11, nrow=11)
  r2 <- focal(r1 >= 250, m, sum, na.rm=TRUE, pad=TRUE)
  r3 <- focal(r1 <= 3, m, sum, na.rm=TRUE, pad=TRUE)
  f_pix <- sum(r2[] >= 200 & r2[] <= 1000) * res(r2)[1]^2
  contour <- sum(r3[] == 0) * res(r3)[1]^2
  return((f_pix/contour)*100)  
}



setwd(dir_cnt)
cnt <- lapply(list_cnt, get_value)
cnt <- as.numeric(cnt)
cnt
setwd(dir_exp)
exp <- lapply(list_exp, get_value)
exp <- as.numeric(exp)

paste(round(mean(cnt),4), "±", round(sd(cnt),4), sep = "")
paste(round(mean(exp),4), "±", round(sd(exp),4), sep = "")
df <- data.frame(head(cnt,10), exp)
names(df) <- c("Control", "Indinol")

t.test(df$Control, df$Indinol)

plot(Control~Indinol, data = df)
abline(lm(Control~Indinol, data = df), col = "red")


df <- melt(df)
names(df) <- c("Group", "Intensity")
ggplot(data=df) + geom_boxplot(aes(x = Group, 
                                   y = Intensity, 
                                   color = Group, 
                                   fill = Group), alpha = 0.5) + 
                                   geom_jitter(aes(x = Group, y = Intensity)) + 
                                   theme_bw()
