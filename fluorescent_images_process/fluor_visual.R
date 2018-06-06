library(jpeg)
dir_cnt <- "~/Documents/immunofluorescence_images/AhR_Tj_control+IND/Контроль//"
dir_exp <- "~/Documents/immunofluorescence_images/AhR_Tj_control+IND/Индинол//"


list_cnt <- grep("tif", list.files(dir_cnt), value = TRUE)
jpg_cnt <- grep("jpg", list.files(dir_cnt), value = TRUE)

list_exp <- grep("tif", list.files(dir_exp), value = TRUE)
jpg_exp <- grep("jpg", list.files(dir_exp), value = TRUE)

setwd(dir_exp)
for (f in seq(1:11)){
    print(f)
    r1 <- raster(list_exp[f])
    j <- readJPEG(jpg_exp[f])
    m <- matrix(2 * pointDistance(expand.grid(-5:5, -5:5),
                                  c(0, 0), lonlat=FALSE),
                                  ncol=11, nrow=11)
    r2 <- focal(r1 >= 250, m, sum, na.rm=TRUE, pad=TRUE)
    r3 <- focal(r1 <= 5, m, sum, na.rm=TRUE, pad=TRUE)
    unf_pix <- sum(r1[] >= 0 & r1[] <= 255) * res(r1)[1]^2
    f_pix <- sum(r2[] >= 200 & r2[] <= 1000) * res(r2)[1]^2
    contour <- sum(r3[] == 0) * res(r3)[1]^2
    print(res(r1)[1]^2)
    print(res(r2)[1]^2)
    print(res(r3)[1]^2)
    png(paste(f, ".png", sep = ""))
    par(mfrow = c(2,2)) 
    plot(r1)
    plot(r2)
    plot(r3)
    plot(1:2, type='n')
    rasterImage(j,1,1,2,2)
    dev.off()
}

foca
