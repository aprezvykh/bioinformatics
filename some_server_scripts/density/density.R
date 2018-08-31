x <- read.table("~/R_scripts/density/dens.txt")
d <- density(x$V1)
plot(d)
