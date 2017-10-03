experimental <- as.data.frame(read.csv("~/motoneurons compare/exp_cpm.csv"))
motoneurons <- as.data.frame(read.csv("~/motoneurons compare/exp_moto.csv"))

experimental$avg <- (experimental$mouse_ntg_1_176.counts +
                    experimental$mouse_ntg_1_177.counts +
                    experimental$mouse_ntg_1_178.counts +
                    experimental$mouse_ntg_1_213.counts +
                    experimental$mouse_ntg_1_225.counts)/5
motoneurons$avg <- (motoneurons$moto_las_1.counts +
                    motoneurons$moto_las_2.counts)/2

experimental$avg2 <- motoneurons$avg
experimental$mouse_ntg_1_176.counts <- NULL
experimental$mouse_ntg_1_178.counts <- NULL
experimental$mouse_ntg_1_225.counts <- NULL
experimental$mouse_ntg_1_177.counts <- NULL
experimental$mouse_ntg_1_213.counts <- NULL
rownames(experimental) <- experimental$X
experimental$X <- NULL
names(experimental) <- c("experimental", "moto")
