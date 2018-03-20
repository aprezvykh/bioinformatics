a <- read.xlsx("~/counts/AIKAR/results/control_early-aikar_early/Results edgeR.xlsx", sheetIndex = 3)
b <- read.xlsx("~/counts/AIKAR/results/control_late-aikar_late/Results edgeR.xlsx", sheetIndex = 3)
f <- read.xlsx("~/counts/AIKAR/results/control_early-control_late/Results edgeR.xlsx", sheetIndex = 3)



outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}



i <- intersect(a$NA., b$NA.)

c1 <- outersect(a$NA., i)
c2 <- outersect(b$NA., i)

c3 <- intersect(f$NA., c1)
c4 <- intersect(f$NA., c2)

c5 <- outersect(c1, c3)
c6 <- outersect(c2, c4)


d <- a[(a$NA. %in% c5),]
e <- b[(b$NA. %in% c6),]



setwd("~/counts/AIKAR/results/")
write.xlsx(d, file = "cel_filtered_results_diffexpression.xlsx", sheetName = "early", append = TRUE)
write.xlsx(e, file = "cel_filtered_results_diffexpression.xlsx", sheetName = "late", append = TRUE)
getwd()
