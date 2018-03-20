late <- read.xlsx("~/counts/AIKAR.elim.2/results/control_late-aikar_late/Results edgeR.xlsx", sheetIndex = 3)
early <- read.xlsx("~/counts/AIKAR.elim.3.early/results/control_early-aikar_early/Results edgeR.xlsx", sheetIndex = 3)
control <- read.xlsx("~/counts/AIKAR.elim.control/results/control_early-control_late/Results edgeR.xlsx", sheetIndex = 3)


i <- intersect(control$NA., early$NA.)

outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

a <- early[(early$NA. %in% i),]
b <- control[(control$NA. %in% i),]
c <- data.frame(a$NA., a$logFC, b$logFC)

c$index <- c$a.logFC*c$b.logFC
c <- c[which(c$index < 0),]
elim.early <- c$a.NA.
early.flt <- outersect(early$NA., elim.early)
early.flt.final <- early[early.flt,]

setwd("~/counts/AIKAR/")
write.xlsx(early.flt.final, file = "filtered diffexpression ET_tt_LFC0,5.xlsx", sheetName = "control_early-aikar_early", append = TRUE)


i <- intersect(control$NA., late$NA.)


a <- late[(late$NA. %in% i),]
b <- control[(control$NA. %in% i),]
c <- data.frame(a$NA., a$logFC, b$logFC)

c$index <- c$a.logFC*c$b.logFC
c <- c[which(c$index < 0),]
elim.late <- c$a.NA.
late.flt <- outersect(late$NA., elim.late)
late.flt.final <- late[late.flt,]

setwd("~/counts/AIKAR/")
write.xlsx(late.flt.final, file = "filtered diffexpression ET_tt_LFC0,5.xlsx", sheetName = "control_late-aikar_late", append = TRUE)
