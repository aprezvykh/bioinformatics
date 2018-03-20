a <- read.xlsx("~/counts/AIKAR.dmel.full/results/CONTROL_LARVAE-5MM_LARVAE/Results edgeR.xlsx", sheetIndex = 3)
b <- read.xlsx("~/counts/AIKAR.dmel.full/results/5MM_LARVAE-10MM_LARVAE/Results edgeR.xlsx", sheetIndex = 3)

int <- intersect(a$NA., b$NA.)


a.com <- a[(a$NA. %in% int),]
b.com <- b[(b$NA. %in% int),]

df <- data.frame(a.com$NA., a.com$logFC, b.com$logFC, a.com$PValue, b.com$PValue, a.com$FDR, b.com$FDR, a.com$symbol, a.com$name, a.com$term)
write.xlsx(df, file = "AIKAR_dose_dependent.xlsx", sheetName = "all")
