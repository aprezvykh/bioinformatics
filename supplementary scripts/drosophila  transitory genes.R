library(xlsx)
kn <- read.xlsx("~/counts/dr_multimap/results/K-N/Results edgeR.xlsx", sheetIndex = 3)
nf <- read.xlsx("~/counts/dr_multimap/results/N-F/Results edgeR.xlsx", sheetIndex = 3)
fm <- read.xlsx("~/counts/dr_multimap/results/F-M/Results edgeR.xlsx", sheetIndex = 3)

i <- intersect(kn$NA., nf$NA.)
ii <- intersect(i, fm$NA.)

kn_com <- kn[(kn$NA. %in% ii),]
nf_com <- nf[(nf$NA. %in% ii),]
fm_com <- fm[(fm$NA. %in% ii),]

a <- data.frame(kn_com$NA., kn_com$logFC, nf_com$logFC, fm_com$logFC, kn_com$symbol, kn_com$name)


names(a) <- c("flybase_id", "N_vs_K", "F_vs_N", "F24_vs_F", "symbol", "name")
rownames(a) <- a$symbol
a$symbol <- NULL
a$name <- NULL
a$flybase_id <- NULL
a <- as.data.frame(a)
a <- t(a)
ggplot() + geom_bar(data = a, aes(x = rownames(a), y = Hsp70Aa, stat = "identity", color = Hsp70Aa))

dev.off()
