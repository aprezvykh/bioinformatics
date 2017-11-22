tr_up <- read.xlsx("~/trends.xlsx", sheetIndex = 1, header = TRUE)
tr_down <- read.xlsx("~/trends.xlsx", sheetIndex = 2)
names(tr_up) <- c("gene", "co1", "co2", "symbol", "name")
names(tr_down) <- c("gene", "co1", "co2", "symbol")
rownames(tr_up) <- tr_up$gene
rownames(tr_down) <- tr_down$gene
df <- read.xlsx("~/counts_ens/2_late_tg_vs_ctrl_tg/Results diffexpression.xlsx", sheetIndex = 3)
df$log2FoldChange <- df$log2FoldChange*(-1)
rownames(df) <- df$NA.

sub_down <- df[rownames(tr_down),]
sub_down <- sub_down[complete.cases(sub_down), ]

sub_up <- df[rownames(tr_up),]
sub_up <- sub_up[complete.cases(sub_up), ]
