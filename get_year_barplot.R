####get_year_timeline
df <- read.delim("~/ncbi_scrap/memory/acc.txt", header = F, stringsAsFactors = F)
dfa <- strsplit(df$V1, " ") 
df <- unlist(dfa)
l <- df[grep("/", df, fixed = T)]
l <- as.Date(l)
l <- substring(l,1,4)
frq <- table(t(l))
barplot(frq)



