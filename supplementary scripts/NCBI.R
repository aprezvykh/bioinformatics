install.packages("rentrez")
library("rentrez")
dis <- data.frame()
for (a in final$symbol){
  r_search <- entrez_search(db="medgen", term=paste(a))
  searchid <- r_search$id
  for (b in searchid){
    summary <- entrez_summary(db="medgen", id=paste(b))
    terms <- summary$title
  }
  df <- c(a, terms)
  dis <- rbind(df, dis)
  
}
dis <- data.frame()
for (a in final$symbol){
  r_search <- entrez_search(db="omim", term=paste(a))
  s <- r_search$count
  df <- data.frame(a, s)
  dis <- rbind(df, dis)
}

dis <-subset(dis, s>0)
entrez_dbs()
r_search <- entrez_search(db="omim", term="Cyp26b1")
r_search$ids
summary <- entrez_summary(db="omim", id=608428)
summary$title