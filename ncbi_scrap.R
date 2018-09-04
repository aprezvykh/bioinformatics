#!/usr/bin/Rscript
library(stringr)
setwd("~/ncbi_scrap/")
medline = function(file_name){
  lines <- readLines(file_name)
  medline_records <- list()
  key <- 0
  record <- 0
  for(line in lines){
    header <- sub(" {1,20}", "", substring(line, 1, 4))
    value <- sub("^.{6}", "", line)
    if(header == "" & value == ""){
      next
    }
    else if(header == "PMID"){
      record = record + 1
      medline_records[[record]] <- list()
      medline_records[[record]][header] <- value
    }
    else if(header == "" & value != ""){
      medline_records[[record]][key] <- paste(medline_records[[record]][key], value)
    }
    else{
      key <- header
      if(is.null(medline_records[[record]][key][[1]])){
        medline_records[[record]][key] <- value
      }
      else { 
        medline_records[[record]][key] <- paste(medline_records[[record]][key], value, sep=";")
      }
    }
  }
  return(medline_records)
}
data <- medline("medline.txt")
saveRDS(data, "medline.rds")

data <- readRDS("medline.rds")

tags <- as.character(lapply(data, function(x){x$OT}))
years <- as.character(lapply(data, function(x){x$DP}))
df <- data.frame(tags,years,stringsAsFactors = F)
tags_str <- strsplit(df$tags, ";")
years_str <- strsplit(df$years, " ") 
df$year_only <- as.character(lapply(years_str, function(x)x[1]))

#df[which(is.null(df$tags)),] <- "nodata"
all <- data.frame()
for(f in 1:length(n)){
    un <- unique(df$year_only)
    sub <- df[which(df$year_only == un[f]),]$tags
    sub <- lapply(sub,function(x){strsplit(x, ";")})
    sub <- as.character(unlist(sub, recursive = T))
    frq <- data.frame(table(t(sub)))
    frq <- frq[order(frq$Freq, decreasing = T),]
        if(nrow(frq<10)){
          print("nodata!")
        }
    df <- data.frame(frq[seq(1:10),], f)
    all <- rbind(df, all)
}
