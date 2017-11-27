tg_12_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/glia/tg12-glia.csv")
tg_12_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/moto/tg12-moto.csv")
tg_12_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg2/other/tg12-other.csv")

tg_13_glia <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/glia/tg13-glia.csv")
tg_13_moto <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/moto/tg13-moto.csv")
tg_13_others <- read.csv("~/counts/ALS Mice/filtering fc=1 + FDR/Tg1-Tg3/other/tg13-other.csv")


lip <- read.xlsx("~/counts/ALS Mice/Lipids.xlsx", sheetIndex = 1)
getwd()
df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_12_glia$symbol, ignore.case = TRUE)
  df <- tg_12_glia[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg12-glia", append = TRUE)
}

df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_12_moto$symbol, ignore.case = TRUE)
  df <- tg_12_moto[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg12-moto", append = TRUE)
}

df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_12_others$symbol, ignore.case = TRUE)
  df <- tg_12_others[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg12-others", append = TRUE)
}

################



df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_13_glia$symbol, ignore.case = TRUE)
  df <- tg_13_glia[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg13-glia", append = TRUE)
}

df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_13_moto$symbol, ignore.case = TRUE)
  df <- tg_13_moto[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg13-moto", append = TRUE)
}

df <- data.frame()
lipid <- data.frame()

for(f in lip$found){
  a <- grep(paste(f), tg_13_others$symbol, ignore.case = TRUE)
  df <- tg_13_others[a,]
  lipid <- rbind(df, lipid)
  
}

if(nrow(lipid) > 0){
  write.xlsx(lipid, file = "Association.xlsx", sheetName = "Tg13-others", append = TRUE)
}
