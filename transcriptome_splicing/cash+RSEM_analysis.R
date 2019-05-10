library(edgeR)
library(limma)
library(digest)
col.pan <- colorpanel(100, "blue", "white", "red")
dir <- "~/transcriptomes/reads/intron_retention/corr.all/"
sas <- read.csv("~/transcriptomes/reads/intron_retention/all.rsem/all.csv")
cash <- read.csv("~/transcriptomes/reads/intron_retention/all_cash.final.csv")

setwd(dir)
g <- grep("corr", list.files(dir), value = T)
r <- grep("res", list.files(dir), value = T)

corr <- data.frame()
res <- data.frame()

for(f in g){
  ms <- read.csv(f, stringsAsFactors = F)
  ms$st <- f
  corr <- rbind(ms, corr)
}

for(f in r){
  ms <- read.csv(f, stringsAsFactors = F)
  ms$st <- f
  res <- rbind(ms, res)
}


View(table(t(cash$ens)))
View(cash)
cash[grep("ENSMUSG00000026034", cash$ens),]

###табл по всем сплайс ивентам
###не корелляция
###описать все сплайс ивенты
###поделить барплот на ивенты
#последнюю картинку разбить на 3 штуки
#экспрессию цветом, экспгруппу - формой
###конкретные гены в динамике со сплайсингом
### кластеризация по К - средним
###написать мазину по поводу экзонов

#подробнее удержание интрона тг2-тг3 - где? на что может влияеть? в рамке или где?
#то же самое по микроэкзонам - РАМКА
###принципиально другие РНАсеки - печень, да все что угодно
### гены/сплайс изоформы - в мотонейронах или микроглии? ВАЖНО!!!
###
###общие гены в разных статтестах
###написать в vast-tools!!!
###сравнть с мазинскими данными
###написать мазину?
###читать статьи! 

df <- corr
ggplot(data=corr) + geom_point(aes(x = logCPM, y = dpsi, color = st)) + 
  geom_smooth(data=corr[corr$dpsi>0,],aes(x = logCPM,y = dpsi),se = F) +
  geom_smooth(data=corr[corr$dpsi<0,],aes(x = logCPM,y = dpsi),se = F) +
  theme_bw() + 
  ggtitle("Trends in dpsi-CPM")


df <- data.frame(corr$st, corr$type)
names(df) <- c("st", "type")
ggplot(data=df) + geom_bar(aes(x = st, y = type, color = type),stat = "identity")


df <- data.frame(table(t(corr[corr$type == "IR",]$st)))
ggplot(data=df) + geom_bar(aes(x = Var1, y = Freq), stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("IR events, corellate with expression")

df <- corr[corr$type == "IR",]
ggplot(data=df) + geom_point(aes(x = logFC, y = dpsi, color = st, size = abund)) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  geom_text(aes(x = logFC, y = dpsi, label = symbol), hjust = 1, vjust = 1) + 
  ggtitle("Retained introns - dPSI")
 
df <- corr[corr$type == "Cassette",]
ggplot(data=df) + geom_point(aes(x = logFC, y = dpsi, color = st, size = abund, shape = micro)) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  ggtitle("Cassette - dPSI")

df <- corr
ggplot(data=df) + geom_point(aes(x = logFC, y = dpsi, color = st, size = abund)) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  ggtitle("All events")

df <- corr[corr$micro == "micro",]
df <- df[grep("Cassette|MXE|Cassette_multi", df$type),]
ggplot(data=df) + geom_point(aes(x = logFC, y = dpsi, color = st, size = abund)) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  ggtitle("Microexons only")

df <- corr[corr$micro == "normal",]
df <- df[grep("Cassette|MXE|Cassette_multi", df$type),]
ggplot(data=df) + geom_point(aes(x = logFC, y = dpsi, color = st, size = abund)) + 
  theme_bw() + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  ggtitle("Normal exons only")


###

df <- data.frame(table(t(cash$st)))
ggplot(data=df) + geom_bar(aes(x = Var1, y = Freq, fill = Var1), stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("All alternative spliced events")

df <- cash
ggplot(data=df[df$SplicingType == "IR",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + 
  theme_bw() + 
  ggtitle("Intron retention events")

ggplot(data=df[df$SplicingType == "Cassette",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + 
  theme_bw() + 
  ggtitle("Cassete exons")
  
ggplot(data=df) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + 
  theme_bw() + 
  ggtitle("All events")

ggplot(data=df[df$micro == "micro" & df$SplicingType == "Cassette",]) + geom_boxplot(aes(x = st, y = delta_PSI, fill = st)) + 
  theme_bw() + 
  ggtitle("Cassette microexons")

##moto-glia-ir-supressed

ir.genes <- corr[corr$dpsi>0 & corr$logFC<0 & corr$type == "IR",]
moto <- read.csv("~/motoneuron genes.csv")
glia <- read.csv("~/glial genes.csv")

im <- intersect(ir.genes$ensgene, moto$compare.glia.X)
ig <- intersect(ir.genes$ensgene, glia$compare.glia.X)

moto.ir.supressed <- ir.genes[ir.genes$ensgene %in% im,]
glia.ir.supressed <- ir.genes[ir.genes$ensgene %in% ig,]


###microexonic

me.genes <- corr[grep("Cassette|MXE", corr$type),]

im <- intersect(me.genes$ensgene, moto$compare.glia.X)
ig <- intersect(me.genes$ensgene, glia$compare.glia.X)

moto.me <- me.genes[me.genes$ensgene %in% im,]
glia.me <- me.genes[me.genes$ensgene %in% ig,]

####moto and glia

ima <- intersect(corr$ensgene, moto$compare.glia.X)
iga <- intersect(corr$ensgene, glia$compare.glia.X)

moto.all <- corr[corr$ensgene %in% ima,]
glia.all <- corr[corr$ensgene %in% iga,]

sas.src <- sas
rownames(sas.src) <- sas.src$X
sas.src$X <- NULL
sas.src$SRR646651.rsem.isoforms.results <- NULL

model <- c("fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", 
           "fus", "fus", "fus", "fus", "fus",
           "tdp", "tdp", "tdp", "tdp", 
           "tdp", "tdp", "tdp", "tdp", 
           "sod-moto", "sod-moto",
           "sod-moto", "sod-moto",
           "sod-glia", "sod-glia", 
           "sod-glia", "sod-glia", "sod-glia", 
           "sod-glia", "sod-glia", "sod-glia", "sod-glia", 
           "sod-glia", "sod-glia")

age <- c(60, 60, 60, 60, 60, 
         120, 120, 120, 120, 120, 
         60, 60, 60, 60, 60, 
         90, 90, 90, 90, 
         120, 120, 120, 120, 120,
         38, 38, 38, 38, 
         38, 38, 38, 38, 
         90, 90, 
         90, 90, 
         65, 65, 
         100, 100, 100, 
         150, 150, 150, 150, 
         65, 65)

transgenity <- c("wt", "wt", "wt", "wt", "wt", 
                 "wt", "wt", "wt", "wt", "wt", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "wt", "wt", "wt", "wt", 
                 "wt", "wt", 
                 "tg", "tg", 
                 "tg", "tg", 
                 "tg", "tg", "tg", 
                 "tg", "tg", "tg", "tg", 
                 "wt", "wt") 

sample.id <- c(1,2,3,4,5,
               6,7,8,9,10,
               11,12,13,14,15,
               16,17,18,19,
               20,21,22,23,25,
               1,2,3,4,
               5,6,7,8,
               1,2,
               3,4,
               1,2,
               3,4,5,
               6,7,8,9,
               10,11)

metadata <- data.frame(model = model,
                       age = age,
                       transgenity = transgenity,
                       index = sample.id)

metadata$final <- paste(model, "_", age, "_", transgenity, "_", sample.id, sep = "")
names(sas) <- metadata$final
###linear model











