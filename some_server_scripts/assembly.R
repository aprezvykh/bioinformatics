library(dplyr)
source("http://bioconductor.org/biocLite.R")
biocLite("pachterlab/sleuth")
library("sleuth")
base_dir <- "~/transcriptomes/reads/kallisto/Kallisto_output/"
setwd(base_dir)
sample_id <- grep("dm", list.files(base_dir), value = T)
base_dir <- "~/transcriptomes/reads/kallisto/"


kal_dirs <- sapply(sample_id, function(id)file.path(base_dir, "Kallisto_output", id, "Quants"))
s2c <- read.table("samples_info.txt", header = TRUE, stringsAsFactors=FALSE)

s2c <- dplyr::select(s2c, sample = sample_name,
                     condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

so <- sleuth_prep(s2c, ~condition, num_cores = 32)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table <- sleuth_results(so,'reduced:full', test_type = 'lrt')
write.table(results_table,file="Out.txt",sep="\t")
sleuth_live(so)
