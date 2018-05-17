library(Rsamtools,quietly = TRUE)
setwd("~/genome_spades/")
bam <- scanBam("aligned.bam")

.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}
check_pos <- function(x){
  if (intToBits(x)[3] == 1){
    return(F)
  } else if (intToBits(x)[5] != 1){
    return(T)
  } else {
    return(F)
  }
}
check_neg <- function(x){
  if (intToBits(x)[5] == 1){
    return(T)
  } else {
    return(F)
  }
}


bam_field <- names(bam[[1]])
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field
dim(bam_df)


sequence_name = c("scaffold_12928")

test <- subset(bam_df, rname == sequence_name)
sequence_name_neg <- bam_df[bam_df$rname == sequence_name &
                      apply(as.data.frame(bam_df$flag), 1, check_neg),
                    'pos'
                    ]

sequence_name_pos <- bam_df[bam_df$rname == sequence_name &
                      apply(as.data.frame(bam_df$flag), 1, check_pos),
                    'pos'
                    ]
sequence_name_neg_density <- density(sequence_name_neg)
sequence_name_pos_density <- density(sequence_name_pos)

sequence_name_neg_density$y <- sequence_name_neg_density$y * -1

plot(sequence_name_pos_density,
     ylim = range(c(sequence_name_neg_density$y, sequence_name_pos_density$y)),
     main = paste("Coverage of ", sequence_name, sep = ""),
     xlab = paste(sequence_name, sep = ""),
     col = 'blue',
     type='h'
)
lines(sequence_name_neg_density, type='h', col = 'red')





