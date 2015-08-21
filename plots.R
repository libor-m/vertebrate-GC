file <- "DanRer_all.txt"
file <- "MusMus_all.txt"
d <- read.delim(file, header=F)
colnames(d)<-c("seq", "ord", "A", "C", "G", "T")
d$sum <- d$A + d$C + d$G + d$T

library(ggplot2)
# this doesn't make much sense, only shows artifacts due to N letters in the sequences
ggplot(d, aes(A+T, C+G)) + geom_point() + facet_wrap(~seq)

ggplot(d, aes((C+G)/sum)) + geom_histogram() + facet_wrap(~seq, scales="free_y", ncol=5)

ggplot(d, aes((C+G)/sum)) + geom_histogram()

wrap <- facet_wrap(~seq, scales="free_x", ncol=5)
ggplot(d, aes(ord, (C+G)/sum)) + geom_line() + wrap

ggplot(d[d$seq == "11",], aes(ord, (C+G)/sum)) + geom_line()
