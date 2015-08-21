#####
# interactive session file to plot GC profiles
#####
library(dplyr)   # data wrangling
library(gtools)  # mixedsort for chromosome order
library(ggplot2) # plotting

# choose one input file
file <- "DanRer_all.txt"
file <- "MusMus_all.txt"
file <- "LepOcu1_10k.txt"

# load the file and rename columns
d <- read.delim(file, header=F)
colnames(d)<-c("seq", "ord", "A", "C", "G", "T")
d$sum <- d$A + d$C + d$G + d$T

# sort the chromosomes in 'chromosome' order
d$seq <- factor(d$seq, levels=mixedsort(levels(d$seq)))

# get rid of too small chromosomes
# filter for gar
d <- filter(d, seq != "JH591410.1", seq != "MT", seq != "LG29")
# danio
d <- filter(d, seq != "MT", seq != "Zv9_scaffold3530")
# mouse
d <- filter(d, seq != "MT")

# order the data 
d <- arrange(d, seq, ord)

# 
# multi resolution plot
wrap <- facet_wrap(~seq, scales="free_x", ncol=5) 
ggplot(d, aes(ord, (C+G)/sum)) + geom_line() + wrap + ylim(0.25, 0.75)
ggsave("gar_10k.pdf", width=20, height=14, units="in")
ggsave("gar_10k.png", width=14, height=10, dpi=100, units="in")

# double plot - once for print, once for viewing
# A4 paper is 8.3 x 11.7 inches
myplot <- function (name) {
    ggsave(paste(name, "pdf", sep="."), width=20, height=14, units="in")
    ggsave(paste(name, "png", sep="."), width=14, height=10, dpi=100, units="in")
}

####
# this block sums 10k input data and creates another 
# resolution
####
#number of consecutive rows to merge
n_merge <- 100
counts <- d %.% group_by(seq) %.% summarize(rows=n())
# grouping - create sequence restarting on each chromosome
# and being repeated for n_merge times - so group_by is applied on desired rows
t <- sapply(counts$rows, function(n_max) rep(1:(n_max/n_merge + 1), each=n_merge)[1:n_max])
d$grp <- do.call(c, t)
dm <- d %.% group_by(seq, grp) %.% 
    summarize(ord=min(ord), A=sum(A), C=sum(C), G=sum(G), T=sum(T), sum=sum(sum)) %.%
    arrange(seq, ord)

ggplot(dm, aes(ord, (C+G)/sum)) + geom_line() + wrap + ylim(0.25, 0.75)
myplot(paste("gar_", n_merge*10, "k", sep=""))


###
# spare parts - other types of plots
###
# this doesn't make much sense, only shows artifacts due to N letters in the sequences
# scatter plots at-gc for each chromosome
ggplot(d, aes(A+T, C+G)) + geom_point() + facet_wrap(~seq)

# histogram of gc content for each chromosome
ggplot(d, aes((C+G)/sum)) + geom_histogram() + facet_wrap(~seq, scales="free_y", ncol=5)

# gc content over whole genome
ggplot(d, aes((C+G)/sum)) + geom_histogram()

# single chromosome
ggplot(d[d$seq == "11",], aes(ord, (C+G)/sum)) + geom_line()
