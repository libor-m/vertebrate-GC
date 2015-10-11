#
# plot 'wiggle' track along the whole genome nicely
# using ggplot
#
# to accomodate for varying contig/scaffold lengths
# all the chromosomes are concatenated together and 
# further split into blocks of the same length
#

library(plyr)
library(dplyr)
library(gtools)
library(ggplot2)

sortchrom <- function(df) df %>% mutate(seq=seq %>% factor(levels=seq %>% unique %>% mixedsort))

# sort chromosomes by size
chrom_sort_size <- function(d)
  d %>%
    group_by(seq) %>%
    summarise(size=as.numeric(max(pos))) %>%
    arrange(desc(size)) %>%
    ungroup

# sort chromosomes by 'natural' order
# use factor levels for the ordering,
# but go back to character vector afterwards
# not to disrupt the following code
chrom_sort_alnum <- function(d)
  d %>%
    group_by(seq) %>%
    summarise(size=as.numeric(max(pos))) %>%
    ungroup %>%
    sortchrom %>%
    arrange(seq) %>%
    mutate(seq=as.character(seq))

# summarise shorter blocks to 100k blocks
plot_pseudo <- function (d, title, bases_per_pseudo=2.5e8, bases_per_block=1e5, chromsort=chrom_sort_size) {
  d %>%
    group_by(seq, pos=round_any(pos, bases_per_block, floor)) %>%
    summarise(A=sum(A), C=sum(C), G=sum(G), T=sum(T)) %>%
    mutate(sum = A + C + G + T, GC = (C + G) / (sum + 1)) %>%
    ungroup ->
    d100k
  
  # reorder the data decreasingly by the contig size
  d100k %>% 
    chromsort %>%
    filter(size >= 2*bases_per_block) %>%
    mutate(seq_end=cumsum(size),
           alternate=ifelse(row_number() %% 2, "A", "B")) ->
    seqs_sorted
  
  # create a map alternating along the order of the chromosomes
  altmap <- seqs_sorted$alternate
  names(altmap) <- seqs_sorted$seq
  
  # add faceting variables to scaffold endpoints
  seqs_sorted %>%
    mutate(pseudo_seq_num=(seq_end %/% bases_per_pseudo) + 1,
           pseudo_seq0=paste(title, pseudo_seq_num),
           pseudo_seq=factor(pseudo_seq0, levels=unique(pseudo_seq0)),
           pseudo_pos=seq_end %% bases_per_pseudo) ->
    sbs_pseudo
  
  # push the data through the final pipeline
  d100k %>%
    # get rid of singletons
    group_by(seq) %>%
    filter(n() > 1) %>%
    ungroup %>%
    
    # sort the data in the desired order
    arrange(factor(seq, levels=seqs_sorted$seq), pos) %>% 
    
    # wrap the concatenated positions with modulo
    mutate(catpos=seq_along(pos) * bases_per_block,
           pseudo_seq_num=(catpos %/% bases_per_pseudo) + 1,
           pseudo_seq0=paste(title, pseudo_seq_num),
           pseudo_seq=factor(pseudo_seq0, levels=unique(pseudo_seq0)),
           pseudo_pos=catpos %% bases_per_pseudo,
           alt=altmap[seq]) %>%
    
    # plot it
    ggplot(aes(pseudo_pos, GC)) + 
    geom_line(aes(group=seq, colour=alt)) +
    geom_point(y=0.7, size=4, colour="red", shape="|", data=sbs_pseudo) +
    # geom_text(aes(label=seq), y=0.7, hjust=1.5, colour="red", data=sbs_pseudo %>% filter(size > 1.5e7)) +
    # facet_wrap(~pseudo_seq, ncol=1) +
    facet_grid(pseudo_seq ~ .) +
    ylim(c(.3,.7)) +
    xlim(c(0, bases_per_pseudo)) +
    ggtitle(title) +
    scale_color_manual(values=c("A" = "#000000", "B" = "#777777"), guide=F)
  
}
