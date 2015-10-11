library(ggplot2)
library(plyr) # only because of round_any()
library(dplyr)

#### load the GC profile data ####

sources <- c(
  "Ciona-100k.tsv"="Sea squirt (Tunicate)",
  "Branchio-100k.tsv"="Lancelet",
  "Coelacanth-10k.tsv"="Coelacanth",
  "DanRer_all.txt"="Zebrafish",
  "Petromyzon-10k.tsv"="Lamprey",
  "EShark-1k.tsv"="Elephant shark",
  "LepOcu_all.txt"="Spotted gar",
  "MusMus_all.txt"="House mouse",
  "Human-100k.tsv"="Human"
)

dfile <- names(sources[3])
d <- read.delim(paste0('profiles/', dfile), col.names=c("seq", "pos", "A", "C", "G", "T")) %>%
  mutate(sum = A + C + G + T, GC = (C + G) / sum)

#### GC profile along chromosomes/scaffolds ####

d %>% 
  group_by(seq) %>%
  summarise(count=n()) %>%
  arrange(desc(count)) %>%
  head(n=20) %>%
  .$seq ->
  top20seqs

# profiles along top 20 longest scaffolds
d %>% 
  filter(seq  %in% top20seqs) %>% 
  ggplot(aes(pos, GC)) +
  geom_line() +
  facet_wrap(~seq) +
  ylim(c(0.3, 0.7))

# profiles, summarize in 100k windows
d %>% 
  filter(seq  %in% top20seqs) %>% 
  group_by(seq, pos=round_any(pos, 100000, floor)) %>%
  summarise(A=sum(A), C=sum(C), G=sum(G), T=sum(T)) %>%
  mutate(sum = A + C + G + T, GC = (C + G) / (sum+1)) %>% 
  ggplot(aes(pos, GC)) +
  geom_line() +
  facet_wrap(~seq) +
  ylim(c(0.3, 0.7))

#### GC profile for genomes with many short scaffolds ####

# for species with only small scaffolds
# join all the data into pseudo blocks of given length (250 Mb currently)
source('chromoplot.R')

load.plot.save <- function(fn, bases_per_pseudo=2.5e8) {
  d <- read.delim(paste0('profiles/', fn), 
                  col.names=c("seq", "pos", "A", "C", "G", "T")) %>%
    mutate(sum = A + C + G + T, GC = (C + G) / sum)
  
  d %>% plot_pseudo(sources[fn], bases_per_pseudo)
  
  d %>% 
    group_by(seq) %>%
    summarise(size=as.numeric(max(pos + bases_per_block))) %>%
    {(sum(.$size) %/% bases_per_pseudo) + 1} ->
    npanels

  # when there is too few panels, do not make them too high
  # do not go over A4 height on the other side
  h <- min(280, npanels * 40 + 50)
  
  ggsave(paste0('results/mondsee/pseudo-', sources[fn], '.pdf'), width=190, height=h, units = "mm")
}

# run for all input files
sapply(names(sources), load.plot.save)

#### test rescale function ####
block_rescale <- function(d, bases_per_block) {
  d %>%
    group_by(seq, pos=round_any(pos, bases_per_block, floor)) %>%
    summarise(A=sum(A), C=sum(C), G=sum(G), T=sum(T)) %>%
    mutate(sum = A + C + G + T, GC = (C + G) / (sum+1)) %>%
    ungroup()
}

ggplot(data.frame(), aes(pos, GC)) +
  geom_line(data=d %>% filter(seq  %in% top20seqs), colour="gray") +
  geom_line(data=d %>% filter(seq  %in% top20seqs) %>% block_rescale(1e5), colour="red") +
  facet_wrap(~seq) +
  ylim(c(0.3, 0.7))
ggsave('results/mondsee/block-scaling-test.pdf', width=280, height=190, units="mm")

#### GC summary for more species ####

d %>% 
  ggplot(aes(GC)) + 
  geom_histogram()

# read profile data, 
# add a filename to each row
# calculate stats
read.name <- function(fn) {
  read.delim(paste0("profiles/", fn), 
             col.names=c("seq", "pos", "A", "C", "G", "T"),
             colClasses=c("seq"="factor")) %>%
    mutate(sum = A + C + G + T, GC = (C + G) / sum, source=fn)
}

sources <- c(
  "Ciona-100k.tsv"="Sea squirt (Tunicate)",
  "Branchio-100k.tsv"="Lancelet",
  "Coelacanth-10k.tsv"="Coelacanth",
  "DanRer_all.txt"="Zebrafish",
  "Petromyzon-10k.tsv"="Lamprey",
  "EShark-1k.tsv"="Elephant shark",
  "LepOcu_all.txt"="Spotted gar",
  "MusMus_all.txt"="House mouse",
  "Human-100k.tsv"="Human"
)

# load all tables, each of them has the source in 
# source column
tabs <- lapply(names(sources), read.name)
da <- bind_rows(tabs)

# do summaries for all genomes

da %>%
  mutate(source=factor(source, levels=names(sources), labels=sources)) %>%
  ggplot(aes(GC)) + 
  geom_density(colour=NA, fill="#444444") +
  facet_wrap(~source) + 
  xlim(c(0.3, 0.6)) +
  xlab("GC fraction") +
  ylab("block density") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
ggsave('results/mondsee/gc-summary.pdf', width=7.4, height=5, units="in")

# spare parts
# normalized histogram
#geom_histogram(aes(y=..count../sum(..count..))) +
