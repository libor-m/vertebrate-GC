library(ggplot2)
library(plyr) # only because of round_any()
library(dplyr)

#### GC profile along chromosomes ####

dfile <- "Branchio-100k.tsv"
dfile <- "profiles/Petromyzon-10k.tsv"
dfile <- "profiles/Coelacanth-10k.tsv"
d <- read.delim(dfile, col.names=c("seq", "ord", "A", "C", "G", "T")) %>%
  mutate(sum = A + C + G + T, GC = (C + G) / sum)

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
  ggplot(aes(ord, GC)) +
  geom_line() +
  facet_wrap(~seq) +
  ylim(c(0.3, 0.7))

# profiles, summarize in 100k windows
d %>% 
  filter(seq  %in% top20seqs) %>% 
  group_by(seq, ord=round_any(ord, 100000, floor)) %>%
  summarise(A=sum(A), C=sum(C), G=sum(G), T=sum(T)) %>%
  mutate(sum = A + C + G + T, GC = (C + G) / (sum+1)) %>% 
  ggplot(aes(ord, GC)) +
  geom_line() +
  facet_wrap(~seq) +
  ylim(c(0.3, 0.7))

#### GC profile for genomes with many short scaffolds ####

# for species with only small scaffolds
# join all the data into pseudo blocks of given length (250 Mb currently)

bases_per_block <- 1e5
blocks_per_pseudo <- 2.5e8 / bases_per_block

# reorder the data decreasingly by the contig size
d %>% 
  group_by(seq) %>%
  summarise(size=as.numeric(max(ord + 1))) %>%
  arrange(desc(size)) %>%
  mutate(seq_end=cumsum(size)) ->
  seqs_by_size

d %>%
  # summarise shorter blocks to 100k blocks
  group_by(seq, ord=round_any(ord, bases_per_block, floor)) %>%
  summarise(A=sum(A), C=sum(C), G=sum(G), T=sum(T)) %>%
  mutate(sum = A + C + G + T, GC = (C + G) / (sum+1)) %>%
  ungroup() %>%
  # start with the longest sequences
  arrange(factor(seq, levels=seqs_by_size$seq), ord) %>% 
  # create pseudo sequences, each containing required number of bases
  mutate(pseudo_seq_num=rep(1:30, length.out=n(), each=blocks_per_pseudo), 
         pseudo_seq0=paste("pseudo", pseudo_seq_num),
         pseudo_seq=factor(pseudo_seq0, levels=unique(pseudo_seq0))) %>% 
  group_by(pseudo_seq) %>%
  mutate(pseudo_ord=(seq_along(pseudo_seq) - 1) * bases_per_block) %>%
  # plot it
  ggplot(aes(pseudo_ord, GC)) + 
  geom_line() +
  facet_wrap(~pseudo_seq, ncol=1) +
  ylim(c(.3,.7))

#### GC summary for more species ####

d %>% 
  ggplot(aes(GC)) + 
  geom_histogram()

# read profile data, 
# add a filename to each row
# calculate stats
read.name <- function(fn) {
  read.delim(paste0("profiles/", fn), 
             col.names=c("seq", "ord", "A", "C", "G", "T"),
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
