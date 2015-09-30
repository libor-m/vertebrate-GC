# version 4 utilizes the (ensembl) soft-masked genomes
# distinguishing between 'real' dna and TEs

library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd('c:/work/gar/')

read.name <- function(fn, datadir) {
  read.delim(paste0(datadir, "/", fn), 
             col.names=c("seq", "pos", "A", "C", "G", "T", "a", "c", "g", "t"),
             colClasses=c("seq"="factor")) %>%
    mutate(source=fn, 
           real_sum = A + C + G + T, 
           real_GC = (C + G) / real_sum, 
           mask_sum = a + c + g + t,
           mask_GC = (c + g) / mask_sum,
           total_sum = A + C + G + T + a + c + g + t, 
           total_GC = (C + c + G + g) / total_sum
           )
}

sources <- c(
  "Ciona_intestinalis.KH.dna_sm.toplevel.fa.gz.tsv" = "Sea squirt (Tunicate)",
  "Latimeria_chalumnae.LatCha1.dna_sm.toplevel.fa.gz.tsv" = "Latimeria",
  "Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz.tsv" = "Zebrafish",
  "Petromyzon_marinus.Pmarinus_7.0.dna_sm.toplevel.fa.gz.tsv" = "Lamprey",
  "Takifugu_rubripes.FUGU4.dna_sm.toplevel.fa.gz.tsv" = "Fugu",
  "Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz.tsv" = "Tetraodon",
  "Gasterosteus_aculeatus.BROADS1.dna_sm.toplevel.fa.gz.tsv" = "Stickelback",
  "Gadus_morhua.gadMor1.dna_sm.toplevel.fa.gz.tsv" = "Cod",
  "Lepisosteus_oculatus.LepOcu1.dna_sm.toplevel.fa.gz.tsv" = "Spotted gar",
  "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz.tsv" = "Mouse",
  "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.tsv" = "Human"
)

tabs <- lapply(names(sources), read.name, "profiles-v4")
da <- bind_rows(tabs)

# check how many 'short' blocks we have
da %>%
  mutate(source=factor(source, levels=names(sources), labels=sources)) %>%
  ggplot(aes(total_sum)) + 
  geom_density(colour=NA, fill="#444444") +
  facet_wrap(~source, ncol=3, scales="free")
  
# plot of total GC
da %>%
  mutate(source=factor(source, levels=names(sources), labels=sources)) %>%
  filter(total_sum > 5000) %>%
  ggplot(aes(total_GC)) + 
  geom_density(colour=NA, fill="#444444") +
  facet_wrap(~source, ncol=3) + 
  xlim(c(0.3, 0.6)) +
  xlab("total GC fraction") +
  ylab("block density") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
ggsave('results/profiles-v4/gc-summary.pdf', width=7.4, height=5, units="in")

# plot total gc on coarser scale
da %>%
  mutate(source=factor(source, levels=names(sources), labels=sources)) %>%
  group_by(source, seq, pos=round_any(pos, 5e6, floor)) %>%
  summarise(total_sum=sum(A + a + C + c + G + g + T + t),
            total_GC=sum(C + c + G + g) / total_sum) %>%
  filter(total_sum > 2500000) %>%
  ggplot(aes(total_GC)) + 
  geom_density(colour=NA, fill="#444444") +
  facet_wrap(~source, ncol=3) + 
  xlim(c(0.3, 0.6)) +
  xlab("total GC fraction") +
  ylab("block density") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())

ggsave('results/profiles-v4/gc-summary-5M.pdf', width=7.4, height=5, units="in")
# the result is, that scaling up the window 
# does not change the profiles much - that could mean that the blocks are larger than 1M..

# separate profiles for masked and normal areas
# (i. e. isolate variable 'masked' from column names)
# and calculate few more handy columns
da %>% 
  select(-contains("GC"), -contains("sum")) %>%
  gather(base, count, A:t) %>%
  mutate(compartment=ifelse(base %in% c('A', 'C', 'G', 'T'), "normal", "masked")) %>%
  mutate(base=toupper(base)) %>%
  spread(base, count) %>%
  mutate(sum=A + C + G + T,
         GC=(C + G) / sum,
         source=factor(source, levels=names(sources), labels=sources)) ->
  dm

# check the chunk sizes
dm %>%
  ggplot(aes(sum, fill=compartment)) +
  geom_histogram(position="identity", alpha=0.7) +
  facet_wrap(~source, scale="free_y") +
  ggtitle("length of blocks")
ggsave('results/profiles-v4/block-len.pdf', width=7.4, height=5, units="in")

# number of masked/normal regions
dm %>%
  filter(sum > 1000) %>%
  ggplot(aes(source, fill=compartment)) +
  geom_bar(position="dodge") +
  coord_flip() +
  ggtitle("Count of 10k blocks")
ggsave('results/profiles-v4/block-count.pdf', width=6, height=7.4, units="in")

# gc histograms
dm %>%
  filter(sum > 1000) %>%
  ggplot(aes(GC, fill=compartment)) +
  geom_density(colour=NA, alpha=0.7) +
  facet_wrap(~source, ncol=3) +
  xlim(c(0.3, 0.6)) +
  xlab("GC fraction") +
  ylab("block density") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
ggsave('results/profiles-v4/gc-masked.pdf', width=7.4, height=5, units="in")

# assembly snakes
# ( too many short contigs influence the bins )
dm %>% 
  group_by(source, seq) %>%
  summarise(len=max(pos)) %>%
  group_by(source) %>%
  arrange(desc(len)) %>%
  mutate(len=as.numeric(len),
         cumlen=cumsum(len), 
         cumrel=cumlen/max(cumlen), 
         nx=cut(100*cumrel, 0:100, labels=1:100)) %>%
  group_by(source, nx) %>% 
  summarize(mlen=min(len)) %>%
  mutate(nx=as.numeric(nx)) ->
  dnx

dnx %>%
  ggplot(aes(nx, mlen, colour=source)) + 
  geom_line(size=2) +
  geom_vline(xintercept=50, colour="gray") + 
  geom_vline(xintercept=90, colour="gray") +
  xlab("N(x)") + 
  ylab("base pairs") + 
  ggtitle("N(x) values for different assemblies")
ggsave('results/profiles-v4/assembly-snakes.pdf', width=7.4, height=5, units="in")

# profiles along chromosomes
# run plot_pseudo from profiles-v2.R
#
plot.save <- function(species, din, prefix, bases_per_pseudo=2.5e8, bases_per_block=1e5) {
  din %>%
    filter(source == species) ->
    d
  
  d %>% plot_pseudo(species, bases_per_pseudo)
  
  d %>% 
    group_by(seq) %>%
    summarise(size=as.numeric(max(pos + bases_per_block))) %>%
    {(sum(.$size) %/% bases_per_pseudo) + 1} ->
    npanels

  # when there is too few panels, do not make them too high
  # do not go over A4 height on the other side
  h <- min(280, npanels * 40 + 50)
  
  ggsave(paste0(prefix, species, '.pdf'), width=190, height=h, units = "mm")
}

# test plots
plot.save("Fugu", dm, prefix='profile-')
dm %>% filter(source == "Tetraodon") %>% plot_pseudo("Tetraodon")

# plot all species, all bases
sapply(levels(dm$source), plot.save, 
       dm, 
       prefix='results/profiles-v4/pseudo/pseudo-')

# try to plot only the unmasked areas
sapply(levels(dm$source), plot.save, 
       dm %>% filter(compartment == 'normal'), 
       prefix='results/profiles-v4/pseudo-unmask/pseudo-unmask-')

sapply(levels(dm$source), plot.save, 
       dm %>% filter(compartment == 'masked'), 
       prefix='results/profiles-v4/pseudo-mask/pseudo-mask-')
