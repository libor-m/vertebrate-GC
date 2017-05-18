library(tidyverse)

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

c("5323414.bgz.tsv" = "Cod",
  "fr3.fa.bgz.tsv" = "Fugu",
  "GCF_000721915.3_Eluc_V3_genomic.fna.bgz.tsv" = "Pike",
  "O_niloticus_UMD1.fasta.bgz.tsv" = "Tilapia",
  "Oryzias_latipes.MEDAKA1.dna_sm.toplevel.fa.bgz.tsv" = "Medaka") ->
 sources

# load "data all"
lapply(names(sources), read.name, "data/profiles-2017") %>%
    bind_rows ->
    da

# separate profiles for masked and normal areas
# (i. e. isolate variable 'masked' from column names)
# and calculate few more handy columns
wide_to_long <- function(d, source_names)
    d %>%
    select(-contains("GC"), -contains("sum")) %>%
    gather(base, count, A:t) %>%
    mutate(compartment=ifelse(base %in% c('A', 'C', 'G', 'T'), "normal", "masked")) %>%
    mutate(base=toupper(base)) %>%
    spread(base, count) %>%
    mutate(sum=A + C + G + T,
           GC=(C + G) / sum,
           source=factor(source, levels=names(source_names), labels=source_names))

da %>% wide_to_long(sources) -> dm

source('chromoplot.R')

# improved plot.save which:
# - tries to keep the height of the facet constant
# - uses the plot data to find number of facets
# - filters out too short scaffolds
plot.save2 <- function(species,
                       din,
                       prefix="",
                       bases_per_pseudo=5e8,
                       bases_per_block=1e5,
                       min_scaffold=3e5,
                       panel_height=25) {
    din %>%
        filter(source == species) ->
        d
    
    d %>%
        group_by(seq) %>%
        summarise(len=max(pos)) %>%
        ungroup %>%
        filter(len > min_scaffold) %>%
        .$seq ->
        long_seqs
    
    d %>%
        filter(seq %in% long_seqs,
               !grepl("rand", seq)) %>%
        plot_pseudo(species,
                    bases_per_pseudo,
                    chromsort=chrom_sort_alnum) ->
        the_plot
    
    npanels <- the_plot$data$pseudo_seq %>% levels %>% length
    
    # total height =
    #  npanels * panel_height +
    #  (npanels-1) * panel_spacing +
    #  extra
    # where panel_spacing and extra should be constant
    # measuring few outputs:
    # panel_spacing seems to be 1.2 mm
    # extra is 11.1 mm for the top margin, 13.1 for the bottom margin
    h <- npanels * panel_height + (npanels - 1) * 1.2 + 24.2
    
    ggsave(paste0(prefix, species, '.pdf'),
           plot=the_plot,
           width=290, height=h, units = "mm")
}


plot.save2("Cod", dm, prefix = "results/profiles-2017/profile-", min_scaffold=1e6)

plot.save2("Fugu", dm, prefix="results/profiles-2017/profile-", min_scaffold=1e6)
plot.save2("Pike", dm, prefix="results/profiles-2017/profile-", min_scaffold=1e6)
plot.save2("Tilapia", dm, prefix="results/profiles-2017/profile-", min_scaffold=1e6)
plot.save2("Medaka", dm, prefix="results/profiles-2017/profile-", min_scaffold=1e6)

