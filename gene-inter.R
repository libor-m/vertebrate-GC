library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)

##
# repeated code into functions
##
load_data <- function(path)
    bind_rows(read_tsv(paste0(path, ".genes.bed.hist.tsv"),
                       col_names = c("seq", "start", "end", "idx", "base", "count"),
                       col_types = "ciiici") %>%
                  mutate(compartment="gene"),
              read_tsv(paste0(path, ".genes.bed.complement.hist.tsv"),
                       col_names = c("seq", "start", "end", "idx", "base", "count"),
                       col_types = "ciiici") %>%
                  mutate(compartment="intergenic"))

long2wide <- function(d)
    d %>%
        spread(base, count, fill = 0) %>%
        mutate(masked = a + g + c + t + n,
               gc_masked = g + c,
               gc_total = G + C + g + c,
               total = a + A + g + G + c + C + t + T + n, 
               normal_total = A + C + G + T,
               gc_normal = G + C) 

densplot <- function(d)
    d %>%
      ggplot(aes(gc_total/total, fill=compartment)) +
      geom_density(alpha=0.6, colour=NA) +
      xlim(0, 0.8)

load_plot_save <- function(path) {
    load_data(path) %>%
        long2wide %>%
        densplot
    ggsave(paste0('results/gc-dens/gc-dens-', 
                  strsplit(path, "/")[[1]] %>% tail(1),
                  ".pdf"),
           width = 6,
           height = 5)
}
    

# test it with danio
load_data("data-genomes/danio/Danio_rerio") -> d_danio
long2wide(d_danio) -> d_danio_wide
d_danio_wide %>% densplot
ggsave('results/gc-dens-Danio_rerio.pdf', width = 6, height = 5)

# all at once
load_plot_save("data-genomes/danio/Danio_rerio")
load_plot_save("data-genomes/b-belcheri/Branchiostoma")
load_plot_save("data-genomes/tetraodon/Tetraodon_nigroviridis")
load_plot_save("data-genomes/stickleback/Gasterosteus_aculeatus")
