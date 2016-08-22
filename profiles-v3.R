setwd('c:/work/gar/')

read.name <- function(fn, datadir) {
  read.delim(paste0(datadir, "/", fn),
             col.names=c("seq", "pos", "A", "C", "G", "T"),
             colClasses=c("seq"="factor")) %>%
    mutate(sum = A + C + G + T, GC = (C + G) / sum, source=fn)
}

sources <- c(
  "Ciona_intestinalis.KH.gz.tsv" = "Sea squirt (Tunicate)",
  "Branchiostoma_floridae_v2.0.assembly.fasta.tsv" = "Lancelet",
  "Latimeria_chalumnae.LatCha1.dna_sm.toplevel.fa.gz.tsv" = "Latimeria",
  "Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz.tsv" = "Zebrafish",
  "Petromyzon_marinus.Pmarinus_7.0.dna_sm.toplevel.fa.gz.tsv" = "Lamprey",
  "calMil1.fa.gz.tsv" = "Elephant shark",
  "Takifugu_rubripes.FUGU4.dna_sm.toplevel.fa.gz.tsv" = "Fugu",
  "Gasterosteus_aculeatus.BROADS1.dna_sm.toplevel.fa.gz.tsv" = "Stickelback",
  "Gadus_morhua.gadMor1.dna_sm.toplevel.fa.gz.tsv" = "Cod",
  "BADN01.1.fsa_nt.gz.tsv" = "Tuna",
  "LepOcu1.fa.gz.tsv" = "Spotted gar",
  "Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz.tsv" = "Mouse",
  "Homo_sapiens.GRCh38.dna_sm.toplevel.fa.tsv" = "Human"
)

tabs <- lapply(names(sources), read.name, "data/profiles-v3")
da <- bind_rows(tabs)

da %>%
  mutate(source=factor(source, levels=names(sources), labels=sources)) %>%
  ggplot(aes(GC)) +
  geom_density(colour=NA, fill="#444444") +
  facet_wrap(~source, ncol=3) +
  xlim(c(0.3, 0.6)) +
  xlab("GC fraction") +
  ylab("block density") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
ggsave('results/profiles-v3/gc-summary.pdf', width=7.4, height=5, units="in")

