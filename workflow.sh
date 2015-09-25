# base_counts.py needs .fai
# index the newcoming genomes
parallel -j1 "/opt/samtools-0.1.18/samtools faidx genomes/{}/*.gz" ::: cod fugu stickelback tuna

# some fixes with human assembly
(<Homo_sapiens.GRCh38.dna_sm.toplevel.fa.fai-orig head -22;cat Homo_sapiens.GRCh38.dna_sm.toplevel.fa.fai-XY) > Homo_sapiens.GRCh38.dna_sm.toplevel.fa.fai

# check if all fai are ok (more than one record)
find genomes -name '*.fai'| sort -t'/' -k2,2 | xargs wc -l

# run base_counts.py for all fai
find genomes -name '*.fai'|parallel -u "python base_counts.py {.} 10000 > data/{/.}.tsv"

# quick check - compare number of 10k blocks with size of the genome
# chromosomes sizes from .fai
find genomes -name '*.fai'|xargs mawk '{print FILENAME "\t" $1 "\t" $2}'|cut -d'/' -f3 | sed 's/.fa.gz.fai//' > data/chrom-sizes.tsv
# count of data blocks
wc -l data/* | head -n -1 | mawk '{print $2 "\t" $1}' | cut -d'/' -f2 | sed 's/.fa.gz.tsv//' > data/block-counts.tsv

# summary of chromosome sizes
<data/chrom-sizes.tsv mawk '{m[$1] += $3} END{for (x in m) {print x, m[x]}}'

