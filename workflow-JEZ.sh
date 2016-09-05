# additional work for revisions in JEZ

#
# download genomes
#
OUT=data-genomes/chicken
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-85/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna_sm.toplevel.fa.gz
EOF

OUT=data-genomes/stickleback
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-85/fasta/gasterosteus_aculeatus/dna/Gasterosteus_aculeatus.BROADS1.dna_sm.toplevel.fa.gz
ftp://ftp.ensembl.org/pub/release-85/fasta/gasterosteus_aculeatus/cds/Gasterosteus_aculeatus.BROADS1.cds.all.fa.gz
ftp://ftp.ensembl.org/pub/release-85/gff3/gasterosteus_aculeatus/Gasterosteus_aculeatus.BROADS1.85.gff3.gz
EOF

OUT=data-genomes/tetraodon
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-85/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.dna_sm.toplevel.fa.gz
ftp://ftp.ensembl.org/pub/release-85/fasta/tetraodon_nigroviridis/cds/Tetraodon_nigroviridis.TETRAODON8.cds.all.fa.gz
ftp://ftp.ensembl.org/pub/release-85/gff3/tetraodon_nigroviridis/Tetraodon_nigroviridis.TETRAODON8.85.gff3.gz
EOF

OUT=data-genomes/danio
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-85/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna_sm.toplevel.fa.gz
ftp://ftp.ensembl.org/pub/release-85/fasta/danio_rerio/cds/Danio_rerio.GRCz10.cds.all.fa.gz
ftp://ftp.ensembl.org/pub/release-85/gff3/danio_rerio/Danio_rerio.GRCz10.85.gff3.gz
EOF

# get chinese lancelet (not in Ensembl)
# urls captured from 'Inspect' in Chrome browser
# http://genome.bucm.edu.cn/lancelet/download_data.php
# weird chinese double slashed links..
OUT=data-genomes/b-belcheri
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
http://genome.bucm.edu.cn/download//datas//Branchiostoma.belcheri_v18h27.r3_ref_cds.fa.gz
http://genome.bucm.edu.cn/download//datas//Branchiostoma.belcheri_v18h27.r3_ref_genome.softmasked.fa.gz
http://genome.bucm.edu.cn/download//datas//Branchiostoma.belcheri_v18h27.r3_ref_protein.fa.gz
http://genome.bucm.edu.cn/download//datas//Branchiostoma.belcheri_v18h27.r3_ref_annotation.gff3.gz
EOF

# change the naming so it matches the Ensembl names..
mv $OUT/Branchiostoma.belcheri_v18h27.r3_ref_genome.softmasked.fa.bgz $OUT/Branchiostoma.belcheri_v18h27.r3.dna_sm.toplevel.fa.bgz

#
# index it with faidx
#
find data-genomes -name '*.dna_sm.*.bgz' |
  parallel -j1 samtools faidx {}

#
# calculate base counts
#
OUT=data/profiles-v4
mkdir -p $OUT
find data-genomes -name '*.dna_sm.*.fai' |
  parallel -u "python base_counts.py -q --window 10000 {.} > $OUT/{/.}.tsv"

#
# danio, tetraodon, stickleback
# calculate GC in regions
#
find data-genomes -name '*.gff3.bgz' |
    parallel "zcat {} | mawk -F'\t' '(\$3 == \"gene\")' | bedtools sort | bgzip > {.}.genes.bgz"

find data-genomes -name '*.gff3.genes.bgz' |
    parallel tabix -p gff {}


# use globbing to fing matching genome..
rmext () { echo ${1%%.*} ;}
export -f rmext

# test the mlutiple substitutions / subshells it first;)
find data-genomes -name '*.gff3.genes.bgz' |
    parallel "echo \$( rmext {} ).*.fai"

# bedtools complement dislikes our sorting..
find data-genomes -name '*.gff3.genes.bgz' | head -1 |
    parallel "bedtools complement -i {} -g <( <\$( rmext {} ).*.fai cut -f-2 | sort ) | bedtools sort | bgzip > {.}.complement.bgz"

#
# the same with bedops, hope they'll cooperate
# http://bedops.uwencode.org/forum/index.php?topic=19.0
# `bedops -c -L file.bed -`
#
# use globbing to fing matching genome..
rmext () { echo ${1%%.*} ;}
export -f rmext

# merge the gene regions
find data-genomes -name '*.gff3.bgz' |
    parallel "zcat {} | mawk -F'\t' '(\$3 == \"gene\")' | gff2bed | sort-bed - | bedops -m - | bgzip > \$( rmext {} ).genes.bed.bgz"

# index
find data-genomes -name '*.genes.bed.bgz' |
    parallel tabix -p bed {}

# complement
find data-genomes -name '*.genes.bed.bgz' |
    parallel "<{} zcat | bedops -c -L - <( <\$( rmext {} ).*.fai mawk '{print \$1, \$2, \$2+1}' ) | sort-bed - | bgzip > {.}.complement.bgz"

# index
find data-genomes -name '*.complement.bgz' |
    parallel tabix -p bed {}

# do the histograms for all the regions
find data-genomes -name '*.genes.*.bgz' |
  parallel "python base_hist_region.py \$( rmext {} ).*.toplevel.fa.bgz {} > {.}.hist.tsv"

##
# junkyard
##

# pysam is probably broken for accessing compressed fasta
python base_hist_region.py LepOcu1.fa.gz gtf/sorted-genes.gz > LepOcu-hist-genes.csv

# extract the stats
python base_hist_region.py LepOcu1.fa gtf/sorted-genes.gz > LepOcu-hist-genes.csv

