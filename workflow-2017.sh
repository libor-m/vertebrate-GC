# more GC profiles to compare with Laurent's work

OUT=data-genomes/tilapia
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-88/fasta/oreochromis_niloticus/dna/Oreochromis_niloticus.Orenil1.0.dna_sm.toplevel.fa.gz
EOF

OUT=data-genomes/tilapia2
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
http://cichlid.umd.edu/download/O_niloticus_UMD1/O_niloticus_UMD1.fasta.gz
EOF

OUT=data-genomes/medaka
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ensembl.org/pub/release-88/fasta/oryzias_latipes/dna/Oryzias_latipes.MEDAKA1.dna_sm.toplevel.fa.gz
EOF

# this is too slow..
# OUT=data-genomes/fugu
# mkdir -p $OUT
# parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
# http://www.fugu-sg.org/downloads/data/fugu5_chromosomes.fa.gz
# EOF

OUT=data-genomes/fugu
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
http://hgdownload.cse.ucsc.edu/goldenPath/fr3/bigZips/fr3.fa.gz
EOF

OUT=data-genomes/pike
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/721/915/GCF_000721915.3_Eluc_V3/GCF_000721915.3_Eluc_V3_genomic.fna.gz
EOF

OUT=data-genomes/cod
mkdir -p $OUT
parallel "wget -q -O - {} | zcat | bgzip > $OUT/{/.}.bgz" <<EOF
https://ndownloader.figshare.com/files/5323414
EOF

#
# index it with faidx
#
find data-genomes -name '*.bgz' |
  parallel -j1 samtools faidx {}

# setup env with anaconda python
# skip if already done
conda create -n gar python=2
. activate gar
conda install -c bioconda pysam

# activate anyways
. activate gar

#
# calculate base counts
#
OUT=data/profiles-2017
mkdir -p $OUT
find data-genomes -name '*.fai' |
  parallel -u "python base_counts.py -q --window 10000 {.} > $OUT/{/.}.tsv"
