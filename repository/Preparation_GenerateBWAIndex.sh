#!/bin/bash

##########################################################################################
##
## Preparation_GenerateBWAIndex.sh
##
## Uses a downloaded .fasta-file to generate the BWA index needed for contig-aware mapping.
##
##########################################################################################

species=$1
version=$2

if [ $species = 'Mouse' ]; then
	sed -ri "s/>CM[0-9\.]* Mus musculus chromosome ([0-9XY]*).*/>\1/g" ref/$version/$version.fna
elif [ $species = 'Human' ]; then
	sed -ri "s/>CM[0-9\.]* Homo sapiens chromosome ([0-9XY]*).*/>\1/g" ref/$version/$version.fna
fi
samtools faidx ref/$version/$version.fna
samtools faidx ref/$version/$version.fna $(grep -P -o '^>[A-Z].*\.\d*' ref/$version/$version.fna | sed 's/>//g') > ref/$version/$version"_alt.fna"
samtools faidx ref/$version/$version.fna $(grep -P -o '>[0-9XY]+' ref/$version/$version.fna | sed 's/>//g') > ref/$version/$version"_primary.fna"

fasta_to_fastq.pl ref/$version/$version"_alt.fna" > ref/$version/haplotypes.fastq

awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ref/$version/$version"_primary.fna" | sed 's/^>//' | awk '{if(NR%2) printf "%s\t1\t", $0; else print $0}' > ref/$version/$version.canonical_chromosomes.bed
bgzip ref/$version/$version.canonical_chromosomes.bed
tabix -p bed ref/$version/$version.canonical_chromosomes.bed.gz


mkdir -p ref/$version/primary_index
bwa index -p ref/$version/primary_index/bwa_index -a bwtsw ref/$version/$version"_primary.fna"

bwa mem ref/$version/primary_index/bwa_index ref/$version/haplotypes.fastq | \
	samtools sort -O sam - | \ 
	samtools view -h -F 0x800 - > ref/$version/alt_mapping.sam

mkdir -p ref/$version/bwa_index
grep '^@SQ' ref/$version/alt_mapping.sam > ref/$version/header
grep '^[^@]' ref/$version/alt_mapping.sam > ref/$version/data

paste <(cut -f 1 ref/$version/data) \
	  <(yes 0 | head -n $(wc -l < ref/$version/data)) \
	  <(cut -f 3-4 ref/$version/data) \
	  <(yes 255 | head -n $(wc -l < ref/$version/data)) \
	  <(cut -f 6 ref/$version/data) \
	  <(yes '*' | head -n $(wc -l < ref/$version/data)) \
	  <(yes 0 | head -n $(wc -l < ref/$version/data)) \
	  <(yes 0 | head -n $(wc -l < ref/$version/data)) \
	  <(yes '*' | head -n $(wc -l < ref/$version/data)) \
	  <(yes '*' | head -n $(wc -l < ref/$version/data)) \
	  <(grep -P -o 'NM:i:\d*' ref/$version/data) | \
	cat ref/$version/header - > ref/$version/bwa_index/$version.alt

bwa index -p ref/$version/bwa_index/$version -a bwtsw ref/$version/$version.fna

rm ref/$version/$version"_primary.fna"
rm ref/$version/$version"_alt.fna"
rm ref/$version/haplotypes.fastq
rm ref/$version/alt_mapping.sam
rm -r ref/$version/primary_index
rm ref/$version/data
rm ref/$version/header