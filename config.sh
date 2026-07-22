#!/bin/bash
##########################################################################################
##
## config.sh
##
##########################################################################################

repository_dir=$(whereis MoCaSeq.sh | awk '{print $2}' | sed 's/\/[^\/]*$//')/repository
genomes_dir=$(whereis MoCaSeq.sh | awk '{print $2}' | sed 's/\/[^\/]*$//')/ref

discvrseq_dir=/fs/pool/pool-mann-projects/MaxZ/conda_envs/genomeseq/bin
vep_dir=/fs/pool/pool-mann-projects/MaxZ/conda_envs/genomeseq/share/ensembl-vep-96.3-0

if [ $species = 'Mouse' ]; then
	genome_dir=$genomes_dir/GRCm38.p6
	snp_file=$genome_dir/MGP.v5.snp_and_indels.exclude_wild.vcf.gz
	genome_file=$genome_dir/GRCm38.p6.fna
	genomeindex_dir=$genome_dir/bwa_index/GRCm38.p6
	interval_file=$genome_dir/GRCm38.SureSelect_Mouse_All_Exon_V1.bed.list
	bammatcher_file=$genome_dir/GRCm38.bammatcher_bash.conf
	snpeff_version=GRCm38.86
	microsatellite_file=$genome_dir/GRCm38.p6.microsatellites
	callregions_file=$genome_dir/GRCm38.p6.canonical_chromosomes.bed.gz
	CGC_file=$genome_dir/GRCm38.Census_allMon_Jan_15_11_46_18_2018_mouse.tsv
	TruSight_file="NULL"
	chromosomes=19
	gcWig_file=$genome_dir/GRCm38.p6.gc.20000.wig
	mapWig_file=$genome_dir/GRCm38.p6.map.20000.wig
	exons_file=$genome_dir/GRCm38.SureSelect_Mouse_All_Exon_V1.bed
	centromere_file="NULL"
	varregions_file=$genome_dir/GRCm38.AgilentProbeGaps.txt
	gencode_file_exons=$genome_dir/GRCm38.gencode_M20_Exons.rds
	gencode_file_genes=$genome_dir/GRCm38.gencode_M20_Genes.rds
	vepdata_dir=$genome_dir/VEP
	dbsnp_file=$genome_dir/
	cosmiccoding_file=$genome_dir/
	cosmicnoncoding_file=$genome_dir/
	clinvar_file=$genome_dir/
	dbnsfp_file=$genome_dir/

elif [ $species = 'Human' ]; then
	genome_dir=$genomes_dir/GRCh38.p14
	snp_file=$genome_dir/GRCh38.p14.dbSNP.vcf.gz
	genome_file=$genome_dir/GRCh38.p14.fna
	genomeindex_dir=$genome_dir/bwa_index/GRCh38.p14
	interval_file=$genome_dir/GRCh38.SureSelectXT_Human_All_Exon_V8_Plus_NCV.bed
	bammatcher_file=$genome_dir/GRCh38.bammatcher_bash.conf
	snpeff_version=GRCh38.86
	microsatellite_file=$genome_dir/GRCh38.p14.microsatellites
	callregions_file=$genome_dir/GRCh38.p14.canonical_chromosomes.bed.gz
	CGC_file="NULL"
	TruSight_file="NULL"
	chromosomes=22
	gcWig_file=$genome_dir/GRCh38.p14.gc.20000.wig
	mapWig_file=$genome_dir/GRCh38.p14.map.20000.wig
	exons_file=$genome_dir/GRCh38.SureSelectXT_Human_All_Exon_V8.bed
	centromere_file=$genome_dir/GRCh38.p14.centromeres.bed
	varregions_file="NULL"
	gencode_file_exons=$genome_dir/GRCh38.p14.gencode.v46.Exons.rds
	gencode_file_genes=$genome_dir/GRCh38.p14.gencode.v46.Genes.rds
	vepdata_dir=$genome_dir/VEP
	dbsnp_file=$genome_dir/GRCh38.p14.dbSNP.vcf.gz
	cosmiccoding_file=$genome_dir/GRCh38.p14.cosmicCoding.vcf.gz
	cosmicnoncoding_file=$genome_dir/GRCh38.p14.cosmicNonCoding.vcf.gz
	clinvar_file=$genome_dir/GRCh38.p14.clinvar.vcf.gz
	dbnsfp_file=$genome_dir/GRCh38.p14.dbNSFP.txt.gz
fi
