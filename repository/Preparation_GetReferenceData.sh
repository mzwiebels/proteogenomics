#!/bin/bash

##########################################################################################
##
## Preparation_GetReferenceDataMouse.sh
##
## Main routine for the download of all reference data needed for the WES and WGS workflows.
##
##########################################################################################

config_file=$1
temp_dir=$2
species=$3
version=$4


#reading configuration from $config_file
. $config_file

mkdir -p ref
mkdir -p ref/$version
mkdir -p ref/$version/VEP

if [ $species = 'Mouse' ]; then
	GCAdir=$(curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/ | grep -o "\S*$version\S*")
	GCAfile=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/$GCAdir/$GCAdir"_genomic.fna.gz"
	regionsFile=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/$GCAdir/$GCAdir"_assembly_structure/genomic_regions_definitions.txt"

	VEPfile="ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/mus_musculus_vep_96_"${version%.*}".tar.gz"
	taxName=mus_musculus

	gencodeDir=$(curl -s https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/ | grep -oP 'gencode\.vM\d*' | uniq)
	gencodeFile=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/$gencodeDir.annotation.gtf.gz

	chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y"
elif [ $species = 'Human' ]; then
	GCAdir=$(curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/ | grep -o "\S*$version\S*")
	GCAfile=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$GCAdir/$GCAdir"_genomic.fna.gz"
	regionsFile=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/$GCAdir/$GCAdir"_assembly_structure/genomic_regions_definitions.txt"

	VEPfile="ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/homo_sapiens_vep_96_"${version%.*}".tar.gz"
	taxName=homo_sapiens

	gencodeDir=$(curl -s https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/ | grep -oP 'gencode\.v\d*' | uniq)
	gencodeFile=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/$gencodeDir.annotation.gtf.gz

	chromosomes="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
fi

echo '---- Get reference data ----' | tee ref/$version/GetReferenceData.txt
echo '---- Generate reference data for Version '$version' ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

#rerouting STDERR to report file
exec 2>> ref/$version/GetReferenceData.txt

echo '---- Copying over files from repository ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

if [ $species = 'Mouse' ]; then
	cp $repository_dir/../data/GRCm38.Census_allMon_Jan_15_11_46_18_2018_mouse.tsv ref/$version

	cp $repository_dir/../data/GRCm38.bammatcher_bash.conf ref/$version

	cp $repository_dir/../data/GRCm38.AgilentProbeGaps.txt ref/$version
	cp $repository_dir/../data/GRCm38.Genecode_M20_Exons.rds ref/$version
	cp $repository_dir/../data/GRCm38.Genecode_M20_Genes.rds ref/$version

	cp $repository_dir/../data/GRCm38.RefFlat ref/$version
elif [ $species = 'Human' ]; then
	cp $repository_dir/../data/GRCh38.bammatcher_bash.conf ref/$version

	#cp $repository_dir/../data/GRCm38.Genecode_M20_Exons.rds ref/$version
	#cp $repository_dir/../data/GRCm38.Genecode_M20_Genes.rds ref/$version
fi


echo '---- Downloading reference genome ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

wget -nv -c -P ref/$version $GCAfile
gunzip ref/$version/$(basename $GCAfile)

wget -nv -c -P ref/$version $regionsFile
grep '^CEN' ref/$version/genomic_regions_definitions.txt | cut -f 1,3,4 | sed -r 's/CEN([0-9XY]*)/\1/g' > ref/$version/$version.centromeres.bed
rm ref/$version/genomic_regions_definitions.txt


echo '---- Generate BWA Index ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

sh $repository_dir/Preparation_GenerateBWAIndex.sh $species $version


echo '---- Generate sequence dictionary ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

picard CreateSequenceDictionary \
	O=ref/$version/$version.dict \
	R=ref/$version/$version.fna

echo '---- Generate exons covered by SureSelect ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt
# Download specific .zip-file from https://earray.chem.agilent.com/suredesign/search.htm
# Attention: Most use old versions (mm9/hg19) -> they need liftover
# use "_Regions.bed" for further work - this covers all regions which are targeted in this kit
# Liftover with https://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core using default settings
# Rename to .bed and move to main reference directory
# Included in the data-directory is a version which has already been lifted over - nothing more to do but generating the sequence dictionary
cp $repository_dir/../data/$(basename $exons_file) ref/$version/
picard BedToIntervalList \
	I=ref/$version/$(basename $exons_file) \
	O=ref/$version/$(basename $exons_file).list \
	SD=ref/$version/$version.dict


echo '---- Downloading reference genome (for VEP) ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

wget -nv -c -P ref/$version/VEP $VEPfile
tar -xzf ref/$version/VEP/$(basename $VEPfile)
mv $taxName ref/$version/VEP
rm ref/$version/VEP/$(basename $VEPfile)


if [ $species = 'Mouse' ]; then
	echo '---- Generate customized Sanger DB ----' | tee -a ref/$version/GetReferenceData.txt
	date | tee -a ref/$version/GetReferenceData.txt

	sh $repository_dir/Preparation_GenerateSangerMouseDB.sh $version $temp_dir
	rm -rf $temp_dir

	wget -nv -c -P ref/"$version"/ ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz
	wget -nv -c -P ref/"$version"/ ftp://ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz.tbi
	mv ref/"$version"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz ref/"$version"/MGP.v6.snp_and_indels.vcf.gz
	mv ref/"$version"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz.tbi ref/"$version"/MGP.v6.snp_and_indels.vcf.gz.tbi
	rm -r ref/"$version"/ftp-mouse.sanger.ac.uk/REL-1807-SNPs_Indels/
	rm -r ref/"$version"/ftp-mouse.sanger.ac.uk/

	bcftools view -s ^CAST_EiJ,SPRET_EiJ,PWK_PhJ,WSB_EiJ,MOLF_EiJ,ZALENDE_EiJ,LEWES_EiJ --min-ac=1 --no-update ref/"$version"/MGP.v6.snp_and_indels.vcf.gz -O z -o ref/"$version"/MGP.v6.snp_and_indels.exclude_wild.vcf.gz
	tabix -p vcf ref/"$version"/MGP.v6.snp_and_indels.exclude_wild.vcf.gz

	rm "ref/"$version"/MGP.v6.snp_and_indels.vcf.gz"
	rm "ref/"$version"/MGP.v6.snp_and_indels.vcf.gz.tbi"
elif [ $species = 'Human' ]; then
	echo '---- Downloading human SNP data (from dbSNP) ----' | tee -a ref/$version/GetReferenceData.txt
	date | tee -a ref/$version/GetReferenceData.txt

	wget -nv -c -P ref/$version https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
	wget -nv -c -P ref/$version https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi
	mv ref/$version/00-common_all.vcf.gz ref/$version/$version.dbSNP.vcf.gz
	mv ref/$version/00-common_all.vcf.gz.tbi ref/$version/$version.dbSNP.vcf.gz.tbi
fi


echo '---- Downloading Gencode Annotation ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

wget -nv -c -P ref/$version $gencodeFile
gunzip ref/$version/$(basename $gencodeFile) | grep -v "^#" > ref/$version/$version.$gencodeDir.gtf

Rscript $repository_dir/Preparation_AnnotateGencode.R ref/$version/$version.$gencodeDir.gtf $version $gencodeDir
rm ref/$version/$(basename $gencodeFile .gz) ref/$version/$version.$gencodeDir.gtf


echo '---- Generate reference data for msisensor  ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

msisensor scan -d ref/$version/$version.fna -o ref/$version/$version.microsatellites


echo '---- Optional for WES: Generating reference data for CopywriteR ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

Rscript $repository_dir/Preparation_GenerateCopywriterReferences.R $species $version


echo '---- Optional for WGS: Generating reference data for HMMCopy ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt

generateMap.pl -o ref/$version/$version.fna.map.bw -b -w 150 -c $chromosomes ref/$version/$version.fna
generateMap.pl -o ref/$version/$version.fna.map.bw -w 150 -c $chromosomes ref/$version/$version.fna

for resolution in 1000 10000 20000 50000; 
do
	gcCounter -w $resolution -c $chromosomes ref/$version/$version.fna > ref/$version/$version.gc.$resolution.wig
	mapCounter -w $resolution -c $chromosomes ref/$version/$version.fna.map.bw > ref/$version/$version.map.$resolution.wig
done
rm ref/$version/*.ebwt

echo '---- Finished generating reference data ----' | tee -a ref/$version/GetReferenceData.txt
date | tee -a ref/$version/GetReferenceData.txt
echo 'DONE' | tee -a ref/$version/GetReferenceData.txt

exit 0