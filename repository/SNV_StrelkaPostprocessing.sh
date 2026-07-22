#!/bin/bash

##########################################################################################
##
## SNV_StrelkaPostprocessing.sh
##
## Postprocessing for Strelka.
##
##########################################################################################

name=$1
species=$2
config_file=$3
filtering=$4

source $config_file

echo '---- Strelka Postprocessing I (Indel size selection, filtering) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gunzip $name/results/Strelka/Strelka/results/variants/somatic.snvs.vcf.gz
gunzip $name/results/Strelka/Strelka/results/variants/somatic.indels.vcf.gz

cat $name/results/Strelka/Strelka/results/variants/somatic.snvs.vcf | \
SnpSift filter \
	"( ( ( FILTER = 'PASS' ) & ( ALT = 'T' ) & ( GEN[NORMAL].TU[0] <= 1 ) & ( GEN[TUMOR].TU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].TU[0] ) ) \
	| ( ( FILTER = 'PASS' ) & ( ALT = 'A' ) & ( GEN[NORMAL].AU[0] <= 1 ) & ( GEN[TUMOR].AU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].AU[0] ) ) \
	| ( ( FILTER = 'PASS' ) & ( ALT = 'C' ) & ( GEN[NORMAL].CU[0]<= 1 ) & ( GEN[TUMOR].CU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].CU[0] ) ) \
	| ( ( FILTER = 'PASS' ) & ( ALT = 'G' ) & ( GEN[NORMAL].GU[0] <= 1 ) & ( GEN[TUMOR].GU[0] >= 2 ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].GU[0] ) ) )" \
	> $name/results/Strelka/$name.str.snp.postprocessed.vcf

cat $name/results/Strelka/Strelka/results/variants/somatic.indels.vcf | \
SnpSift filter \
	"( ( FILTER = 'PASS' ) & ( GEN[TUMOR].DP >= 5 ) & ( GEN[NORMAL].DP >= 5 ) & ( GEN[NORMAL].TIR[0] <= 1) & \
	( GEN[TUMOR].TIR[0] >= 2 ) & ( GEN[TUMOR].DP * 0.05 <= GEN[TUMOR].TIR[0] ) )" \
	> $name/results/Strelka/$name.str.indel.filtered.vcf

gatk --java-options '-Xmx256G' SelectVariants --max-indel-size 10 \
	-V $name/results/Strelka/"$name".str.indel.filtered.vcf \
	-O $name/results/Strelka/"$name".str.indel.postprocessed.vcf


echo '---- Strelka Postprocessing II (Filtering out known SNV/Indel using dbSNP or the Sanger Mouse database) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

for method in indel snp; 
	do
	(
		bgzip $name/results/Strelka/"$name".str.$method.postprocessed.vcf
		tabix -p vcf $name/results/Strelka/"$name".str.$method.postprocessed.vcf.gz

		if [ $filtering = 'all' ] || [ $filtering = 'hard' ]; then
			bcftools isec -C -c none -O z -w 1 \
				-o $name/results/Strelka/"$name".str.$method.postprocessed.snp_removed.vcf.gz \
				$name/results/Strelka/"$name".str.$method.postprocessed.vcf.gz \
				$snp_file
		elif [ $filtering = 'none' ]; then
			cp $name/results/Strelka/"$name".str.$method.postprocessed.vcf.gz \
			$name/results/Strelka/"$name".str.$method.postprocessed.snp_removed.vcf.gz
		fi

		mv $name/results/Strelka/"$name".str.$method.postprocessed.snp_removed.vcf.gz \
			$name/results/Strelka/"$name".Strelka.$method.vcf.gz
		gunzip $name/results/Strelka/"$name".Strelka.$method.vcf.gz
	) &
done

wait


echo '---- Strelka Postprocessing III (Extracting allele frequencies) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'Human' ]; then
	species=homo_sapiens
	assembly=GRCh38
elif [ $species = 'Mouse' ]; then
	species=mus_musculus
	assembly=GRCm38
fi
cache_dir=$genome_dir/VEP/

for method in snp indel; do
	(
	vep --cache --species $species \
		-i $name/results/Strelka/$name.Strelka.$method.vcf \
		-o $name/results/Strelka/$name.Strelka.$method.vep.vcf \
		--fasta $genome_file --assembly $assembly \
		--offline --no_progress --no_stats \
		--buffer_size 5000 --sift b --ccds --uniprot --hgvs \
		--symbol --numbers --domains --gene_phenotype --canonical \
		--protein --biotype --uniprot --tsl --pubmed --variant_class \
		--shift_hgvs 1 --check_existing --total_length --allele_number \
		--no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele \
		--pick_order canonical,tsl,biotype,rank,ccds,length --format vcf \
		--fork 4 --cache_version 96 --polyphen b --af --af_1kg --af_esp \
		--af_gnomad --force_overwrite --dir_cache $cache_dir

	vcf2maf.pl \
		--input-vcf $name/results/Strelka/$name.Strelka.$method.vcf \
		--output-maf $name/results/Strelka/$name.Strelka.$method.vep.maf \
		--vep-path $vep_dir --vep-data $vepdata_dir --vep-overwrite \
		--ncbi-build $assembly --ref-fasta $genome_file \
		--vcf-normal-id NORMAL --vcf-tumor-id TUMOR \
		--tumor-id $name --species $species

	maf2vcf.pl \
		--input-maf $name/results/Strelka/$name.Strelka.$method.vep.maf \
		--output-dir $name/results/Strelka/ \
		--output-vcf $name/results/Strelka/$name.Strelka.$method.maf.vcf \
		--ref-fasta $genome_file --per-tn-vcfs 

	sed -i "s/$name\\tNORMAL/TUMOR\\tNORMAL/g" $name/results/Strelka/$name.Strelka.$method.maf.vcf
	bgzip $name/results/Strelka/$name.Strelka.$method.maf.vcf
	tabix -p vcf $name/results/Strelka/$name.Strelka.$method.maf.vcf.gz
	) &
done

wait
bcftools concat -a $name/results/Strelka/$name.Strelka.snp.maf.vcf.gz $name/results/Strelka/$name.Strelka.indel.maf.vcf.gz > $name/results/Strelka/$name.Strelka.vcf


echo '---- Strelka Postprocessing IV (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'homo_sapiens' ]; then
	SnpSift annotate \
		$dbsnp_file $name/results/Strelka/$name.Strelka.vcf \
		> $name/results/Strelka/$name.Strelka.ann1.vcf

	SnpSift annotate \
		$cosmiccoding_file $name/results/Strelka/$name.Strelka.ann1.vcf \
		> $name/results/Strelka/$name.Strelka.ann2.vcf

	SnpSift annotate \
		$cosmicnoncoding_file $name/results/Strelka/$name.Strelka.ann2.vcf \
		> $name/results/Strelka/$name.Strelka.ann3.vcf

	SnpSift annotate \
		$clinvar_file $name/results/Strelka/$name.Strelka.ann3.vcf \
		> $name/results/Strelka/$name.Strelka.ann4.vcf

	SnpSift DbNSFP \
		-db $dbnsfp_file -v -f MetaLR_pred,MetaSVM_pred,SIFT_pred,PROVEAN_pred \
		$name/results/Strelka/$name.Strelka.ann4.vcf \
		> $name/results/Strelka/"$name".Strelka.ann5.vcf

	snpEff $snpeff_version -canon \
		-csvStats $name/results/Strelka/$name.Strelka.annotated.vcf.stats \
		$name/results/Strelka/"$name".Strelka.ann5.vcf \
		> $name/results/Strelka/$name.Strelka.annotated.vcf

	cat $name/results/Strelka/$name.Strelka.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Strelka/$name.Strelka.annotated.one.vcf

	SnpSift extractFields \
		$name/results/Strelka/$name.Strelka.annotated.one.vcf \
		CHROM POS REF ALT "GEN[TUMOR].AD[0]" "GEN[TUMOR].AD[1]" \
		"GEN[NORMAL].AD[0]" "GEN[NORMAL].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
		ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF GENOME_SCREEN_SAMPLE_COUNT \
		SAMPLE_COUNT CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred dbNSFP_PROVEAN_pred \
		> $name/results/Strelka/"$name".Strelka.txt

elif [ $species = 'mus_musculus' ]; then
	snpEff $snpeff_version -canon \
		-csvStats $name/results/Strelka/$name.Strelka.annotated.vcf.stats \
		$name/results/Strelka/$name.Strelka.vcf \
		> $name/results/Strelka/$name.Strelka.annotated.vcf

	cat $name/results/Strelka/$name.Strelka.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Strelka/$name.Strelka.annotated.one.vcf

	SnpSift extractFields \
		$name/results/Strelka/$name.Strelka.annotated.one.vcf \
		CHROM POS REF ALT "GEN[TUMOR].AD[0]" "GEN[TUMOR].AD[1]" \
		"GEN[NORMAL].AD[0]" "GEN[NORMAL].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
		ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		> $name/results/Strelka/$name.Strelka.txt
fi

awk 'BEGIN{FS=OFS="\t"} FNR==1{print $1,$2,$3,$4,"GEN[TUMOR].AF",substr($0, index($0,$5)),"DP"}' $name/results/Strelka/$name.Strelka.txt > $name/results/Strelka/$name.Strelka.head
awk 'BEGIN{FS=OFS="\t"} NR>1{print $0,$NF=$5+$6}' $name/results/Strelka/$name.Strelka.txt > $name/results/Strelka/$name.Strelka.cut
awk 'BEGIN{FS=OFS="\t"} {$5=$6/$NF"\t"$5}1' $name/results/Strelka/$name.Strelka.cut > $name/results/Strelka/$name.Strelka.tmp
cat $name/results/Strelka/$name.Strelka.head $name/results/Strelka/$name.Strelka.tmp > $name/results/Strelka/$name.Strelka.txt

sh $repository_dir/SNV_CleanUp.sh $name Strelka