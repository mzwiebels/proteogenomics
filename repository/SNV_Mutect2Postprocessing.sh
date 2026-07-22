#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Postprocessing.sh
##
## Postprocessing for Mutect2.
##
##########################################################################################

name=$1
species=$2
config_file=$3
filtering=$4

. $config_file

echo '---- Mutect2 Postprocessing I (OrientationFilter, Indel size selection, filtering) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gatk --java-options '-Xmx256G' LearnReadOrientationModel \
	-I $name/results/Mutect2/$name.m2.f1r2.tar.gz \
	-O $name/results/Mutect2/$name.m2.read-orientation-model.tar.gz

gatk --java-options '-Xmx256G' FilterMutectCalls \
	--variant $name/results/Mutect2/"$name".m2.vcf \
	--output $name/results/Mutect2/"$name".m2.filt.vcf \
	--reference $genome_file \
	--ob-priors $name/results/Mutect2/$name.m2.read-orientation-model.tar.gz

gatk --java-options '-Xmx256G' SelectVariants --max-indel-size 10 \
	-V $name/results/Mutect2/$name.m2.filt.vcf \
	-output $name/results/Mutect2/$name.m2.filt.selected.vcf

if [ $filtering = 'all' ]; then
	cat $name/results/Mutect2/$name.m2.filt.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.05) & \
		( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 5 ) & \
		( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 5 ) & \
		(GEN[Tumor].AD[1] >= 2) & (GEN[Normal].AD[1] <= 1) )" \
		> $name/results/Mutect2/$name.m2.postprocessed.vcf
elif [ $filtering = 'hard' ]; then
	cat $name/results/Mutect2/$name.m2.filt.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS') & (GEN[Tumor].AF >= 0.1) & \
		( ( GEN[Tumor].AD[0] + GEN[Tumor].AD[1]) >= 10 ) & \
		( ( GEN[Normal].AD[0] + GEN[Normal].AD[1]) >= 10 ) & \
		(GEN[Tumor].AD[1] >= 3) & (GEN[Normal].AD[1] = 0) )" \
		> $name/results/Mutect2/$name.m2.postprocessed.vcf
elif [ $filtering = 'none' ]; then
	cat $name/results/Mutect2/$name.m2.filt.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS' ) )" \
		> $name/results/Mutect2/$name.m2.postprocessed.vcf
fi

echo '---- Mutect2 Postprocessing II (Filtering out known SNV/Indel using dbSNP or the Sanger Mouse database) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

bgzip -f $name/results/Mutect2/$name.m2.postprocessed.vcf
tabix -p vcf $name/results/Mutect2/$name.m2.postprocessed.vcf.gz

if [ $filtering = 'all' ] || [ $filtering = 'hard' ]; then
	bcftools isec -C -c none -O z -w 1 \
		-o $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
		$name/results/Mutect2/$name.m2.postprocessed.vcf.gz \
		$snp_file
elif [ $filtering = 'none' ]; then
	cp $name/results/Mutect2/$name.m2.postprocessed.vcf.gz $name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz
fi

bcftools norm -m -any \
	$name/results/Mutect2/$name.m2.postprocessed.snp_removed.vcf.gz \
	-O z -o $name/results/Mutect2/$name.Mutect2.vcf.gz

gunzip -f $name/results/Mutect2/$name.Mutect2.vcf.gz


echo '---- Mutect2 Postprocessing III (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'Human' ]; then
	SnpSift annotate \
		$dbsnp_file $name/results/Mutect2/$name.Mutect2.vcf \
		> $name/results/Mutect2/$name.Mutect2.ann1.vcf

	SnpSift annotate \
		$cosmiccoding_file $name/results/Mutect2/$name.Mutect2.ann1.vcf \
		> $name/results/Mutect2/$name.Mutect2.ann2.vcf

	SnpSift annotate \
		$cosmicnoncoding_file $name/results/Mutect2/$name.Mutect2.ann2.vcf \
		> $name/results/Mutect2/$name.Mutect2.ann3.vcf

	SnpSift annotate \
		$clinvar_file $name/results/Mutect2/$name.Mutect2.ann3.vcf \
		> $name/results/Mutect2/$name.Mutect2.ann4.vcf

	SnpSift DbNSFP \
		-db $dbnsfp_file -v -f MetaLR_pred,MetaSVM_pred,SIFT_pred,PROVEAN_pred \
		$name/results/Mutect2/$name.Mutect2.ann4.vcf \
		> $name/results/Mutect2/$name.Mutect2.ann5.vcf

	snpEff $snpeff_version -canon \
		-csvStats $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats \
		$name/results/Mutect2/$name.Mutect2.ann5.vcf \
		> $name/results/Mutect2/$name.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.Mutect2.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Mutect2/$name.Mutect2.annotated.one.vcf

	SnpSift extractFields \
		$name/results/Mutect2/$name.Mutect2.annotated.one.vcf \
		CHROM POS REF ALT "GEN[Tumor].AF" "GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \
		"GEN[Normal].AD[0]" "GEN[Normal].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
		ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF GENOME_SCREEN_SAMPLE_COUNT \
		SAMPLE_COUNT CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred dbNSFP_PROVEAN_pred \
		> $name/results/Mutect2/$name.Mutect2.txt

elif [ $species = 'Mouse' ]; then
	snpEff $snpeff_version -canon \
		-csvStats $name/results/Mutect2/$name.Mutect2.annotated.vcf.stats \
		$name/results/Mutect2/$name.Mutect2.vcf \
		> $name/results/Mutect2/$name.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.Mutect2.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Mutect2/$name.Mutect2.annotated.one.vcf

	SnpSift extractFields \
		$name/results/Mutect2/$name.Mutect2.annotated.one.vcf \
		CHROM POS REF ALT "GEN[Tumor].AF" "GEN[Tumor].AD[0]" "GEN[Tumor].AD[1]" \
		"GEN[Normal].AD[0]" "GEN[Normal].AD[1]" ANN[*].GENE  ANN[*].EFFECT \
		ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		> $name/results/Mutect2/$name.Mutect2.txt
fi

gatk --java-options '-Xmx256G' IndexFeatureFile \
	-F $name/results/Mutect2/"$name".m2.filt.vcf

gatk --java-options '-Xmx256G' IndexFeatureFile \
	-F $name/results/Mutect2/"$name".Mutect2.vcf


if [ $species = 'Mouse' ]; then
	java -Xmx256G -jar $discvrseq_dir/DISCVRSeq-1.07.jar VariantQC \
		-R $genome_file \
		-V $name/results/Mutect2/"$name".m2.filt.vcf \
		-O $name/results/Mutect2/"$name".m2.filt.vcf.html \
		-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L X -L Y

	java -Xmx256G -jar $discvrseq_dir/DISCVRSeq-1.07.jar VariantQC \
		-R $genome_file \
		-V $name/results/Mutect2/"$name".Mutect2.vcf \
		-O $name/results/Mutect2/"$name".Mutect2.vcf.html \
		-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L X -L Y
elif [ $species = 'Human' ]; then
	java -Xmx256G -jar $discvrseq_dir/DISCVRSeq-1.07.jar VariantQC \
		-R $genome_file \
		-V $name/results/Mutect2/"$name".m2.filt.vcf \
		-O $name/results/Mutect2/"$name".m2.filt.vcf.html \
		-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y

	java -Xmx256G -jar $discvrseq_dir/DISCVRSeq-1.07.jar VariantQC \
		-R $genome_file \
		-V $name/results/Mutect2/"$name".Mutect2.vcf \
		-O $name/results/Mutect2/"$name".Mutect2.vcf.html \
		-L 1 -L 2 -L 3 -L 4 -L 5 -L 6 -L 7 -L 8 -L 9 -L 10 -L 11 -L 12 -L 13 -L 14 -L 15 -L 16 -L 17 -L 18 -L 19 -L 20 -L 21 -L 22 -L X -L Y
fi

sh $repository_dir/SNV_CleanUp.sh $name Mutect2 MS
