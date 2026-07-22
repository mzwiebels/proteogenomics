#!/bin/bash

##########################################################################################
##
## SNV_Mutect2Postprocessing.sh
##
## Postprocessing for Mutect2 in single-sample mode.
##
##########################################################################################

name=$1
species=$2
config_file=$3
type=$4
filtering=$5
artefact_type=$6

. $config_file

echo '---- Mutect2 SS Postprocessing I (OrientationFilter, Indel size selection, filtering) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

gatk --java-options '-Xmx256G' FilterMutectCalls \
	--variant $name/results/Mutect2/$name.$type.m2.vcf \
	--output $name/results/Mutect2/$name.$type.m2.filt.vcf \
	--reference $genome_file

SnpSift extractFields \
	$name/results/Mutect2/"$name".$type.m2.filt.vcf \
	CHROM POS REF ALT "GEN["$type"].AF" "GEN["$type"].AD[0]" \
	"GEN["$type"].AD[1]" MMQ[1] MBQ[1] \
	> $name/results/Mutect2/$name.$type.Mutect2.Positions.txt

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.filt.vcf $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
elif [ $artefact_type = 'CT' ]; then
	gatk --java-options '-Xmx256G' FilterByOrientationBias \
		-V $name/results/Mutect2/$name.$type.m2.filt.vcf \
		-P $name/results/QC/$name.$type.bam.artifacts.pre_adapter_detail_metrics \
		--artifact-modes C/T --output $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
elif [ $artefact_type = 'GT' ]; then
	gatk --java-options '-Xmx256G' FilterByOrientationBias \
		-V $name/results/Mutect2/$name.$type.m2.filt.vcf \
		-P $name/results/QC/$name.$type.bam.artifacts.pre_adapter_detail_metrics \
		--artifact-modes G/T --output $name/results/Mutect2/$name.$type.m2.filt.AM.vcf
fi

if [ $artefact_type = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.filt.AM.vcf $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
elif [ $artefact_type = 'CT' ] || [ $artefact_type = 'GT' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.AM.vcf | \
	SnpSift filter \
		"( ( ( FILTER = 'PASS'  ) & (exists GEN[$type].OBP) & \
		(GEN[$type].OBP <= 0.05) ) | ( ( FILTER = 'PASS' ) ) )" \
		> $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf
fi

gatk --java-options '-Xmx256G' SelectVariants --max-indel-size 10 \
	-V $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.vcf \
	-output $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf

if [ $filtering = 'all' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS') & (GEN[$type].AF >= 0.05) & \
		(GEN[$type].AD[1] >= 2) & (GEN[$type].AD[0] + GEN[$type].AD[1] >= 5) )" \
		> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
elif [ $filtering = 'hard' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS') & (GEN[$type].AF >= 0.1) & \
		(GEN[$type].AD[1] >= 2) & (GEN[$type].AD[0] + GEN[$type].AD[1] >= 10) )" \
		> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
elif [ $filtering = 'none' ]; then
	cat $name/results/Mutect2/$name.$type.m2.filt.AM.filtered.selected.vcf | \
	SnpSift filter \
		"( ( FILTER = 'PASS' ) )" \
		> $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
fi

echo '---- Mutect2 SS Postprocessing II (Filtering out known SNV/Indel using dbSNP or the Sanger Mouse database) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

bgzip -f $name/results/Mutect2/$name.$type.m2.postprocessed.vcf
tabix -p vcf $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz

if [ $filtering = 'all' ] || [ $filtering = 'hard' ]; then
	bcftools isec -C -c none -O z -w 1 \
		-o $name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz \
		$name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz \
		$snp_file
elif [ $filtering = 'none' ]; then
	cp $name/results/Mutect2/$name.$type.m2.postprocessed.vcf.gz \
	$name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz
fi

bcftools norm -m -any \
	$name/results/Mutect2/$name.$type.m2.postprocessed.snp_removed.vcf.gz \
	-O z -o $name/results/Mutect2/$name.$type.Mutect2.vcf.gz

gunzip $name/results/Mutect2/$name.$type.Mutect2.vcf.gz


echo '---- Mutect2 SS Postprocessing III (Annotate calls) ----' | tee -a $name/results/QC/$name.report.txt
echo "$(date) \t timestamp: $(date +%s)" | tee -a $name/results/QC/$name.report.txt

if [ $species = 'Human' ]; then
	SnpSift annotate \
		$dbsnp_file $name/results/Mutect2/$name.$type.Mutect2.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf

	SnpSift annotate \
		$cosmiccoding_file $name/results/Mutect2/$name.$type.Mutect2.ann1.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf

	SnpSift annotate \
		$cosmicnoncoding_file $name/results/Mutect2/$name.$type.Mutect2.ann2.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf

	SnpSift annotate \
		$clinvar_file $name/results/Mutect2/$name.$type.Mutect2.ann3.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.ann4.vcf

	SnpSift DbNSFP \
		-db $dbnsfp_file -v -f MetaLR_pred,MetaSVM_pred,SIFT_pred,PROVEAN_pred \
		$name/results/Mutect2/$name.$type.Mutect2.ann4.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.ann5.vcf

	snpEff $snpeff_version -canon \
		-csvStats $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats \
		$name/results/Mutect2/$name.$type.Mutect2.ann5.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf

	SnpSift extractFields \
		$name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf \
		CHROM POS REF ALT "GEN[$type].AF" "GEN[$type].AD[0]" \
		"GEN[$type].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
		ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		dbNSFP_MetaLR_pred dbNSFP_MetaSVM_pred ID G5 AC AN AF GENOME_SCREEN_SAMPLE_COUNT \
		SAMPLE_COUNT CLNDN CLNSIG CLNREVSTAT dbNSFP_SIFT_pred dbNSFP_PROVEAN_pred \
		> $name/results/Mutect2/$name.$type.Mutect2.txt

elif [ $species = 'Mouse' ]; then
	snpEff $snpeff_version -canon \
		-csvStats $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf.stats \
		$name/results/Mutect2/$name.$type.Mutect2.vcf \
		> $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf

	cat $name/results/Mutect2/$name.$type.Mutect2.annotated.vcf | \
		vcfEffOnePerLine.pl \
		> $name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf

	SnpSift.jar extractFields \
		$name/results/Mutect2/$name.$type.Mutect2.annotated.one.vcf \
		CHROM POS REF ALT "GEN[$type].AF" "GEN[$type].AD[0]" \
		"GEN[$type].AD[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT \
		ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		> $name/results/Mutect2/$name.$type.Mutect2.txt
fi

sh $repository_dir/SNV_CleanUp.sh $name Mutect2 SS $type