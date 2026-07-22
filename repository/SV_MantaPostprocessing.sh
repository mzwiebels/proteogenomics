#!/bin/bash

##########################################################################################
##
## SNV_MantaPostprocessing.sh
##
## Postprocessing for Manta.
##
##########################################################################################

name=$1
config_file=$2
runmode=$3
type=$4

. $config_file

if [ $runmode = 'MS' ]; then
	cp $name/results/Manta/results/variants/somaticSV.vcf.gz $name/results/Manta/$name.man.vcf.gz
	gunzip $name/results/Manta/$name.man.vcf.gz
	cat $name/results/Manta/$name.man.vcf | SnpSift filter "( ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].SR[1] + GEN[Tumor].SR[0] ) * 0.05 <= GEN[Tumor].SR[1] ) ) | ( ( FILTER = 'PASS' ) & ( ( GEN[Tumor].PR[1] + GEN[Tumor].PR[0] ) * 0.05 <= GEN[Tumor].PR[1] ) ) )" \
		> $name/results/Manta/$name.Manta.vcf

	bgzip $name/results/Manta/$name.Manta.vcf
	tabix -p vcf $name/results/Manta/$name.Manta.vcf.gz

	bcftools stats $name/results/Manta/$name.Manta.vcf.gz > $name/results/Manta/$name.Manta.vcf.gz.stats
	gunzip $name/results/Manta/$name.Manta.vcf.gz

	snpEff $snpeff_version -canon \
		-csvStats $name/results/Manta/$name.Manta.annotated.vcf.stats \
		$name/results/Manta/$name.Manta.vcf
		> $name/results/Manta/$name.Manta.annotated.vcf

	cat $name/results/Manta/$name.Manta.annotated.vcf | vcfEffOnePerLine.pl > $name/results/Manta/$name.Manta.annotated.one.vcf
	SnpSift extractFields $name/results/Manta/$name.Manta.annotated.one.vcf CHROM POS REF ALT "GEN[Tumor].SR[0]" "GEN[Tumor].SR[1]" "GEN[Tumor].PR[0]" "GEN[Tumor].PR[1]" "GEN[Normal].SR[0]" "GEN[Normal].SR[1]" "GEN[Normal].PR[0]" "GEN[Normal].PR[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		> $name/results/Manta/$name.Manta.txt

elif [ $runmode = 'SS' ]; then
	cp $name/results/Manta/results/variants/diploidSV.vcf.gz $name/results/Manta/$name.man.vcf.gz
	gunzip $name/results/Manta/$name.man.vcf.gz
	cat $name/results/Manta/$name.man.vcf | SnpSift filter "( ( FILTER = 'PASS' ) )" > $name/results/Manta/$name.$type.Manta.vcf

	bgzip $name/results/Manta/$name.$type.Manta.vcf
	tabix -p vcf $name/results/Manta/$name.$type.Manta.vcf.gz

	bcftools stats $name/results/Manta/$name.$type.Manta.vcf.gz > $name/results/Manta/$name.$type.Manta.vcf.gz.stats
	gunzip $name/results/Manta/$name.$type.Manta.vcf.gz

	snpEff $snpeff_version -canon \
		-csvStats $name/results/Manta/$name.Manta.annotated.vcf.stats \
		$name/results/Manta/$name.$type.Manta.vcf \
		> $name/results/Manta/$name.$type.Manta.annotated.vcf

	cat $name/results/Manta/$name.$type.Manta.annotated.vcf | vcfEffOnePerLine.pl > $name/results/Manta/$name.$type.Manta.annotated.one.vcf
	SnpSift extractFields $name/results/Manta/$name.$type.Manta.annotated.one.vcf CHROM POS REF ALT "GEN[$type].SR[0]" "GEN[$type].SR[1]" "GEN[$type].PR[0]" "GEN[$type].PR[1]" ANN[*].GENE  ANN[*].EFFECT ANN[*].IMPACT ANN[*].FEATUREID ANN[*].HGVS_C ANN[*].HGVS_P \
		> $name/results/Manta/$name.$type.Manta.txt
fi

sh $repository_dir/SNV_CleanUp.sh $name Manta