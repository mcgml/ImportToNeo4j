#!/bin/bash
set -euo pipefail

#Description: Script for importing variants into variantDB
#Author: Matthew Lyon
#Status: Development
#Mode: BY_SAMPLE/BY_COHORT
#Date: 25/05/2016
#Version: 1.1

PATH=$PATH:/share/apps/bigWigToWig-distros


#check args
if [ "$#" -ne 1 ]; then
    echo "Usage <VCF>"
fi

#check index is present
if [ ! -f "$1".idx ]; then
	echo VCF not indexed
	exit
fi

#import variants and genotypes into DB
echo importing variants to DB
~/jre1.8.0_71/bin/java -Xmx16g -jar ./import2neo4j/ImportToNeo4j-1.0.4.jar \
"$1" \
graph.db


#sort imported variants VCF
echo sorting imported variants
java -jar /share/apps/picard-tools-distros/picard-tools-1.129/picard.jar SortVcf \
I=imported.vcf \
O=imported.sorted.vcf \
CREATE_INDEX=false \
SD=/data/db/human/gatk/2.8/b37/human_g1k_v37.dict


#annotate variants
echo annotating imported variants
perl ./ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-v \
-i imported.sorted.vcf \
--format vcf \
-o imported.sorted.vep.vcf \
--vcf \
--force_overwrite \
--species homo_sapiens \
--assembly GRCh37 \
--everything \
--cache \
--dir /data/diagnostics/apps/VariantDatabase/ensembl-tools-release-82/scripts/variant_effect_predictor/vep \
--offline \
--shift_hgvs 1 \
--fork 8 \
--cache_version 82 \
-custom /data/db/human/GERP/All_hg19_RS.bw,GERP,bigwig \
-custom /data/db/human/phyloP/hg19.100way.phyloP100way.bw,phyloP,bigwig \
-custom /data/db/human/phastCons/hg19.100way.phastCons.bw,phastCons,bigwig \
--no_stats


if [ -f imported.sorted.vep.vcf ]; then

	echo annotations found

	#add  dbSNPId
	/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
	-T VariantAnnotator \
	-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
	-V imported.sorted.vep.vcf \
	-o imported.sorted.vep.dbsnp.vcf \
	--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
	-L imported.sorted.vep.vcf \
	-ip 300 \
	--resourceAlleleConcordance \
	-dt NONE
        
        # bgzip and tabix
        bgzip imported.sorted.vep.dbsnp.vcf
        tabix imported.sorted.vep.dbsnp.vcf.gz
        
        # add correct frequencies using vcfanno
        /share/apps/vcfanno/vcfanno_linux64 af_anno.toml imported.sorted.vep.dbsnp.vcf.gz  > imported.sorted.vep.af.vcf
	
	# Remove format field which we don't need and causes an error
	awk  ' BEGIN { OFS = "\t" } { print $1,$2,$3,$4,$5,$6,$7,$8 }' imported.sorted.vep.af.vcf > imported.sorted.vep.af.fixed.vcf

        #index vcf
	/share/apps/igvtools-distros/igvtools_2.3.75/igvtools index imported.sorted.vep.af.fixed.vcf	


	#import annotations
	~/jre1.8.0_71/bin/java -Xmx16g -jar ./import2neo4j/ImportToNeo4j-1.0.4.jar \
	imported.sorted.vep.af.fixed.vcf \
	graph.db \
	-a	



	#clean up
	rm imported.sorted.vep.vcf
	rm imported.sorted.vep.vcf.idx
	rm imported.sorted.vep.dbsnp.vcf.gz
	rm imported.sorted.vep.dbsnp.vcf.gz.tbi
        rm imported.sorted.vep.af.fixed.vcf
        rm imported.sorted.vep.af.fixed.vcf.idx
	rm imported.sorted.vep.af.vcf
	rm imported.sorted.vep.dbsnp.vcf.idx
else
	echo annotations not found

fi

#clean up
rm imported.vcf
rm imported.sorted.vcf
