#import variants and genotypes into DB
java -Xmx16g -jar ./ImportToNeo4j.jar \
"$1" \
data

#sort imported variants VCF
/share/apps/igvtools-distros/igvtools_2.3.57/igvtools sort imported.vcf imported.sorted.vcf

#annotate variants
perl ~/ensembl-tools-release-82/scripts/variant_effect_predictor/variant_effect_predictor.pl \
-v \
-i imported.sorted.vcf \
--format vcf \
-o imported.sorted.vep.vcf \
--vcf \
--force_overwrite \
--species homo_sapiens \
--assembly GRCh37 \
--everything \
--fork 12 \
--cache \
--offline \
--shift_hgvs 1 \
--cache_version 82 \
-custom /data/db/human/GERP/All_hg19_RS.bw,GERP,bigwig \
-custom /data/db/human/phyloP/hg19.100way.phyloP100way.bw,phyloP,bigwig \
-custom /data/db/human/phastCons/hg19.100way.phastCons.bw,phastCons,bigwig \
--no_stats

#add population frequencies and dbSNPId
java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.4-46/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T VariantAnnotator \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V imported.sorted.vep.vcf \
-o imported.sorted.vep.af.vcf \
--resource:kGPhase3 /data/db/human/1kg/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf \
--resource:exac /data/db/human/ExAC/ExAC.r0.3.sites.vep.vcf \
-E kGPhase3.EAS_AF \
-E kGPhase3.EUR_AF \
-E kGPhase3.AFR_AF \
-E kGPhase3.AMR_AF \
-E kGPhase3.SAS_AF \
-E exac.AC_AFR \
-E exac.AC_AMR \
-E exac.AC_EAS \
-E exac.AC_FIN \
-E exac.AC_NFE \
-E exac.AC_OTH \
-E exac.AC_SAS \
-E exac.AN_AFR \
-E exac.AN_AMR \
-E exac.AN_EAS \
-E exac.AN_FIN \
-E exac.AN_NFE \
-E exac.AN_OTH \
-E exac.AN_SAS \
-L imported.sorted.vep.vcf \
-ip 300 \
-dt NONE

#import annotations
java -Xmx16g -jar ./ImportToNeo4j.jar \
imported.sorted.vep.af.vcf \
data \
-a
