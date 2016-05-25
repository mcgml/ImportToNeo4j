#Description: Script for annotating VCF with metadata
#Author: Matthew Lyon
#Status: RELEASE
#Mode: BY_SAMPLE/BY_COHORT
#Date: 25/05/2016
#Version: 1

#add metadata to VCF
resultsFolder="http://10.59.210.245/results"
tissue="Blood"
assay="AgilentHaloPlex"
pipelineName="AgilentHaloPlex"
pipelineVersion="1"

#get file basename
filename=$(echo "$1" | cut -d. -f1)

#print VCF headers to new file
grep '^##' "$1" > "$filename"_meta.vcf

#loop over variable files
for i in $(ls *.variables); do

	#load variables into scope
	. "$i"

	#add metadata
	echo \#\#SAMPLE\=\<ID\="$SampleID",Tissue\="$tissue",WorklistId\="$ExperimentName",SeqId\="$RunID",Assay\="$assay",PipelineName\="$pipelineName",PipelineVersion\="$pipelineVersion",RemoteBamFilePath\="$resultsFolder"/"$RunID"/"$SampleID"/"$RunID"_"$SampleID".bam,RemoteVcfFilePath="$resultsFolder"/"$RunID"/"$RunID"_Variants_Filtered.vcf\> >> "$filename"_meta.vcf

done

#add variant calls
grep -v '^##' "$1" >> "$filename"_meta.vcf

#validate and index VCF
java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.4-46/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$filename"_meta.vcf \
-dt NONE
