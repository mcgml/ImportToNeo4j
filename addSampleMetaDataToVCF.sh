#add metadata to VCF
resultsFolder="http://10.59.210.245/results"

#get file basename
filename=$(echo "$1" | cut -d. -f1)

#print VCF headers to new file
grep '^##' "$1" > "$filename"_meta.vcf

#loop over variable files
for i in $(ls *.variables); do

	#load variables into scope
	. "$i"

	#add metadata
	echo \#\#SAMPLE\=\<ID\="$SampleID",Tissue\=Blood,WorklistId\="$ExperimentName",SeqId\="$RunID",Assay\=TruSightOne,PipelineName\=IlluminaTruSightOne,PipelineVersion\=1,RemoteBamFilePath\="$resultsFolder"/"$RunID"/"$SampleID"/"$RunID"_"$SampleID".bam,RemoteVcfFilePath="$resultsFolder"/"$RunID"/"$RunID"_Variants_Filtered.vcf\> >> "$filename"_meta.vcf

done

#add variant calls
grep -v '^##' "$1" >> "$filename"_meta.vcf

#validate and index VCF
~/jre1.8.0_71/bin/java -Xmx2g -jar /share/apps/GATK-distros/GATK_3.4-46/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$filename"_meta.vcf \
-dt NONE
