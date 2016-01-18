#add metadata to VCF

#get file basename
filename=$(echo "$1" | cut -d. -f1)

#print VCF headers to new file
grep '^##' "$1" > "$filename"_meta.vcf

#loop over variable files
for i in $(ls *.variables); do

	#load variables into scope
	. "$i"

	#add metadata
	echo "##SAMPLE=<ID=$sampleId,Tissue=$tissue,WorklistId=$worklistId,SeqId=$seqId,Assay=$assay,PipelineName=$pipelineName,PipelineVersion=$pipelineVersion,RemoteBamFilePath=$remoteBamFilePath,RemoteVcfFilePath=$remoteVcfFilePath>" >> "$filename"_meta.vcf

done

#add variant calls
grep -v '^##' "$1" >> "$filename"_meta.vcf

#validate and index VCF
java -jar /share/apps/gatk-distros/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$filename"_meta.vcf \
-dt NONE
