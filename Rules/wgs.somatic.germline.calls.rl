rule somatic_wgs_germline_calls:
     input:  bams=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output: "{x}.g.vcf"
     params: genome=config['references'][pfamily]['GENOME'],regions="exome_targets.bed",knowns=config['references'][pfamily]['KNOWNVCF2'],snpsites=config['references'][pfamily]['SNPSITES'],gatk=config['bin'][pfamily]['GATK'],rname="pl:germcalls"
     threads: 4
     shell:  "module load GATK/3.8-0; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bams} --emitRefConfidence GVCF --use_jdk_inflater --use_jdk_deflater --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp {params.snpsites} -nct {threads}  -o {output}"