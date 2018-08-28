rule gatk_haplotype_caller:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-0; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]}"


