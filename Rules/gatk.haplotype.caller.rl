rule gatk_haplotype_caller:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".realign.bam"
     output:"{x}.g.vcf"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 8     
     shell: "{params.gatk} -T HaplotypeCaller -R {params.genome} -I {input} --read_filter BadCigar --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]}"


