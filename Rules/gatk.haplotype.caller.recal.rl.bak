rule gatk_haplotype_caller_recal:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: gres="lscratch:100",gatk=config['bin']['GATK'],genome=config['references']['GENOME'],snpsites=config['references']['SNPSITES'],rname="pl:hapcall"
     threads: 8     
     shell: "{params.gatk} -T HaplotypeCaller -R {params.genome} -I {input} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]}"


