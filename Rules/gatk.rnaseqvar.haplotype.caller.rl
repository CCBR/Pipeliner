rule gatk_rnaseqvar_haplotype_caller:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:rnahapcall"
     threads: 2
     shell: "module load GATK/3.6; GATK -m 24g HaplotypeCaller -R {params.genome} -I {input} --read_filter BadCigar -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]}"


