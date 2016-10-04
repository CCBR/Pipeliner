rule somatic_germline_calls:
     input:  bams=expand("{s}"+".recal.bam",s=samples),
             targets="exome_targets.bed"
     output: filtered="germline_snps.vcf"
     params: genome=config['references'][pfamily]['GENOME'],regions="exome_targets.bed",knowns=config['references'][pfamily]['KNOWNVCF'],snpsites=config['references'][pfamily]['SNPSITES'],gatk=config['bin'][pfamily]['GATK'],rname="pl:germcalls"
     threads: 32
     shell:  "ls *recal.bam > inbams.list; {params.gatk} -T UnifiedGenotyper -R {params.genome} -I inbams.list -L {params.knowns} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --dbsnp {params.snpsites} -nt {threads} -o {output.filtered}"