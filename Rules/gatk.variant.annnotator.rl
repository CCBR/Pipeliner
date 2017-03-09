rule gatk_variant_annotator:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".realign.bam",
            vcf=lambda wildcards: config['project']['units'][wildcards.x]+".g.vcf"
     output: "{x}.annot.g.vcf"
     params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:varannot"
     shell: "{params.gatk} -T VariantAnnotator \
    -R {params.genome} \
    -I {input.bam} \
    -o {output} \
    -A Coverage \
    -A ChromosomeCounts \
    -A HaplotypeScore \
    -A MappingQualityRankSumTest \
    -A MappingQualityZero \
    -A RMSMappingQuality \
    -A SpanningDeletions \
    -A TandemRepeatAnnotator \
    -A FisherStrand \
    -A QualByDepth \
    -A BaseQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A StrandOddsRatio \
    --variant {input.vcf} \
    -L {input.vcf} \
    --dbsnp {params.snpsites};"
    
