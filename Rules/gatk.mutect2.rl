rule gatk_mutect2:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: vcf="mutect2_out/{x}.paired.vcf",
#            vcfFilter="mutect2_out/{x}.PASS.vcf",
            vcfRename="mutect2_out/{x}.FINAL.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2"
    threads: 2
    shell: "{params.gatk} -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L {params.targets} -nct {threads} -o {output.vcf}; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcfRename}; sed -i 's/NORMAL/{params.normalsample}/g' {output.vcfRename}; sed -i 's/TUMOR/{params.tumorsample}/g' {output.vcfRename}"