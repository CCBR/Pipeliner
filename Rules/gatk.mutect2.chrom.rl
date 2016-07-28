rule gatk_chrom_mutect2:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: vcfchrom="mutect2_out/chrom_files/{x}_{CHROMS}.vcf", CHROMS=config['references'][pfamily]['CHROMS']
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=config['references'][pfamily]['REFFLAT'],rname="pl:chrom.mutect2"
    shell: "{params.gatk} -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar {params.knowns} -L {CHROMS} -o {output.vcfchrom}"