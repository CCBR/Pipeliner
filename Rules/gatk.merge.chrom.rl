if config['project']['annotation'] == "hg19":
  rule gatk_merge_chrom:
    input: vcf1="combined_1.vcf",
           vcf2="combined_2.vcf",
           vcf3="combined_3.vcf",
           vcf4="combined_4.vcf",
           vcf5="combined_5.vcf",
           vcf6="combined_6.vcf",
           vcf7="combined_7.vcf",
           vcf8="combined_8.vcf",
           vcf9="combined_9.vcf",     
           vcf10="combined_10.vcf",
           vcf11="combined_11.vcf",
           vcf12="combined_12.vcf",
           vcf13="combined_13.vcf",
           vcf14="combined_14.vcf",
           vcf15="combined_15.vcf",
           vcf16="combined_16.vcf",
           vcf17="combined_17.vcf",
           vcf18="combined_18.vcf",
           vcf19="combined_19.vcf",
           vcf20="combined_20.vcf",
           vcf21="combined_21.vcf",
           vcf22="combined_22.vcf",
           vcfX="combined_X.vcf",
           vcfY="combined_Y.vcf",
           vcfMT="combined_MT.vcf",
           targets=config['project']['workpath']+"/exome_targets.bed"
    output: vcf="combined.vcf",
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],dir=config['project']['workpath'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:merge.mutect2",
    threads: 4
    shell: "module load GATK/3.7; GATK -m 64G CombineVariants -R {params.genome} --filteredrecordsmergetype KEEP_UNCONDITIONAL --assumeIdenticalSamples -o {output.vcf} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} --variant {input.vcf4} --variant {input.vcf5} --variant {input.vcf6} --variant {input.vcf7} --variant {input.vcf8} --variant {input.vcf9} --variant {input.vcf10} --variant {input.vcf11} --variant {input.vcf12} --variant {input.vcf13} --variant {input.vcf14} --variant {input.vcf15} --variant {input.vcf16} --variant {input.vcf17} --variant {input.vcf18} --variant {input.vcf19} --variant {input.vcf20} --variant {input.vcf21} --variant {input.vcf22} --variant {input.vcfX} --variant {input.vcfY} --variant {input.vcfMT} -nt {threads}"

elif config['project']['annotation'] == "hg38":
  rule gatk_merge_chrom:
    input: vcf1="combined_1.vcf",
           vcf2="combined_2.vcf",
           vcf3="combined_3.vcf",
           vcf4="combined_4.vcf",
           vcf5="combined_5.vcf",
           vcf6="combined_6.vcf",
           vcf7="combined_7.vcf",
           vcf8="combined_8.vcf",
           vcf9="combined_9.vcf",     
           vcf10="combined_10.vcf",
           vcf11="combined_11.vcf",
           vcf12="combined_12.vcf",
           vcf13="combined_13.vcf",
           vcf14="combined_14.vcf",
           vcf15="combined_15.vcf",
           vcf16="combined_16.vcf",
           vcf17="combined_17.vcf",
           vcf18="combined_18.vcf",
           vcf19="combined_19.vcf",
           vcf20="combined_20.vcf",
           vcf21="combined_21.vcf",
           vcf22="combined_22.vcf",
           vcfX="combined_X.vcf",
           vcfY="combined_Y.vcf",
           targets=config['project']['workpath']+"/exome_targets.bed"
    output: vcf="combined.vcf",
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],dir=config['project']['workpath'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:merge.mutect2",
    threads: 4
    shell: "module load GATK/3.7; GATK -m 64G CombineVariants -R {params.genome} --filteredrecordsmergetype KEEP_UNCONDITIONAL --assumeIdenticalSamples -o {output.vcf} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} --variant {input.vcf4} --variant {input.vcf5} --variant {input.vcf6} --variant {input.vcf7} --variant {input.vcf8} --variant {input.vcf9} --variant {input.vcf10} --variant {input.vcf11} --variant {input.vcf12} --variant {input.vcf13} --variant {input.vcf14} --variant {input.vcf15} --variant {input.vcf16} --variant {input.vcf17} --variant {input.vcf18} --variant {input.vcf19} --variant {input.vcf20} --variant {input.vcf21} --variant {input.vcf22} --variant {input.vcfX} --variant {input.vcfY} --variant {input.vcfMT} -nt {threads}"


elif config['project']['annotation'] == "mm10":
  rule gatk_merge_chrom:
    input: vcf1="combined_1.vcf",
            vcf2="combined_2.vcf",
            vcf3="combined_3.vcf",
            vcf4="combined_4.vcf",
            vcf5="combined_5.vcf",
            vcf6="combined_6.vcf",
            vcf7="combined_7.vcf",
            vcf8="combined_8.vcf",
            vcf9="combined_9.vcf",
            vcf10="combined_10.vcf",
            vcf11="combined_11.vcf",
            vcf12="combined_12.vcf",
            vcf13="combined_13.vcf",
            vcf14="combined_14.vcf",
            vcf15="combined_15.vcf",
            vcf16="combined_16.vcf",
            vcf17="combined_17.vcf",
            vcf18="combined_18.vcf",
            vcf19="combined_19.vcf",
            vcfX="combined_X.vcf",
            vcfY="combined_Y.vcf",
            targets=config['project']['workpath']+"/exome_targets.bed"
    output: vcf="combined.vcf",
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],dir=config['project']['workpath'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:merge.mutect2",
    threads: 4
    shell: "module load GATK/3.7; GATK -m 64G CombineVariants -R {params.genome} --filteredrecordsmergetype KEEP_UNCONDITIONAL --assumeIdenticalSamples -o {output.vcf} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} --variant {input.vcf4} --variant {input.vcf5} --variant {input.vcf6} --variant {input.vcf7} --variant {input.vcf8} --variant {input.vcf9} --variant {input.vcf10} --variant {input.vcf11} --variant {input.vcf12} --variant {input.vcf13} --variant {input.vcf14} --variant {input.vcf15} --variant {input.vcf16} --variant {input.vcf17} --variant {input.vcf18} --variant {input.vcf19} --variant {input.vcfX} --variant {input.vcfY}  --variant {input.vcfMT} -nt {threads}"