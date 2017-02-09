if config['project']['annotation'] == "hg19":
  rule gatk_merge_mutect2_chrom:
    input: vcf1=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_1.vcf",
           vcf2=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_2.vcf",
           vcf3=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_3.vcf",
           vcf4=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_4.vcf",
           vcf5=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_5.vcf",
           vcf6=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_6.vcf",
           vcf7=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_7.vcf",
           vcf8=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_8.vcf",
           vcf9=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_9.vcf" ,     
           vcf10=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_10.vcf",
           vcf11=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_11.vcf",
           vcf12=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_12.vcf",
           vcf13=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_13.vcf",
           vcf14=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_14.vcf",
           vcf15=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_15.vcf",
           vcf16=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_16.vcf",
           vcf17=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_17.vcf",
           vcf18=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_18.vcf",
           vcf19=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_19.vcf",
           vcf20=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_20.vcf",
           vcf21=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_21.vcf",
           vcf22=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_22.vcf",
           vcfX=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_X.vcf",
           vcfY=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_Y.vcf"
    output: vcf=config['project']['workpath']+"/mutect2_out/{x}_mutect2.vcf",
            snps=config['project']['workpath']+"/mutect2_out/{x}.SNPs.vcf",
            indels=config['project']['workpath']+"/mutect2_out/{x}.INDELs.vcf",
            flagsnps=config['project']['workpath']+"/mutect2_out/{x}.flaggedSNPs.vcf",
            flagindels=config['project']['workpath']+"/mutect2_out/{x}.flaggedINDELs.vcf",
            flagged=config['project']['workpath']+"/mutect2_out/{x}.flagged.vcf",
            filtered=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],rname="pl:merge.mutect2"
    shell: "module load GATK/3.6; GATK -m 64G CombineVariants -R {params.genome} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf} --variant {input.vcf1} --variant {input.vcf2} --variant {input.vcf3} --variant {input.vcf4} --variant {input.vcf5} --variant {input.vcf6} --variant {input.vcf7} --variant {input.vcf8} --variant {input.vcf9} --variant {input.vcf10} --variant {input.vcf11} --variant {input.vcf12} --variant {input.vcf13} --variant {input.vcf14} --variant {input.vcf15} --variant {input.vcf16} --variant {input.vcf17} --variant {input.vcf18} --variant {input.vcf19} --variant {input.vcf20} --variant {input.vcf21} --variant {input.vcf22} --variant {input.vcfX} --variant {input.vcfY}; GATK -m 120G SelectVariants -R {params.genome} --variant {output.vcf} -selectType SNP --excludeFiltered -o {output.snps}; GATK -m 120G SelectVariants -R {params.genome} --variant {output.vcf} -selectType INDEL --excludeFiltered -o {output.indels}; GATK -m 48G VariantFiltration -R {params.genome} --variant {output.snps} --filterExpression \"FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\" -o {output.flagsnps}; GATK -m 48G VariantFiltration -R {params.genome} --variant {output.indels} --filterExpression \"FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o {output.flagindels}; GATK -m 48G CombineVariants -R {params.genome} --variant {output.flagsnps} --variant {output.flagindels} -o {output.flagged} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNSORTED; GATK -m 48G SelectVariants -R {params.genome} --variant {output.flagged} --excludeFiltered -o {output.filtered}"


else:
  rule gatk_merge_mutect2_chrom:
    input: vcfchr1=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr1.vcf",
            vcfchr2=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr2.vcf",
            vcfchr3=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr3.vcf",
            vcfchr4=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr4.vcf",
            vcfchr5=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr5.vcf",
            vcfchr6=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr6.vcf",
            vcfchr7=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr7.vcf",
            vcfchr8=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr8.vcf",
            vcfchr9=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr9.vcf",
            vcfchr10=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr10.vcf",
            vcfchr11=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr11.vcf",
            vcfchr12=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr12.vcf",
            vcfchr13=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr13.vcf",
            vcfchr14=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr14.vcf",
            vcfchr15=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr15.vcf",
            vcfchr16=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr16.vcf",
            vcfchr17=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr17.vcf",
            vcfchr18=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr18.vcf",
            vcfchr19=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr19.vcf",
            vcfchrX=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chrX.vcf",
            vcfchrY=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chrY.vcf"
    output: vcf=config['project']['workpath']+"/mutect2_out/{x}.vcf",
            snps=config['project']['workpath']+"/mutect2_out/{x}.SNPs.vcf",
            indels=config['project']['workpath']+"/mutect2_out/{x}.INDELs.vcf",
            flagsnps=config['project']['workpath']+"/mutect2_out/{x}.flaggedSNPs.vcf",
            flagindels=config['project']['workpath']+"/mutect2_out/{x}.flaggedINDELs.vcf",
            flagged=config['project']['workpath']+"/mutect2_out/{x}.flagged.vcf",
            filtered=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],rname="pl:merge.mutect2"
    shell: "module load GATK/3.6; GATK -m 64G CombineVariants -R {params.genome} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf} --variant {input.vcfchr1} --variant {input.vcfchr2} --variant {input.vcfchr3} --variant {input.vcfchr4} --variant {input.vcfchr5} --variant {input.vcfchr6} --variant {input.vcfchr7} --variant {input.vcfchr8} --variant {input.vcfchr9} --variant {input.vcfchr10} --variant {input.vcfchr11} --variant {input.vcfchr12} --variant {input.vcfchr13} --variant {input.vcfchr14} --variant {input.vcfchr15} --variant {input.vcfchr16} --variant {input.vcfchr17} --variant {input.vcfchr18} --variant {input.vcfchr19} --variant {input.vcfchrX} --variant {input.vcfchrY}; GATK -m 120G SelectVariants -R {params.genome} --variant {output.vcf} -selectType SNP --excludeFiltered -o {output.snps}; GATK -m 120G SelectVariants -R {params.genome} --variant {output.vcf} -selectType INDEL --excludeFiltered -o {output.indels}; GATK -m 48G VariantFiltration -R {params.genome} --variant {output.snps} --filterExpression \"FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\" -o {output.flagsnps}; GATK -m 48G VariantFiltration -R {params.genome} --variant {output.indels} --filterExpression \"FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o {output.flagindels}; GATK -m 48G CombineVariants -R {params.genome} --variant {output.flagsnps} --variant {output.flagindels} -o {output.flagged} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNSORTED; GATK -m 48G SelectVariants -R {params.genome} --variant {output.flagged} --excludeFiltered -o {output.filtered}"