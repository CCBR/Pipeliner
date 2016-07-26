rule gatk_merge_mutect2_chrom:
    input: vcf_chrom=expand("mutect2_out/chrom_files/{{x}}_{chr}.vcf", x=pairs, chr=config['references'][pfamily]['CHROMS'])
    output: vcf="mutect2_out/{x}.vcf",
            vcfRename="mutect2_out/{x}.FINAL.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],rname="pl:merge.mutect2"
    run: 
       F=open(input[0],"r")
       fl=F.read().split()
       F.close()
       var=" --variant "+" --variant ".join(fl)
       cmd="{params.gatk} -T CombineVariants -R {params.genome} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}"+var + "; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} -L {params.targets} --excludeFiltered -o {output.vcfRename}; sed -i 's/NORMAL/{params.normalsample}/g' {output.vcfRename}; sed -i 's/TUMOR/{params.tumorsample}/g' {output.vcfRename}"
       shell(cmd)