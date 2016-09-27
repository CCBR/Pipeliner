rule gatk_merge_somatic_vcfs:
    input: mutect=expand("mutect_out/{p}"+".FINAL.vcf",p=pairs),
           strelka=expand("strelka_out/{p}"+".vcf",p=pairs)
    output: mutectvcf="mutect_out/merged_somatic.vcf",
            strelkavcf="strelka_out/merged_somatic.vcf"
    params: regions="exome_targets.bed",gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="CombineVariants"
    run:
        fl=os.popen("ls mutect_out/*.FINAL.vcf").read().split()
        var=" --variant "+" --variant ".join(fl)
        cmd="{params.gatk} -T CombineVariants -R {params.genome} -L {params.regions} --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNIQUIFY -o {output.mutectvcf}"+var
        shell(cmd)
        flb=os.popen("ls strelka_out/*.vcf | grep -v 'snpEff'").read().split()
        varb=" --variant "+" --variant ".join(flb)
        cmd="{params.gatk} -T CombineVariants -R {params.genome} -L {params.regions} --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNIQUIFY -o {output.strelkavcf}"+varb
        shell(cmd)        