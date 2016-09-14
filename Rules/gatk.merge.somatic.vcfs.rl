rule gatk_merge_somatic_vcfs:
    input: expand("mutect_out/{p}"+".FINAL.vcf",p=pairs)
    output: vcf="mutect_out/merged_somatic.vcf",
#            summary="mutect_out/somatic_variants"
    params: regions="exome_targets.bed",gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="CombineVariants"
    run:
        fl=os.popen("ls mutect_out/*.FINAL.vcf").read().split()
        var=" --variant "+" --variant ".join(fl)
        cmd="{params.gatk} -T CombineVariants -R {params.genome} -L {params.regions} --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNIQUIFY -o {output.vcf}"+var
        shell(cmd)