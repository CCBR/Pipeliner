rule gatk_merge_somatic_vcfs:
    input: expand("mutect2/{p}"+".FINAL.vcf",p=pairs)
    output: vcf="mutect2/merged_somatic.vcf",
#            summary="mutect2/somatic_variants"
    params: regions=config['references']['REFFLAT'],gres="lscratch:100",gatk=config['bin']['GATK'],genome=config['references']['GENOME'],snpsites=config['references']['SNPSITES'],rname="CombineVariants"
    run:
        fl=os.popen("ls mutect2/*.FINAL.vcf").read().split()
        var=" --variant "+" --variant ".join(fl)
        cmd="{params.gatk} -T CombineVariants -R {params.genome} -L {paras.regions} --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}"+var
        system(cmd)
#        cmd="module load picard/2.1.1; java -Xmx24g -jar $PICARDJARPATH/picard.jar CollectVariantCallingMetrics INPUT={output.vcf} OUTPUT={output.summary}"
#        system(cmd)