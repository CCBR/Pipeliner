rule snpeff_somatic:
        input: "mutect2/merged_somatic.vcf"
        output: "mutect2/merged_somatic_snpEff.vcf"
        params: snpeff=config['bin']['SNPEFF'],regions=config['references']['REFFLAT'],genome=config['references']['SNPEFF_GENOME'],snpsift=config['bin']['SNPSIFT'],rname="pl:snpeff_somatic"
        shell: "module load snpEff;{params.snpeff} -v {params.genome} -interval {params.regions} -cancer -cancerSamples pairs {input} > {output}"