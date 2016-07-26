rule snpeff:
    input: "exome.recode.vcf"
    output: "exome.snpeff.vcf"
    params: snpeff=config['bin']['SNPEFF'],genome=config['references']['SNPEFF_GENOME'],snpsift=config['bin']['SNPSIFT'],effconfig=config['references']['SNPEFF_CONFIG'],rname="pl:snpeff"
    run: "module load snpEff; {params.snpeff} -v {params.genome} -c {params.effconfig} -stats snpeff.stats.html {input} > exome.snpeff.vcf"