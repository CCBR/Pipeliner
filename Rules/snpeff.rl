rule snpeff:
    input: "exome.recode.vcf"
    output: "exome.snpeff.vcf"
    params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:snpeff"
    run: shell("{params.snpeff} -v -c {params.effconfig} -stats snpeff.stats.html {params.genome} {input} > {output}")