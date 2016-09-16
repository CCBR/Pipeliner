rule snpeff_genomeseq:
    input: "combined.vcf"
    output: "combined.snpeff.vcf"
    params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:snpeff_genome"
    run: shell("{params.snpeff} -v -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input} > {output}")