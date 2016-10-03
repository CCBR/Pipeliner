rule snpeff:
    input: vcf="exome.recode.vcf",
           strict="exome.strictFilter.vcf",
           relax="exome.relaxedFilter.vcf",
           targets="exome_targets.bed"
    output: vcf="exome.snpeff.vcf",
            strict="exome.strictFilter.snpeff.vcf",
            relax="exome.relaxedFilter.snpeff.vcf"
    params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:snpeff"
    run: shell("{params.snpeff} -v -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input.vcf} > {output.vcf}; {params.snpeff} -v -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input.strict} > {output.strict}; {params.snpeff} -v -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input.relax} > {output.relax}")