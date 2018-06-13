rule snpeff_genomeseq:
    input: vcf="combined.vcf",
           strict="combined.strictFilter.vcf",
           relax="combined.relaxedFilter.vcf"
    output: vcf="combined.snpeff.vcf",
            strict="combined.strictFilter.snpeff.vcf",
            relax="combined.relaxedFilter.snpeff.vcf"
    params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:snpeff_genome"
    run: shell("module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input.vcf} > {output.vcf}; java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff.strict.stats.csv -stats snpeff.strict.stats.html {params.genome} {input.strict} > {output.strict}; java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff.relaxed.stats.csv -stats snpeff.relaxed.stats.html {params.genome} {input.relax} > {output.relax}")