rule germline_filter_rnaseqvar:
    input: "combined.vcf"
    output: strict="combined.strictFilter.vcf",
            snpeff="combined.strictFilter.snpeff.vcf",
            strictflagged=temp("combined_strictFlagged.vcf"),
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],effgenome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="rnafilter"
    shell: "{params.gatk} -T VariantFiltration -R {params.genome} --variant {input} -window 35 -cluster 3 --filterExpression \"QD < 2.0 || FS > 30.0\" --filterName \"my_filter\" -o {output.strictflagged}; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.strictflagged} --excludeFiltered -o {output.strict}; {params.snpeff} -v -canon -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.effgenome} {output.strict} > {output.snpeff}; mkdir -p sample_vcfs"