rule germline_filter_rnaseqvar:
    input: "combined.vcf"
    output: strict="targets.strictFilter.vcf",
            strictflagged=temp("combined_strictFlagged.vcf"),
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="rnafilter"
    shell: "{params.gatk} -T VariantFiltration -R {params.genome} --variant {input} -window 35 -cluster 3 --filterExpression \"QD < 2.0 || FS > 60.0 || MQ < 40.0\" --filterName \"my_filter\" -o {output.strictflagged}; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.strictflagged} --excludeFiltered -o {output.strict}"