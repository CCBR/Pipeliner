rule germline_filter_wgs:
    input: "combined.vcf"
    output: strict="combined.strictFilter.vcf",
            relax="combined.relaxedFilter.vcf",
            snps=temp("snps.vcf"),
            indels=temp("indels.vcf"),
            strictsnps=temp("snps_strictFlagged.vcf"),
            strictindels=temp("indels_strictFlagged.vcf"),
            strictflagged=temp("combined_strictFlagged.vcf"),
            relaxsnps=temp("snps_relaxFlagged.vcf"),
            relaxindels=temp("indels_relaxFlagged.vcf"),
            relaxflagged=temp("combined_relaxFlagged.vcf")
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],batch ="-l nodes=1:gpfs -q ccr",rname="filter"
    shell: "{params.gatk} -T SelectVariants -R {params.genome} -V {input} -selectType SNP -o {output.snps}; {params.gatk} -T SelectVariants -R {params.genome} --variant {input} -selectType INDEL -o {output.indels}; {params.gatk} -T VariantFiltration -R {params.genome} --variant {output.snps} --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o {output.strictsnps}; {params.gatk} -T VariantFiltration -R {params.genome} --variant {output.indels} --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o {output.strictindels}; {params.gatk} -T CombineVariants -R {params.genome} --variant {output.strictsnps} --variant {output.strictindels} -o {output.strictflagged} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNSORTED; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.strictflagged} --excludeFiltered -o {output.strict}; {params.gatk} -T VariantFiltration -R {params.genome} --variant {output.snps} --filterExpression "QD < 2.0 || FS > 60.0" --filterName "my_snp_filter" -o {output.relaxsnps}; {params.gatk} -T VariantFiltration -R {params.genome} --variant {output.indels} --filterExpression "QD < 2.0 || FS > 200.0" --filterName "my_indel_filter" -o {output.relaxindels}; {params.gatk} -T CombineVariants -R {params.genome} --variant {output.relaxsnps} --variant {output.relaxindels} -o {output.relaxflagged} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNSORTED; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.relaxflagged} --excludeFiltered -o {output.relax}"