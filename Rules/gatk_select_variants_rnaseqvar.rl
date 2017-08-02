rule gatk_select_variants_rnaseqvar:
    input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam",
           vcf="targets.strictFilter.vcf"
    output: vcf="sample_vcfs/{x}.sample.vcf",
            csvstats="sample_vcfs/{x}.stats.csv",
            htmlstats="sample_vcfs/{x}.stats.html",
            bed="sample_vcfs/{x}.snpeff.bed"
    params: sample=lambda wildcards: config['project']['units'][wildcards.x],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],effgenome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:extract"
    shell: "{params.gatk} -T SelectVariants -R {params.genome} -V {input.vcf} -sn {params.sample} -env -o {output.vcf}; {params.snpeff} -v -c {params.effconfig} -o bed -csvStats {output.csvstats} -stats {output.htmlstats} {params.effgenome} {output.vcf} > {output.bed}"