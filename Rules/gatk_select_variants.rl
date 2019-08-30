rule gatk_select_variants:
    input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam",
           vcf="exome.relaxedFilter.vcf",
    output: vcfgz="sample_vcfs/{x}.sample.vcf.gz",
            csvstats="sample_vcfs/{x}.stats.csv",
            htmlstats="sample_vcfs/{x}.stats.html",
            bed="sample_vcfs/{x}.snpeff.bed",
            stats="sample_vcfs/{x}.bcftools"
    params: sample=lambda wildcards: config['project']['units'][wildcards.x],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets="exome_targets.bed",snpeff=config['bin'][pfamily]['SNPEFF'],effgenome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],exons=config['references'][pfamily]['EXONS'],rname="pl:extract"
    shell: """{params.gatk} -T SelectVariants -R {params.genome} -V {input.vcf} -sn {params.sample} -env -o sample_vcfs/{params.sample}.sample.vcf
           module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v -c {params.effconfig} -o bed -csvStats {output.csvstats} -stats {output.htmlstats} {params.effgenome} sample_vcfs/{params.sample}.sample.vcf > {output.bed}
           module load samtools
           bgzip sample_vcfs/{params.sample}.sample.vcf
           tabix -p vcf {output.vcfgz}
           bcftools stats --exons {params.exons} --fasta-ref {params.genome} --regions-file {params.targets} {output.vcfgz} > {output.stats}"""