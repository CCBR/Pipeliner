rule snpeff:
    input: vcf="combined.vcf",
           strict="exome.strictFilter.vcf",
           relax="exome.relaxedFilter.vcf",
           targets="exome_targets.bed"
    output: vcf="exome.snpeff.vcf",
            strict="exome.strictFilter.snpeff.vcf",
            relax="exome.relaxedFilter.snpeff.vcf"
    params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],exons=config['references'][pfamily]['EXONS'],rname="pl:snpeff"
    shell: """module load snpEff/4.3t
           java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff.stats.csv -stats snpeff.stats.html {params.genome} {input.vcf} > {output.vcf}
           java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff_strictFilter.stats.csv -stats snpeff_strictFilter.stats.html {params.genome} {input.strict} > {output.strict}
           java -Xmx12g -jar $SNPEFF_JAR -v -canon -c {params.effconfig} -csvStats snpeff_relaxFilter.stats.csv -stats snpeff_relaxFilter.stats.html {params.genome} {input.relax} > {output.relax}
           module load samtools
           bcftools stats --exons {params.exons} --fasta-ref {params.ref} --regions-file {input.targets} {input.vcf} > raw_variants.bcftools
           bcftools stats --exons {params.exons} --fasta-ref {params.ref} --regions-file {input.targets} {input.strict} > strict_variants.bcftools
           bcftools stats --exons {params.exons} --fasta-ref {params.ref} --regions-file {input.targets} {input.relax} > relax_variants.bcftools"""