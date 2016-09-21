rule strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf="strelka_out/{x}.vcf",
            outdir="strelka_out/{x}",
            csvstats="strelka_out/{x}.strelka.stats.csv",
            htmlstats="strelka_out/{x}.strelka.stats.html",
#            snpeffvcf="strelka_out/{x}.snpeff.vcf"
    params: dir=config['project']['workpath'],strelkaconfig=config['references'][pfamily]['STRELKA_CONFIG'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets="exome_targets.bed",snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka; $STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl --config={params.strelkaconfig} --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --output-dir={params.dir}/{output.outdir}; cd {output.outdir}; make -j 4; {params.gatk} -T CombineVariants -R {params.genome} --variant results/passed.somatic.snvs.vcf --variant results/passed.somatic.indels.vcf -L {params.dir}/{params.targets} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {params.dir}/{output.vcf}; cd {params.dir}/strelka_out; module load snpEff; {params.snpeff} -v {params.snpeffgenome} -interval {params.dir}/{params.targets} -cancer -csvStats {params.dir}/{output.csvstats} -stats {params.dir}/{output.htmlstats} -cancerSamples pairs {params.dir}/{output.vcf}"