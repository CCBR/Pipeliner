rule strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".realign.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".realign.bam"
    output: vcf="strelka_out/{x}.vcf",
            outdir="strelka_out/{x}"
    params: dir=config['project']['workpath'],strelkaconfig=config['references'][pfamily]['STRELKA_CONFIG'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets=config['references'][pfamily]['REFFLAT'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka; $STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl --config={params.strelkaconfig} --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --output-dir={params.dir}/{output.outdir}; cd {output.outdir}; make -j 4; {params.gatk} -T CombineVariants -R {params.genome} --variant results/passed.somatic.snvs.vcf --variant results/passed.somatic.indels.vcf -L {params.targets} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {params.dir}/{output.vcf}"