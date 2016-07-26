rule strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".realign.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".realign.bam"
    output: vcf="strelka/{x}.vcf",
            dir="strelka/{x}"
    params: dir=config['project']['workpath'],strelkaconfig=config['references']['STRELKA_CONFIG'],gatk=config['bin']['GATK'],genome=config['references']['GENOME'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka; $STRELKA_INSTALL_DIR/bin/configureStrelkaWorkflow.pl --config={params.strelkaconfig} --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --output-dir={params.dir}/{output.dir}; cd {output.dir}; make -j 4; {params.gatk} -T CombineVariants -R {params.genome} --variant results/passed.somatic.snvs.vcf --variant results/passed.somatic.indels.vcf --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {params.dir}/{output.vcf}"