rule wgs_strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
    output: vcf=config['project']['workpath']+"/strelka_out/{x}.vcf",
            outdir=config['project']['workpath']+"/strelka_out/{x}",
            csvstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.csv",
            htmlstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.html",
            out=config['project']['workpath']+"/strelka_out/{x}.snpeff.out"
    params: dir=config['project']['workpath'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka/2.7.1; module load python; configureStrelkaSomaticWorkflow.py --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --runDir={output.outdir}; cd {output.outdir}; ./runWorkflow.py -m local -j {threads}; {params.gatk} -T CombineVariants -R {params.genome} --variant results/passed.somatic.snvs.vcf --variant results/passed.somatic.indels.vcf --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}; cd {params.dir}/strelka_out; module load snpEff; {params.snpeff} -v {params.snpeffgenome} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples {params.dir}/pairs {output.vcf} > {output.out}"