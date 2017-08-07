rule strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets=ancient("exome_targets.bed")
    output: vcf=config['project']['workpath']+"/strelka_out/{x}.vcf",
            vcffiltered=config['project']['workpath']+"/strelka_out/{x}_filtered.vcf",
            outdir=config['project']['workpath']+"/strelka_out/{x}",
            csvstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.csv",
            htmlstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.html",
            out=config['project']['workpath']+"/strelka_out/{x}.snpeff.out"
    params: dir=config['project']['workpath'],strelkaconfig=config['references'][pfamily]['STRELKA_CONFIG'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets="exome_targets.bed",snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka/2.7.1; module load python; configureStrelkaSomaticWorkflow.py --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --runDir={output.outdir} --exome; cd {output.outdir}; ./runWorkflow.py -m local -j {threads}; {params.gatk} -T CombineVariants -R {params.genome} --variant results/variants/somatic.snvs.vcf.gz --variant results/variants/somatic.indels.vcf.gz -L {params.dir}/{params.targets} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcffiltered}; cd {params.dir}/strelka_out; module load snpEff; {params.snpeff} -v {params.snpeffgenome} -interval {params.dir}/{params.targets} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples {params.dir}/pairs {output.vcffiltered} > {output.out}"