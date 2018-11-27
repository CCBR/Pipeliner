rule strelka:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets=ancient("exome_targets.bed"),
           normalbai=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam.bai",
           tumorbai=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam.bai"
    output: vcf=temp(config['project']['workpath']+"/strelka_out/{x}.vcf"),
            vcffiltered=temp(config['project']['workpath']+"/strelka_out/{x}_filtered.vcf"),
            final=config['project']['workpath']+"/strelka_out/{x}_FINAL.vcf",
            outdir=config['project']['workpath']+"/strelka_out/{x}",
            csvstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.csv",
            htmlstats=config['project']['workpath']+"/strelka_out/{x}.strelka.stats.html",
            out=config['project']['workpath']+"/strelka_out/{x}.snpeff.out",
            names=temp(config['project']['workpath']+"/strelka_out/{x}.samples"),
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],dir=config['project']['workpath'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],targets="exome_targets.bed",snpeff=config['bin'][pfamily]['SNPEFF'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:strelka_calls"
    threads: 4
    shell: "module load strelka/2.9.0; module load python; configureStrelkaSomaticWorkflow.py --ref={params.genome} --tumor={params.dir}/{input.tumor} --normal={params.dir}/{input.normal} --runDir={output.outdir} --exome; cd {output.outdir}; ./runWorkflow.py -m local -j {threads}; {params.gatk} -T CombineVariants -R {params.genome} --variant results/variants/somatic.snvs.vcf.gz --variant results/variants/somatic.indels.vcf.gz -L {params.dir}/{params.targets} --assumeIdenticalSamples --filteredrecordsmergetype KEEP_UNCONDITIONAL -o {output.vcf}; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcffiltered}; cd {params.dir}/strelka_out; module load samtools/1.6; echo $'NORMAL {params.normalsample}\nTUMOR {params.tumorsample}' > {output.names}; bcftools reheader -o {output.final} -s {output.names} {output.vcffiltered}; module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v {params.snpeffgenome} -c {params.effconfig} -interval {params.dir}/{params.targets} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples {params.dir}/pairs {output.final} > {output.out}"