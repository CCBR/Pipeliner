rule gatk_realign:
        input:  "{x}.dedup.bam"
        output: re=temp("{x}.realign.bam"),
                int="{x}.fin.bam.intervals"
        params: gatk=config['bin'][pfamily]['GATK'],
                genome=config['references'][pfamily]['GENOME'],
                sam=config['bin'][pfamily]['SAMTOOLS'],
                picard3=config['bin'][pfamily]['PICARD3'],
                knownindels=config['references'][pfamily]['KNOWNINDELS'],rname="pl:realign"
        shell:  "{params.sam} index {input};{params.gatk} -T RealignerTargetCreator -I {input} -R {params.genome} {params.knownindels} -o {output.int}; {params.gatk} -T IndelRealigner -R {params.genome} -I {input} -targetIntervals {output.int} {params.knownindels} -o {output.re}"