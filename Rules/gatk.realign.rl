rule gatk_realign:
        input:  "{x}.dedup.bam"
        output: re="{x}.realign.bam",
                int="{x}.fin.bam.intervals"
        params: gatk=config['bin'][pfamily]['GATK'],
                genome=config['references'][pfamily]['GENOME'],
                sam=config['bin'][pfamily]['SAMTOOLS'],
                picard3=config['bin'][pfamily]['PICARD3'],
                indelsites=config['references'][pfamily]['INDELSITES'],rname="pl:realign"
        shell:  "{params.sam} index {input};{params.gatk} -T RealignerTargetCreator -I {input} -R {params.genome} -known {params.indelsites} -o {output.int}; {params.gatk} -T IndelRealigner -R {params.genome} -I {input} -targetIntervals {output.int} -o {output.re}"
#        shell:  "{params.sam} index {input}; {params.sam} faidx {params.genome};  {params.gatk} -T RealignerTargetCreator -I {input} -R {params.genome} -known {params.indelsites} -o {output.int}; {params.gatk} -T IndelRealigner -R {params.genome} -I {input} -targetIntervals {output.int} -o {output.re}"
