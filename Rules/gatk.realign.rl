rule gatk_realign:
        input:  "{x}.dedup.bam"
        output: re=temp("{x}.realign.bam"),
                int="{x}.fin.bam.intervals"
        params: gatk=config['bin'][pfamily]['GATK'],
                genome=config['references'][pfamily]['GENOME'],
                sam=config['bin'][pfamily]['SAMTOOLS'],
                picard3=config['bin'][pfamily]['PICARD3'],
                knownindels=config['references'][pfamily]['KNOWNINDELS'],rname="pl:realign"
        threads: 2
        shell:  "{params.sam} index {input};module load GATK/3.5-0; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T RealignerTargetCreator -I {input} -R {params.genome} {params.knownindels} -o {output.int} -nt 2; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T IndelRealigner -R {params.genome} -I {input} -targetIntervals {output.int} {params.knownindels} -o {output.re}"