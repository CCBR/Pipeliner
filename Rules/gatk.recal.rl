rule gatk_recal:
      input:  "{x}.realign.bam"
      output: bam="{x}.recal.bam",
              re=temp("{x}_recal_data.bam.grp")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNRECAL'],rname="pl:recal"
      threads: 16
      shell:  "module load GATK/3.8-0; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T BaseRecalibrator -I {input} -R {params.genome} {params.knowns} -nct {threads} -o {output.re} --disable_auto_index_creation_and_locking_when_reading_rods --use_jdk_inflater --use_jdk_deflater; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T PrintReads -R {params.genome} -nct 8 -I {input} --use_jdk_inflater --use_jdk_deflater -BQSR {output.re} -o {output.bam}"
