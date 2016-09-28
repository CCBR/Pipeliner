rule gatk_recal:
      input:  "{x}.realign.bam"
      output: bam="{x}.recal.bam",
              re=temp("{x}_recal_data.bam.grp")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knownrecal=config['references'][pfamily]['KNOWNRECAL']
      threads: 32
      shell:  "{params.gatk} -T BaseRecalibrator -I {input} -R {params.genome} {params.knownrecal} -nct {threads} -o {output.re}; {params.gatk} -T PrintReads -R {params.genome} -I {input} -BQSR {output.re} -o {output.bam}"