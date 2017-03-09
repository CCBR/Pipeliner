rule gatk_splitcigar:
      input:  "{x}.dedup_unsplit.bam"
      output: "{x}.dedup.bam"
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="pl:split"
      shell:  "{params.gatk} -T SplitNCigarReads -I {input} -R {params.genome} -o {output} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"