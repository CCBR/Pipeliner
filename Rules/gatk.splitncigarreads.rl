rule gatk_splitcigar:
      input:  "{x}.dedup_unsplit.bam"
      output: "{x}.dedup.bam"
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="pl:split"
      shell:  "module load GATK/3.6; GATK -m 64g SplitNCigarReads -I {input} -R {params.genome} -o {output} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"