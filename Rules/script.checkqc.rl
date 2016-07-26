rule script_checkqc:
      input:  "{x}.recal.bam"
      output: cov= "{x}.uniqcov",
              covall = "{x}_covall.txt",
              stat = "{x}_covstat.txt",
              ism = "{x}.insert_size_metrics",
              hist = "{x}_PIChist.pdf",
              flag = "{x}.flagstat.txt",
              tbl="{x}.recal.tbl"
      params: covcalc = config['bin'][pfamily]['COVCALC'],
              covfreq = config['bin'][pfamily]['COVFREQ'],
              gatk=config['bin'][pfamily]['GATK'],
              genome=config['references'][pfamily]['GENOME'],
              picard2=config['bin'][pfamily]['PICARD2'],
              pichist=config['bin'][pfamily]['PICHIST'],
              readdist = config['bin'][pfamily]['READDIST'],
              refflat=config['references'][pfamily]['REFFLAT'],rname="pl:checkqc"
      shell:  "coverageBed -abam {input} -hist -b {params.refflat} > {output.cov}; \
              perl {params.covcalc} {output.cov} {output.stat}; \
              grep 'all' {output.cov} > {output.covall}; \
              echo 'source(\"{params.covfreq}\"); RunData(file2 = \"{output.covall}\")' | R --vanilla; \
              {params.picard2} INPUT={input} REFERENCE_SEQUENCE={params.genome} OUTPUT={output.ism} HISTOGRAM_FILE={output.hist}; \
              R --vanilla < {params.pichist} --args {output.ism} {output.hist} {input}; \
              samtools flagstat {input} > {output.flag}; \
              {params.gatk} -T ReadLengthDistribution -I {input} -R {params.genome} -o {output.tbl}; \
              echo 'source(\"{params.readdist}\"); RunHist(file = \"{output.tbl}\")' | R --vanilla;"
